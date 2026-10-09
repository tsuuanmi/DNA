//! Projection of compact sample scientific evidence into the public JSON contract.

use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

use crate::error::{ReportError, Result};
use crate::model::alignment::Orientation;
use crate::model::read_observation::{SangerAttachment, SangerRejection};
use crate::model::reference::Reference;
use crate::model::result::{AlignmentResult, IntervalResult, PeakHeightsResult, ReferenceResult};
use crate::model::sample_evidence::SampleEvidence;
use crate::model::sample_result::{
    ReadRejectionResult, RejectedSampleReadResult, SampleCoverageResult,
    SampleEvidenceProfileResult, SampleEvidenceResult, SampleLocusDifferenceObservationResult,
    SampleLocusDifferenceResult, SampleLocusSupportTopologyResult, SampleOverlapResult,
    SampleProvenanceResult, SampleReadResult, SampleVariantCallResult,
    SampleVariantOppositionResult, SampleVariantResult, SampleVariantSupportResult,
    SampleVariantSupportTopologyResult,
};
use crate::plugin::PluginDescriptor;
use crate::profile::ProfileIdentity;
use crate::report::json::{project_plugins, project_profile};
use crate::report::notation::{self, SampleNotation};
use crate::report::sanger_call;

/// Inputs consumed to build one immutable sample-evidence document.
pub(crate) struct CompletedSampleEvidence {
    pub(crate) sample_id: String,
    pub(crate) reference: Reference,
    pub(crate) profile: ProfileIdentity,
    pub(crate) evidence: SampleEvidence,
    /// Per-read represented calls, present only when the profile declares notation.
    pub(crate) notation: Option<SampleNotation>,
    /// Plugins of the workflow, in execution order.
    pub(crate) plugins: &'static [&'static PluginDescriptor],
    /// Sanger evidence joined to the core sample evidence.
    pub(crate) sanger: SangerSampleEvidence,
}

/// Sanger evidence of a sample's reads, keyed by read content identity. The
/// core sample evidence never carries it (ADR-0069).
pub(crate) struct SangerSampleEvidence {
    pub(crate) reads: BTreeMap<String, SangerAttachment>,
    pub(crate) rejected: BTreeMap<String, SangerRejection>,
}

/// One admitted read's Sanger evidence and selected orientation, in
/// read-registry order.
#[derive(Clone, Copy)]
struct JoinedRead<'a> {
    sanger: &'a SangerAttachment,
    orientation: Orientation,
}

/// Builds `dna.sample_evidence/v10` without filesystem side effects.
pub(crate) fn build(completed: CompletedSampleEvidence) -> Result<SampleEvidenceResult> {
    let CompletedSampleEvidence {
        sample_id,
        reference,
        profile,
        evidence,
        notation,
        plugins,
        sanger,
    } = completed;
    if evidence.reference_sha256 != reference.sequence_sha256 {
        return Err(ReportError::Inconsistent(
            "sample evidence reference identity does not match report reference",
        )
        .into());
    }

    let (read_names, rejected_names) = reviewer_read_names(&evidence)?;
    let rejected_reads = evidence
        .rejected_reads
        .iter()
        .zip(rejected_names)
        .map(|(read, name)| {
            let rejected =
                sanger
                    .rejected
                    .get(&read.input_sha256)
                    .ok_or(ReportError::Inconsistent(
                        "rejected sample read lacks its Sanger evidence",
                    ))?;
            Ok(RejectedSampleReadResult {
                name,
                sha256: read.input_sha256.clone(),
                integrity: crate::report::signal::project_integrity(&rejected.integrity),
                callability: crate::report::callability::project(&rejected.callability),
                rejection: ReadRejectionResult {
                    reason: "callable_calls_below_minimum",
                    callable_calls: rejected.rejection.callable_calls,
                    minimum_callable_calls: rejected.rejection.minimum_callable_calls,
                },
            })
        })
        .collect::<Result<_>>()?;
    let joined =
        evidence
            .reads
            .iter()
            .map(|read| {
                Ok(JoinedRead {
                    sanger: sanger.reads.get(&read.input_sha256).ok_or(
                        ReportError::Inconsistent("sample read lacks its Sanger evidence"),
                    )?,
                    orientation: read.alignment.orientation,
                })
            })
            .collect::<Result<Vec<_>>>()?;
    let notation = notation
        .map(|notation| {
            let identities = evidence
                .reads
                .iter()
                .map(|read| read.input_sha256.as_str())
                .collect::<Vec<_>>();
            notation::project(&reference, &identities, &read_names, notation)
        })
        .transpose()?;
    let reads = evidence
        .reads
        .into_iter()
        .zip(read_names.iter().zip(&joined))
        .map(|(read, (name, joined))| SampleReadResult {
            name: name.clone(),
            sha256: read.input_sha256,
            integrity: crate::report::signal::project_integrity(&joined.sanger.signal.integrity),
            callability: crate::report::callability::project(&joined.sanger.callability),
            alignment: AlignmentResult {
                orientation: read.alignment.orientation,
                callable_bases: read.alignment.callable_bases,
                identity: read.alignment.identity,
                gap_opens: read.alignment.gap_opens,
                unresolved_bases: read.alignment.unresolved_bases,
                masked_bases: read.alignment.masked_bases,
                reference_segments: read
                    .alignment
                    .reference_segments
                    .into_iter()
                    .map(|segment| IntervalResult {
                        start: segment.start_0based,
                        end: segment.end_0based_exclusive,
                    })
                    .collect(),
                callable_reference_segments: read
                    .alignment
                    .callable_reference_segments
                    .into_iter()
                    .map(|segment| IntervalResult {
                        start: segment.start_0based,
                        end: segment.end_0based_exclusive,
                    })
                    .collect(),
                wraps_origin: read.alignment.wraps_origin,
            },
        })
        .collect();

    let coverage = evidence
        .coverage
        .into_iter()
        .map(|segment| SampleCoverageResult {
            reference: IntervalResult {
                start: segment.start_0based,
                end: segment.end_0based_exclusive,
            },
            read_depth: segment.read_depth,
            forward_depth: segment.forward_depth,
            reverse_depth: segment.reverse_depth,
        })
        .collect();

    let overlaps = evidence
        .overlaps
        .into_iter()
        .map(|overlap| {
            Ok(SampleOverlapResult {
                left: read_name(&read_names, overlap.left_read_index)?.into(),
                right: read_name(&read_names, overlap.right_read_index)?.into(),
                shared_positions: overlap.shared_positions,
                comparable_bases: overlap.comparable_bases,
                agreements: overlap.agreements,
                conflicts: overlap.conflicts,
                agreement: overlap.agreement,
                eligible: overlap.eligible,
                exclusion_reasons: overlap.exclusion_reasons,
            })
        })
        .collect::<Result<Vec<_>>>()?;

    let locus_differences = evidence
        .locus_differences
        .into_iter()
        .map(|difference| {
            let observations = difference
                .observations
                .into_iter()
                .map(|observation| {
                    let sanger = observation
                        .call_index_0based
                        .map(|call| -> Result<(u8, bool)> {
                            let read = joined_read(&joined, observation.read_index)?;
                            Ok((
                                call_quality(read.sanger, call)?,
                                read.sanger.in_noisy_region(call),
                            ))
                        })
                        .transpose()?;
                    Ok(SampleLocusDifferenceObservationResult {
                        read: read_name(&read_names, observation.read_index)?.into(),
                        state: observation.state,
                        base: observation.base,
                        quality: sanger.map(|(quality, _)| quality),
                        profile: observation
                            .profile
                            .map(|profile| SampleEvidenceProfileResult {
                                a: profile.weights[0],
                                c: profile.weights[1],
                                g: profile.weights[2],
                                t: profile.weights[3],
                            }),
                        in_noisy_region: sanger.map(|(_, noisy)| noisy),
                    })
                })
                .collect::<Result<Vec<_>>>()?;
            Ok(SampleLocusDifferenceResult {
                position: difference.position_1based,
                reference: difference.reference_base,
                support_topology: SampleLocusSupportTopologyResult {
                    reads: difference.support_topology.reads,
                    forward_reads: difference.support_topology.forward_reads,
                    reverse_reads: difference.support_topology.reverse_reads,
                    reference_reads: difference.support_topology.reference_reads,
                    alternate_reads: difference.support_topology.alternate_reads,
                    unresolved_reads: difference.support_topology.unresolved_reads,
                    deletion_reads: difference.support_topology.deletion_reads,
                    masked_reads: difference.support_topology.masked_reads,
                    profile_reads: difference.support_topology.profile_reads,
                    profile_forward_reads: difference.support_topology.profile_forward_reads,
                    profile_reverse_reads: difference.support_topology.profile_reverse_reads,
                },
                observations,
            })
        })
        .collect::<Result<Vec<_>>>()?;

    let variants = evidence
        .variants
        .into_iter()
        .map(|variant| {
            let support = variant
                .support
                .into_iter()
                .map(|support| {
                    let read = joined_read(&joined, support.read_index)?;
                    Ok(SampleVariantSupportResult {
                        read: read_name(&read_names, support.read_index)?.into(),
                        eligible: support.eligible,
                        exclusion_reasons: support.exclusion_reasons,
                        calls: support
                            .calls
                            .into_iter()
                            .map(|call| {
                                let sanger = sanger_call::evidence(
                                    &read.sanger.calls,
                                    &read.sanger.quality,
                                    read.orientation,
                                    call.call_index_0based,
                                )
                                .map_err(ReportError::CallEvidence)?;
                                Ok(SampleVariantCallResult {
                                    role: call.role,
                                    base: call.base,
                                    peaks: PeakHeightsResult::from(sanger.peak_heights),
                                    quality: sanger.quality,
                                })
                            })
                            .collect::<Result<_>>()?,
                    })
                })
                .collect::<Result<Vec<_>>>()?;
            Ok(SampleVariantResult {
                position: variant.position_1based,
                reference: variant.reference,
                alternate: variant.alternate,
                kind: variant.kind,
                support_topology: SampleVariantSupportTopologyResult {
                    reads: variant.support_topology.reads,
                    eligible_reads: variant.support_topology.eligible_reads,
                    forward_reads: variant.support_topology.forward_reads,
                    reverse_reads: variant.support_topology.reverse_reads,
                    eligible_forward_reads: variant.support_topology.eligible_forward_reads,
                    eligible_reverse_reads: variant.support_topology.eligible_reverse_reads,
                },
                support,
                opposition: SampleVariantOppositionResult {
                    reads: variant
                        .opposition
                        .read_indices
                        .iter()
                        .map(|&index| read_name(&read_names, index).map(str::to_owned))
                        .collect::<Result<Vec<_>>>()?,
                    forward_reads: variant.opposition.forward_reads,
                    reverse_reads: variant.opposition.reverse_reads,
                },
            })
        })
        .collect::<Result<Vec<_>>>()?;

    Ok(SampleEvidenceResult {
        schema_version: "dna.sample_evidence/v10",
        sample_id,
        provenance: SampleProvenanceResult {
            reference: ReferenceResult {
                name: reference.name,
                topology: reference.topology,
                sha256: reference.sequence_sha256,
            },
            configuration_sha256: evidence.configuration_sha256,
            profile: project_profile(profile),
            plugins: project_plugins(plugins),
        },
        reads,
        rejected_reads,
        coverage,
        overlaps,
        locus_differences,
        variants,
        notation,
    })
}

/// Reviewer-facing names of the admitted and the rejected reads, unique
/// across both.
fn reviewer_read_names(evidence: &SampleEvidence) -> Result<(Vec<String>, Vec<String>)> {
    let mut unique = BTreeSet::new();
    let mut name_of = |input_name: &str| -> Result<String> {
        let name = Path::new(input_name)
            .file_stem()
            .and_then(|value| value.to_str())
            .filter(|value| !value.is_empty())
            .ok_or(ReportError::InvalidReadName)?
            .to_owned();
        if !unique.insert(name.clone()) {
            return Err(ReportError::DuplicateReadName { name }.into());
        }
        Ok(name)
    };
    let names = evidence
        .reads
        .iter()
        .map(|read| name_of(&read.input_name))
        .collect::<Result<Vec<_>>>()?;
    let rejected = evidence
        .rejected_reads
        .iter()
        .map(|read| name_of(&read.input_name))
        .collect::<Result<Vec<_>>>()?;
    Ok((names, rejected))
}

fn joined_read<'a>(joined: &[JoinedRead<'a>], index: usize) -> Result<JoinedRead<'a>> {
    Ok(*joined
        .get(index)
        .ok_or(ReportError::MissingRead { index })?)
}

/// Relative quality of call `index`, whose record must carry that index.
fn call_quality(sanger: &SangerAttachment, index: usize) -> Result<u8> {
    Ok(sanger
        .quality
        .per_call
        .get(index)
        .filter(|quality| quality.index_0based == index)
        .ok_or(ReportError::Inconsistent(
            "sample locus call lacks matching quality evidence",
        ))?
        .relative_quality_score)
}

fn read_name(read_names: &[String], index: usize) -> Result<&str> {
    Ok(read_names
        .get(index)
        .map(String::as_str)
        .ok_or(ReportError::MissingRead { index })?)
}
