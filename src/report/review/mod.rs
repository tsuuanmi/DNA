//! The review list of one sample, `dna.review/v1` (PROP-0004).
//!
//! Items and labels are computed from the sample analysis and its call alone,
//! never from a truth: the consensus sites, the reads' observations and the
//! core's trust in each call, and the Sanger scan position and relative
//! quality of each call, joined by read identity and call index. The rules'
//! version is the `review` plugin's version.

mod items;
mod labels;
mod lengths;
mod ranges;
mod reads;
mod second;
mod sites;
#[cfg(test)]
mod test_support;
mod zones;

use self::labels::summarize;
use self::ranges::{clip, inside, merge, segments, span};
use self::reads::{masked_spans, unreliable_spans};
use self::second::own;
use self::zones::zones;

use std::collections::{BTreeMap, BTreeSet};

use crate::config::ReviewConfig;
use crate::error::{ReportError, Result};
use crate::model::consensus_result::StretchResult;
use crate::model::result::{IntervalResult, ReferenceResult};
use crate::model::review_result::{
    ItemKind, LabelResult, MixtureResult, ReviewCallResult, ReviewProvenanceResult,
    ReviewReadResult, ReviewResult, ReviewStatus, SignalProvenanceResult,
};
use crate::model::sample_result::ConventionResult;
use crate::model::sanger_read::{Rejection, SangerRead};
use crate::report::json::{project_plugins, project_profile};
use crate::report::notation::project_empop;
use crate::report::sample::{reviewer_names, reviewer_read_names};
use crate::sample_analysis::windows::{LengthWindow, windows};
use crate::sample_analysis::{
    CallTrust, ConsensusCall, NotationCall, ReviewCheck, ReviewList, SampleAnalysis, TractNote,
    TractNoteRun,
};
use dna_align::model::alignment::{AlignmentColumn, Orientation};
use dna_call::model::called_read::CalledRead;
use dna_call::read_call::CallConfig;
use dna_consensus::ConsensusConfig;
use dna_consensus::decided::not_analysed;
use dna_consensus::model::consensus::ConsensusSite;
use dna_kernel::model::reference::Reference;
use dna_kernel::plugin::{Contract, PluginDescriptor, PluginFamily};
use dna_kernel::profile::{Convention, LengthTract, Profile};
use dna_kernel::read_evidence::{EvidenceReason, ReadEvidence};
use dna_normalize::length_convention::ShortestDominant;
use dna_sample::model::sample_evidence::LocusObservations;
use dna_sanger::run_signal::RunSignal;
use dna_sanger::signal_processing::SignalProcessingConfig;

/// The review list: items and labels over the consensus, the call, and the
/// reads (PROP-0004). Its version is the labelling rules' version.
pub(crate) const PLUGIN: PluginDescriptor = PluginDescriptor {
    id: "review",
    family: PluginFamily::PostCalling,
    version: 33,
    provides: &[Contract::ReviewList],
    requires: &[
        Contract::ReadEvidence,
        Contract::Consensus,
        Contract::CalledVariants,
    ],
    config_sections: &["review"],
};

/// Inputs consumed to build one review document.
pub(crate) struct CompletedReview {
    pub(crate) sample_id: String,
    pub(crate) reference: Reference,
    pub(crate) profile: Profile,
    pub(crate) configuration_sha256: String,
    pub(crate) outcome: Reviewed,
    /// Plugins of the workflow, in execution order.
    pub(crate) plugins: &'static [&'static PluginDescriptor],
    /// The run signal the calls are judged against, with its file's SHA-256.
    pub(crate) signal: Option<(RunSignal, String)>,
    /// Each admitted read's phase mixtures, by content SHA-256 (ADR-0082).
    pub(crate) mixtures: BTreeMap<String, Vec<MixtureResult>>,
}

/// What the sample analysis and its call gave.
pub(crate) enum Reviewed {
    /// Every read was rejected; each with its cause, by content SHA-256.
    Rejected(BTreeMap<String, Rejection>),
    /// At least one read was admitted; the call, or why it failed.
    Analyzed {
        analysis: Box<SampleAnalysis>,
        call: Result<ConsensusCall>,
    },
}

/// Builds `dna.review/v1` without filesystem side effects.
pub(crate) fn build(completed: CompletedReview) -> Result<ReviewResult> {
    let CompletedReview {
        sample_id,
        reference,
        profile,
        configuration_sha256,
        outcome,
        plugins,
        signal,
        mixtures,
    } = completed;
    let regions = clip(&profile.regions, reference.len());
    let provenance = ReviewProvenanceResult {
        reference: ReferenceResult {
            name: reference.name.clone(),
            topology: reference.topology,
            sha256: reference.sequence_sha256.clone(),
        },
        configuration_sha256,
        profile: project_profile(profile.identity().clone()),
        plugins: project_plugins(plugins),
        rules: PLUGIN.version,
        signal: signal
            .as_ref()
            .map(|(signal, sha256)| SignalProvenanceResult {
                run_id: signal.run_id.clone(),
                sha256: sha256.clone(),
                selection: signal.selection,
            }),
        run_context: match &outcome {
            Reviewed::Analyzed { analysis, .. } => analysis
                .run_context()
                .map(crate::model::run_context_result::ContextUseResult::of),
            Reviewed::Rejected(_) => None,
        },
    };
    let (analysis, call) = match outcome {
        Reviewed::Rejected(rejected) => {
            let names = reviewer_names(
                rejected
                    .values()
                    .map(|read| (read.input_name.as_str(), read.input_sha256.as_str())),
            )?;
            let reads = rejected
                .values()
                .zip(names)
                .map(|(read, name)| ReviewReadResult {
                    name,
                    sha256: read.input_sha256.clone(),
                    strand: None,
                    admitted: false,
                    rejection: Some(read.cause),
                    masked: None,
                    unreliable: None,
                    mixtures: None,
                })
                .collect();
            return Ok(ReviewResult {
                schema_version: "dna.review/v1",
                sample_id,
                status: ReviewStatus::Rejected,
                error: None,
                provenance,
                reads,
                calls: None,
                callable: None,
                reference_only: None,
                not_analysed: None,
                empop: None,
                zones: None,
                summary: summarize(&[], ReviewStatus::Rejected, &regions, &[], &[], 0),
                items: Vec::new(),
            });
        }
        Reviewed::Analyzed { analysis, call } => (analysis, call),
    };

    let reference_only = analysis
        .consensus
        .reference_only
        .iter()
        .map(|position| position + 1)
        .filter(|&position| inside(&regions, position))
        .collect::<Vec<_>>();
    let (names, rejected_names) = reviewer_read_names(&analysis.evidence)?;
    let mut reads = analysis
        .evidence
        .reads
        .iter()
        .zip(&names)
        .map(|(read, name)| {
            let missing = || ReportError::Inconsistent("review read lacks its evidence");
            let called = analysis
                .called
                .iter()
                .find(|called| called.input_sha256 == read.input_sha256)
                .ok_or_else(missing)?;
            let sanger = analysis
                .sanger
                .get(&read.input_sha256)
                .ok_or_else(missing)?;
            Ok(ReviewReadResult {
                name: name.clone(),
                sha256: read.input_sha256.clone(),
                strand: Some(read.alignment.orientation),
                admitted: true,
                rejection: None,
                masked: Some(masked_spans(
                    called.alignment.columns(),
                    &sanger.callability.segments,
                    &analysis.trust(&read.input_sha256).ok_or_else(missing)?,
                    &regions,
                )),
                unreliable: Some(unreliable_spans(
                    called.alignment.columns(),
                    &called.evidence,
                    &regions,
                )),
                mixtures: Some(
                    mixtures
                        .get(&read.input_sha256)
                        .cloned()
                        .unwrap_or_default(),
                ),
            })
        })
        .collect::<Result<Vec<_>>>()?;
    for (read, name) in analysis.evidence.rejected_reads.iter().zip(rejected_names) {
        let rejected =
            analysis
                .rejections
                .get(&read.input_sha256)
                .ok_or(ReportError::Inconsistent(
                    "rejected review read lacks its rejection cause",
                ))?;
        reads.push(ReviewReadResult {
            name,
            sha256: read.input_sha256.clone(),
            strand: None,
            admitted: false,
            rejection: Some(rejected.cause),
            masked: None,
            unreliable: None,
            mixtures: None,
        });
    }

    let (status, error, call) = match call {
        Ok(call) => (ReviewStatus::Called, None, Some(call)),
        Err(error) => (ReviewStatus::CallFailed, Some(error.to_string()), None),
    };
    let callable = merge(match &call {
        Some(call) => segments(
            call.callable_segments
                .iter()
                .map(|segment| (segment.start_0based, segment.end_0based_exclusive)),
        ),
        None => segments(
            analysis
                .evidence
                .reads
                .iter()
                .flat_map(|read| &read.alignment.callable_reference_segments)
                .map(span),
        ),
    });
    let notation = call.as_ref().and_then(|call| call.notation.as_deref());
    // The called sample's analysed regions: the reportable regions without
    // the stretches at their edges the call does not analyse.
    let (analysed, not_analysed) = if call.is_some() {
        analysed_regions(&analysis, &regions)
    } else {
        (regions.clone(), Vec::new())
    };
    let items = Sample::new(
        &analysis,
        &names,
        &reference,
        &analysed,
        &callable,
        notation,
        signal.as_ref().map(|(signal, _)| signal),
    )?
    .items()?;
    let calls = call.as_ref().and_then(|call| {
        call.notation.as_ref().map(|calls| {
            calls
                .iter()
                .map(|call| ReviewCallResult {
                    call: call.text.clone(),
                    position: call.position,
                    ordinal: call.ordinal,
                    reference: call.reference,
                    alternate: call.alternate,
                    convention: call.convention.as_ref().map(|convention| ConventionResult {
                        name: convention.name.clone(),
                        label: convention.label.label(),
                    }),
                })
                .collect()
        })
    });
    let windows = profile
        .windows
        .iter()
        .map(|window| (window.start + 1, window.end(), window.name.as_str()))
        .collect::<Vec<_>>();
    let zones = call
        .as_ref()
        .map(|_| zones(&analysed, &windows, &callable, &reads, &items));
    Ok(ReviewResult {
        schema_version: "dna.review/v1",
        sample_id,
        status,
        error,
        provenance,
        reads,
        calls,
        callable: call.as_ref().map(|_| {
            callable
                .iter()
                .map(|&[start, end]| IntervalResult { start, end })
                .collect()
        }),
        summary: summarize(
            &items,
            status,
            &regions,
            zones.as_deref().unwrap_or_default(),
            if call.is_some() { &reference_only } else { &[] },
            not_analysed
                .iter()
                .map(|stretch| stretch.last - stretch.first + 1)
                .sum(),
        ),
        reference_only: call.as_ref().map(|_| reference_only.clone()),
        not_analysed: call.as_ref().map(|_| not_analysed.clone()),
        empop: call.and_then(|call| call.empop).map(project_empop),
        zones,
        items,
    })
}

/// The review list of an analysed sample whose consensus was called, by the
/// rules of `dna review`: each check item's first and last reference
/// positions (1-based) and label names, and each `tract_length` label with
/// its item's span.
pub(crate) fn list(
    analysis: &SampleAnalysis,
    call: &ConsensusCall,
    signal: Option<&RunSignal>,
) -> Result<ReviewList> {
    let reference = &analysis.loaded.reference;
    let regions = clip(&analysis.loaded.profile.regions, reference.len());
    let (names, _) = reviewer_read_names(&analysis.evidence)?;
    let callable =
        merge(segments(call.callable_segments.iter().map(|segment| {
            (segment.start_0based, segment.end_0based_exclusive)
        })));
    let (analysed, _) = analysed_regions(analysis, &regions);
    let items = Sample::new(
        analysis,
        &names,
        reference,
        &analysed,
        &callable,
        call.notation.as_deref(),
        signal,
    )?
    .items()?;
    let mut review = ReviewList::default();
    for item in &items {
        let (first, last) = (item.first.position, item.last.position);
        if item.kind == ItemKind::Check {
            review.checks.push(ReviewCheck {
                first,
                last,
                labels: item
                    .labels
                    .iter()
                    .map(|label| label.kind().name())
                    .collect::<BTreeSet<_>>()
                    .into_iter()
                    .collect(),
            });
        }
        for label in &item.labels {
            if let LabelResult::TractLength { tract, runs, .. } = label {
                review.notes.push(TractNote {
                    first,
                    last,
                    tract: tract.clone(),
                    runs: runs
                        .iter()
                        .map(|run| TractNoteRun {
                            run: run.run,
                            called: run.called,
                            status: run.status,
                            reading: run.reading,
                            convention: run.convention.as_ref().map(|convention| convention.length),
                        })
                        .collect(),
                });
            }
        }
    }
    Ok(review)
}

/// The reportable `regions` (1-based, inclusive) without the stretches at
/// their ends the call does not analyse (SRS-CONS-012), and those
/// stretches. A region left out whole is dropped.
fn analysed_regions(
    analysis: &SampleAnalysis,
    regions: &[[usize; 2]],
) -> (Vec<[usize; 2]>, Vec<StretchResult>) {
    let stretches = not_analysed(&analysis.consensus, regions);
    let left_out = |position: usize| {
        stretches
            .iter()
            .any(|stretch| stretch.first_0based < position && position <= stretch.last_0based + 1)
    };
    let mut analysed: Vec<[usize; 2]> = Vec::new();
    for &[start, end] in regions {
        for position in (start..=end).filter(|&position| !left_out(position)) {
            match analysed.last_mut() {
                Some(last) if last[1] + 1 == position && last[0] >= start => last[1] = position,
                _ => analysed.push([position, position]),
            }
        }
    }
    let stretches = stretches
        .iter()
        .map(|stretch| StretchResult {
            first: stretch.first_0based + 1,
            last: stretch.last_0based + 1,
            reason: stretch.reason.label(),
        })
        .collect();
    (analysed, stretches)
}

/// One sample's evidence, indexed for the review rules.
struct Sample<'a> {
    /// The review plugin's thresholds (`[review]`).
    config: ReviewConfig,
    /// The profile's report conventions: a run length one covers is not
    /// recorded, so its estimate leaves nothing to assure.
    conventions: &'a [Convention],
    sequence: &'a [u8],
    /// The analysed regions: the reportable regions without the stretches
    /// the call does not analyse (SRS-CONS-012).
    regions: &'a [[usize; 2]],
    /// The profile's reportable regions, whose ends are the amplicons'.
    reportable: Vec<[usize; 2]>,
    /// The stretches at the reportable regions' ends the call does not
    /// analyse (SRS-CONS-012), 1-based and inclusive.
    not_analysed: Vec<[usize; 2]>,
    callable: &'a [[usize; 2]],
    /// Reference segments of the admitted reads, merged, 1-based inclusive.
    covered: Vec<[usize; 2]>,
    calls: &'a [NotationCall],
    /// Whether DNA's calls are in the profile notation; without them, a
    /// difference is the consensus's decision.
    notation: bool,
    sites: &'a [ConsensusSite],
    /// The positions (0-based) reads observe only inside their unreliable
    /// spans (SRS-CONS-012).
    unreliable_only: &'a BTreeMap<usize, EvidenceReason>,
    /// Per read, in evidence order: name, content SHA-256, strand, the core's
    /// trust in each call, and the Sanger evidence.
    reads: Vec<Read<'a>>,
    /// Read index by content SHA-256, as consensus sites name reads.
    by_id: BTreeMap<&'a str, usize>,
    /// Every covered position's observations.
    loci: BTreeMap<usize, LocusObservations>,
    /// The run signal the calls are judged against, the sample's own reads
    /// taken out, when one is given.
    signal: Option<RunSignal>,
    /// The second-channel threshold's floor and run-level parameters.
    second: &'a SignalProcessingConfig,
    /// The profile's anchored repeat windows and length tracts, whose run
    /// lengths the reads' tract decodes judge (SRS-API-020).
    length_windows: Vec<LengthWindow>,
    /// The profile's length tracts.
    length_tracts: &'a [LengthTract],
    /// The profile's nomenclature windows and length tracts, 1-based and
    /// inclusive: where a run's length is judged by its own rules.
    repeat_windows: Vec<[usize; 2]>,
    /// The length convention, when configured (`[normalization]`).
    length_convention: Option<ShortestDominant>,
    /// The reference and the core's and the consensus's configuration, for
    /// the consensus's reading of each read's runs (ADR-0074).
    reference: &'a Reference,
    call_config: &'a CallConfig,
    consensus_config: &'a ConsensusConfig,
}

struct Read<'a> {
    name: &'a str,
    /// The read's alignment columns, in reference order.
    columns: &'a [AlignmentColumn],
    sha256: &'a str,
    strand: Orientation,
    trust: Vec<CallTrust>,
    sanger: &'a SangerRead,
    /// The read's calls as its evidence gives them, in read order.
    evidence: &'a ReadEvidence,
    /// The read as the core called it.
    called: &'a CalledRead,
}

impl<'a> Sample<'a> {
    fn new(
        analysis: &'a SampleAnalysis,
        names: &'a [String],
        reference: &'a Reference,
        regions: &'a [[usize; 2]],
        callable: &'a [[usize; 2]],
        notation: Option<&'a [NotationCall]>,
        signal: Option<&RunSignal>,
    ) -> Result<Self> {
        let signal = signal.map(|signal| own(analysis, signal)).transpose()?;
        let reads = analysis
            .evidence
            .reads
            .iter()
            .zip(names)
            .map(|(read, name)| {
                let missing = || ReportError::Inconsistent("review read lacks its evidence");
                let called = analysis
                    .called
                    .iter()
                    .find(|called| called.input_sha256 == read.input_sha256)
                    .ok_or_else(missing)?;
                Ok(Read {
                    name,
                    columns: called.alignment.columns(),
                    evidence: &called.evidence,
                    called,
                    sha256: &read.input_sha256,
                    strand: read.alignment.orientation,
                    trust: analysis.trust(&read.input_sha256).ok_or_else(missing)?,
                    sanger: analysis
                        .sanger
                        .get(&read.input_sha256)
                        .ok_or_else(missing)?,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let by_id = analysis
            .evidence
            .reads
            .iter()
            .enumerate()
            .map(|(index, read)| (read.input_sha256.as_str(), index))
            .collect();
        let covered = merge(segments(
            analysis
                .evidence
                .reads
                .iter()
                .flat_map(|read| &read.alignment.reference_segments)
                .map(span),
        ));
        // The reportable positions outside the analysed `regions`: the
        // stretches the call does not analyse.
        let reportable = clip(&analysis.loaded.profile.regions, reference.len());
        let mut not_analysed: Vec<[usize; 2]> = Vec::new();
        for position in reportable
            .iter()
            .flat_map(|&[start, end]| start..=end)
            .filter(|&position| !inside(regions, position))
        {
            match not_analysed.last_mut() {
                Some(last) if last[1] + 1 == position => last[1] = position,
                _ => not_analysed.push([position, position]),
            }
        }
        let loci = analysis
            .observations(1, reference.len())?
            .into_iter()
            .map(|locus| (locus.position_1based, locus))
            .collect();
        Ok(Self {
            config: analysis.loaded.config.review,
            conventions: &analysis.loaded.profile.conventions,
            sequence: reference.sequence.as_bytes(),
            regions,
            reportable,
            not_analysed,
            callable,
            covered,
            calls: notation.unwrap_or_default(),
            notation: notation.is_some(),
            sites: &analysis.consensus.sites,
            unreliable_only: &analysis.consensus.unreliable_only,
            reads,
            by_id,
            loci,
            signal,
            second: &analysis.loaded.config.sanger.signal_processing,
            length_windows: windows(&analysis.loaded.profile),
            length_tracts: &analysis.loaded.profile.length_tracts,
            repeat_windows: analysis
                .loaded
                .profile
                .windows
                .iter()
                .map(|window| [window.start + 1, window.end()])
                .chain(
                    analysis
                        .loaded
                        .profile
                        .length_tracts
                        .iter()
                        .map(|tract| tract.window),
                )
                .collect(),
            length_convention: analysis.loaded.config.length_convention,
            reference,
            call_config: &analysis.loaded.config.core,
            consensus_config: &analysis.loaded.config.consensus,
        })
    }

    fn read(&self, index: usize) -> Result<&Read<'a>> {
        Ok(self
            .reads
            .get(index)
            .ok_or(ReportError::MissingRead { index })?)
    }

    fn base(&self, position: usize) -> char {
        char::from(self.sequence[position - 1])
    }
}

/// A sequence's runs of one base, in order, with their lengths.
pub(super) fn runs_of(sequence: &str) -> Vec<(char, usize)> {
    let mut runs: Vec<(char, usize)> = Vec::new();
    for base in sequence.chars() {
        match runs.last_mut() {
            Some((last, length)) if *last == base => *length += 1,
            _ => runs.push((base, 1)),
        }
    }
    runs
}

/// A share rounded to three decimals.
pub(super) fn round(share: f64) -> f64 {
    (share * 1000.0).round() / 1000.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shares_round_to_three_decimals() {
        assert_eq!(round(0.153_49), 0.153);
        assert_eq!(round(0.153_5), 0.154);
    }
}
