//! Resolution of variant call mappings to reference-strand call evidence.

use crate::error::CallEvidenceError;
use crate::model::alignment::Orientation;
use crate::model::basecalls::BaseCalls;
use crate::model::nucleotide::is_canonical;
use crate::model::quality::QualityControlResult;
use crate::model::variant::{VariantCallMapping, VariantCallRole};

/// One mapped call's primary base, primary-event peak heights, and relative
/// quality, expressed on the reference strand.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ReferenceCallEvidence {
    pub(crate) base: char,
    pub(crate) peak_heights: [i32; 4],
    pub(crate) quality: u8,
}

/// One variant call mapping that carries public reference-strand evidence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct PublicCall {
    pub(crate) mapping: VariantCallMapping,
    pub(crate) evidence: ReferenceCallEvidence,
}

/// Resolves the public call evidence of one variant.
///
/// A flanking call whose primary base is unresolved (`N`: its strongest channels
/// tie, carry no signal, or all qualify) has no called base to report, so it is
/// omitted rather than fabricated. Supporting calls are canonical by
/// construction and must resolve, and at least one call must remain.
pub(crate) fn resolve_public_calls(
    calls: &BaseCalls,
    quality: &QualityControlResult,
    orientation: Orientation,
    mappings: &[VariantCallMapping],
) -> Result<Vec<PublicCall>, CallEvidenceError> {
    let mut public = Vec::with_capacity(mappings.len());
    for &mapping in mappings {
        match resolve(calls, quality, orientation, mapping.call_index_0based) {
            Ok(evidence) => public.push(PublicCall { mapping, evidence }),
            Err(CallEvidenceError::UnresolvedCall { .. })
                if mapping.role == VariantCallRole::Flanking => {}
            Err(error) => return Err(error),
        }
    }
    if public.is_empty() {
        return Err(CallEvidenceError::NoResolvedCalls);
    }
    Ok(public)
}

/// Resolves call `index` against its call and quality records, requires a
/// canonical primary base with primary-event peak evidence, and projects it
/// onto the reference strand.
fn resolve(
    calls: &BaseCalls,
    quality: &QualityControlResult,
    orientation: Orientation,
    index: usize,
) -> Result<ReferenceCallEvidence, CallEvidenceError> {
    let call = calls
        .calls
        .get(index)
        .ok_or(CallEvidenceError::MissingCall { index })?;
    let score = quality
        .per_call
        .get(index)
        .ok_or(CallEvidenceError::MissingQuality { index })?;
    if call.index_0based != index || score.index_0based != index {
        return Err(CallEvidenceError::IndexMismatch { index });
    }
    if !is_canonical(call.primary) {
        return Err(CallEvidenceError::UnresolvedCall { index });
    }
    let primary = call
        .primary_peak_evidence
        .as_ref()
        .ok_or(CallEvidenceError::MissingPeakEvidence { index })?;
    Ok(ReferenceCallEvidence {
        base: orientation.reference_base(call.primary),
        peak_heights: orientation.reference_peak_heights(primary.channel_heights),
        quality: score.relative_quality_score,
    })
}

#[cfg(test)]
mod tests {
    use crate::model::basecalls::{BaseCall, ChannelPeak, PeakSource, PrimaryPeakEvidence};
    use crate::model::nucleotide::Nucleotide;
    use crate::model::quality::CallQuality;
    use crate::model::variant::VariantCallRole;

    use super::*;

    fn call(index: usize, primary: char, peak: Option<[i32; 4]>) -> BaseCall {
        BaseCall {
            index_0based: index,
            locus_position_0based: index * 4 + 2,
            window_start_0based: index * 4,
            window_end_0based_exclusive: index * 4 + 4,
            peaks: Nucleotide::ALL.map(|base| ChannelPeak {
                base,
                height: 0,
                position_0based: index * 4 + 2,
                source: PeakSource::LocalMaximum,
            }),
            primary_peak_evidence: peak.map(|channel_heights| PrimaryPeakEvidence {
                position_0based: index * 4 + 2,
                channel_heights,
            }),
            primary,
            ambiguity: primary,
            qualifying_channels: Vec::new(),
            vendor_agrees: None,
        }
    }

    fn records(
        calls: Vec<BaseCall>,
        quality_indices: &[usize],
    ) -> (BaseCalls, QualityControlResult) {
        let primary_sequence = calls.iter().map(|call| call.primary).collect();
        let per_call = quality_indices
            .iter()
            .map(|&index| CallQuality {
                index_0based: index,
                penalty: 0,
                relative_quality_score: 40,
                vendor_quality_applies: false,
            })
            .collect();
        (
            BaseCalls {
                calls,
                primary_sequence,
            },
            QualityControlResult {
                per_call,
                trim_start_0based: 0,
                trim_end_0based_exclusive: quality_indices.len(),
                retained_sequence: String::new(),
            },
        )
    }

    #[test]
    fn projects_reverse_reads_onto_the_reference_strand() {
        let (calls, quality) = records(vec![call(0, 'A', Some([10, 20, 30, 40]))], &[0]);

        let forward = resolve(&calls, &quality, Orientation::Forward, 0);
        let reverse = resolve(&calls, &quality, Orientation::Reverse, 0);

        assert_eq!(
            forward,
            Ok(ReferenceCallEvidence {
                base: 'A',
                peak_heights: [10, 20, 30, 40],
                quality: 40,
            })
        );
        assert_eq!(
            reverse,
            Ok(ReferenceCallEvidence {
                base: 'T',
                peak_heights: [40, 30, 20, 10],
                quality: 40,
            })
        );
    }

    fn mapping(role: VariantCallRole, index: usize) -> VariantCallMapping {
        VariantCallMapping {
            role,
            call_index_0based: index,
            reference_position_0based: Some(index),
        }
    }

    #[test]
    fn omits_unresolved_flanks_from_public_calls() {
        let tied = call(0, 'N', None);
        let mixed = call(2, 'N', Some([287, 171, 150, 280]));
        let (calls, quality) = records(
            vec![tied, call(1, 'G', Some([0, 0, 9, 1])), mixed],
            &[0, 1, 2],
        );
        let mappings = [
            mapping(VariantCallRole::Flanking, 0),
            mapping(VariantCallRole::Flanking, 1),
            mapping(VariantCallRole::Flanking, 2),
        ];

        let resolved = resolve_public_calls(&calls, &quality, Orientation::Forward, &mappings);

        assert_eq!(
            resolved,
            Ok(vec![PublicCall {
                mapping: mappings[1],
                evidence: ReferenceCallEvidence {
                    base: 'G',
                    peak_heights: [0, 0, 9, 1],
                    quality: 40,
                },
            }])
        );
    }

    #[test]
    fn requires_evidence_for_supporting_calls_and_at_least_one_public_call() {
        let (calls, quality) = records(vec![call(0, 'N', None), call(1, 'N', None)], &[0, 1]);

        assert_eq!(
            resolve_public_calls(
                &calls,
                &quality,
                Orientation::Forward,
                &[mapping(VariantCallRole::Supporting, 0)],
            ),
            Err(CallEvidenceError::UnresolvedCall { index: 0 })
        );
        assert_eq!(
            resolve_public_calls(
                &calls,
                &quality,
                Orientation::Forward,
                &[
                    mapping(VariantCallRole::Flanking, 0),
                    mapping(VariantCallRole::Flanking, 1),
                ],
            ),
            Err(CallEvidenceError::NoResolvedCalls)
        );
    }

    #[test]
    fn rejects_unresolvable_call_mappings() {
        let (calls, quality) = records(vec![call(0, 'A', Some([1, 0, 0, 0]))], &[0]);
        assert_eq!(
            resolve(&calls, &quality, Orientation::Forward, 1),
            Err(CallEvidenceError::MissingCall { index: 1 })
        );

        let (calls, quality) = records(vec![call(0, 'A', Some([1, 0, 0, 0]))], &[]);
        assert_eq!(
            resolve(&calls, &quality, Orientation::Forward, 0),
            Err(CallEvidenceError::MissingQuality { index: 0 })
        );

        let (calls, quality) = records(vec![call(0, 'A', Some([1, 0, 0, 0]))], &[5]);
        assert_eq!(
            resolve(&calls, &quality, Orientation::Forward, 0),
            Err(CallEvidenceError::IndexMismatch { index: 0 })
        );

        let (calls, quality) = records(vec![call(0, 'N', None)], &[0]);
        assert_eq!(
            resolve(&calls, &quality, Orientation::Forward, 0),
            Err(CallEvidenceError::UnresolvedCall { index: 0 })
        );

        let (calls, quality) = records(vec![call(0, 'A', None)], &[0]);
        assert_eq!(
            resolve(&calls, &quality, Orientation::Forward, 0),
            Err(CallEvidenceError::MissingPeakEvidence { index: 0 })
        );
    }
}
