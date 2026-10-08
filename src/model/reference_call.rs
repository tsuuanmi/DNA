//! Resolution of variant call mappings to reference-strand call evidence.

use crate::error::CallEvidenceError;
use crate::model::alignment::Orientation;
use crate::model::basecalls::BaseCalls;
use crate::model::quality::QualityControlResult;

/// One mapped call's primary base, primary-event peak heights, and relative
/// quality, expressed on the reference strand.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ReferenceCallEvidence {
    pub(crate) base: char,
    pub(crate) peak_heights: [i32; 4],
    pub(crate) quality: u8,
}

/// Resolves call `index` against its call and quality records, requires
/// primary-event peak evidence, and projects it onto the reference strand.
pub(crate) fn resolve(
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
            Err(CallEvidenceError::MissingPeakEvidence { index: 0 })
        );
    }
}
