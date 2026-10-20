//! Sanger evidence joined to core call records by call index.

use crate::error::CallEvidenceError;
use crate::model::sanger_read::SangerCallSignal;
use dna_align::model::alignment::Orientation;
use dna_sanger::model::basecalls::BaseCalls;
use dna_sanger::model::quality::QualityControlResult;

/// Primary-event peak heights on the reference strand and relative quality of
/// one Sanger call.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct SangerCallEvidence {
    pub(super) peak_heights: [i32; 4],
    pub(super) quality: u8,
}

/// Resolves call `index` against its call and quality records and projects its
/// primary-event peak heights onto the reference strand.
pub(super) fn evidence(
    calls: &BaseCalls,
    quality: &QualityControlResult,
    orientation: Orientation,
    index: usize,
) -> Result<SangerCallEvidence, CallEvidenceError> {
    on_reference(signal(calls, quality, index)?, orientation, index)
}

/// The signal of call `index` from its call and quality records, on the read
/// strand.
pub(super) fn signal(
    calls: &BaseCalls,
    quality: &QualityControlResult,
    index: usize,
) -> Result<SangerCallSignal, CallEvidenceError> {
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
    Ok(SangerCallSignal {
        peaks: call
            .primary_peak_evidence
            .as_ref()
            .map(|primary| primary.channel_heights),
        quality: score.relative_quality_score,
    })
}

/// The signal of call `index` with its primary-event peak heights projected
/// onto the reference strand.
pub(super) fn on_reference(
    signal: SangerCallSignal,
    orientation: Orientation,
    index: usize,
) -> Result<SangerCallEvidence, CallEvidenceError> {
    let peaks = signal
        .peaks
        .ok_or(CallEvidenceError::MissingPeakEvidence { index })?;
    Ok(SangerCallEvidence {
        peak_heights: reference_peak_heights(orientation, peaks),
        quality: signal.quality,
    })
}

/// Projects A/C/G/T channel heights from trace strand to reference strand.
const fn reference_peak_heights(orientation: Orientation, peaks: [i32; 4]) -> [i32; 4] {
    match orientation {
        Orientation::Forward => peaks,
        Orientation::Reverse => [peaks[3], peaks[2], peaks[1], peaks[0]],
    }
}

#[cfg(test)]
mod tests {
    use dna_kernel::model::nucleotide::Nucleotide;
    use dna_sanger::model::basecalls::{BaseCall, ChannelPeak, PeakSource, PrimaryPeakEvidence};
    use dna_sanger::model::quality::CallQuality;

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
    fn projects_reverse_peaks_onto_the_reference_strand() {
        let (calls, quality) = records(vec![call(0, 'A', Some([10, 20, 30, 40]))], &[0]);

        assert_eq!(
            evidence(&calls, &quality, Orientation::Forward, 0),
            Ok(SangerCallEvidence {
                peak_heights: [10, 20, 30, 40],
                quality: 40,
            })
        );
        assert_eq!(
            evidence(&calls, &quality, Orientation::Reverse, 0),
            Ok(SangerCallEvidence {
                peak_heights: [40, 30, 20, 10],
                quality: 40,
            })
        );
    }

    #[test]
    fn rejects_unresolvable_call_records() {
        let (calls, quality) = records(vec![call(0, 'A', Some([1, 0, 0, 0]))], &[0]);
        assert_eq!(
            evidence(&calls, &quality, Orientation::Forward, 1),
            Err(CallEvidenceError::MissingCall { index: 1 })
        );

        let (calls, quality) = records(vec![call(0, 'A', Some([1, 0, 0, 0]))], &[]);
        assert_eq!(
            evidence(&calls, &quality, Orientation::Forward, 0),
            Err(CallEvidenceError::MissingQuality { index: 0 })
        );

        let (calls, quality) = records(vec![call(0, 'A', Some([1, 0, 0, 0]))], &[5]);
        assert_eq!(
            evidence(&calls, &quality, Orientation::Forward, 0),
            Err(CallEvidenceError::IndexMismatch { index: 0 })
        );

        let (calls, quality) = records(vec![call(0, 'A', None)], &[0]);
        assert_eq!(
            evidence(&calls, &quality, Orientation::Forward, 0),
            Err(CallEvidenceError::MissingPeakEvidence { index: 0 })
        );
    }
}
