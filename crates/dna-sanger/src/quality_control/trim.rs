//! Trim interval from the callable span (`dna.callable_span_trim/v1`).

use crate::error::{QualityControlError, Result};
use crate::model::basecalls::BaseCalls;
use crate::model::callability::{PhaseState, ReadCallability};
use crate::model::quality::{CallQuality, QualityControlResult};
use crate::model::sanger::Chromatogram;
use crate::quality_control::QualityControlConfig;
use crate::quality_control::penalty;
use crate::quality_control::quality;

/// Scores calls and derives the one retained interval: the read's callable
/// span, widened on each side by up to `context_margin` calls of an adjacent
/// dephased segment. Dephased calls still read the main ladder, so they anchor
/// the alignment next to the span; calls of any other masked state carry no
/// alignment information and are trimmed.
pub(crate) fn analyze(
    trace: &Chromatogram,
    calls: &BaseCalls,
    callability: &ReadCallability,
    config: &QualityControlConfig,
    context_margin: usize,
) -> Result<QualityControlResult> {
    let span_start = callability.callable_start_0based;
    let span_end = callability.callable_end_0based_exclusive;
    if span_start >= span_end || span_end > calls.len() {
        return Err(QualityControlError::InvalidCallableSpan {
            start: span_start,
            end: span_end,
            calls: calls.len(),
        }
        .into());
    }
    let penalties = penalty::calculate(calls, config.penalty_window_size)?;
    let scores = quality::relative_scores(&penalties, config.max_relative_quality_score);
    let dephased = |index: usize| {
        callability
            .segment_at(index)
            .filter(|segment| segment.state == PhaseState::Dephased)
    };
    let trim_start = span_start
        .checked_sub(1)
        .and_then(dephased)
        .map_or(span_start, |segment| {
            span_start
                .saturating_sub(context_margin)
                .max(segment.call_start_0based)
        });
    let trim_end = dephased(span_end).map_or(span_end, |segment| {
        span_end
            .saturating_add(context_margin)
            .min(segment.call_end_0based_exclusive)
    });
    let retained_sequence = calls.primary_sequence[trim_start..trim_end].to_owned();
    let per_call = calls
        .calls
        .iter()
        .enumerate()
        .map(|(index, call)| {
            let vendor_quality = trace
                .vendor
                .qualities
                .as_ref()
                .and_then(|qualities| qualities.get(index))
                .copied();
            CallQuality {
                index_0based: index,
                penalty: penalties[index],
                relative_quality_score: scores[index],
                vendor_quality_applies: vendor_quality.is_some()
                    && call.vendor_agrees == Some(true),
            }
        })
        .collect();
    Ok(QualityControlResult {
        per_call,
        trim_start_0based: trim_start,
        trim_end_0based_exclusive: trim_end,
        retained_sequence,
    })
}

#[cfg(test)]
mod tests {
    use crate::model::basecalls::{BaseCall, ChannelPeak, PeakSource};
    use crate::model::sanger::VendorEvidence;
    use dna_kernel::model::nucleotide::Nucleotide;

    use super::*;

    // SRS-CALL-010, SRS-QC-004
    #[test]
    fn trims_to_the_callable_span_and_applies_matching_vendor_quality() -> Result<()> {
        let locations = [2, 6, 10, 14];
        let calls = BaseCalls {
            calls: locations
                .iter()
                .enumerate()
                .map(|(index, &locus_position)| BaseCall {
                    index_0based: index,
                    locus_position_0based: locus_position,
                    window_start_0based: locus_position.saturating_sub(1),
                    window_end_0based_exclusive: locus_position + 2,
                    peaks: std::array::from_fn(|channel| ChannelPeak {
                        base: Nucleotide::ALL[channel],
                        height: 1,
                        position_0based: locus_position,
                        source: PeakSource::LocusFallback,
                    }),
                    primary_peak_evidence: None,
                    primary: 'A',
                    ambiguity: 'A',
                    qualifying_channels: vec![Nucleotide::A],
                    vendor_agrees: Some(true),
                })
                .collect(),
            primary_sequence: "AAAA".into(),
        };
        let trace = Chromatogram {
            source_name: "synthetic.ab1".into(),
            source_sha256: String::new(),
            channels: std::array::from_fn(|_| vec![0; 16]),
            locus_positions: locations.to_vec(),
            duplicate_loci: 0,
            vendor: VendorEvidence {
                primary: Some("AAAA".into()),
                qualities: Some(vec![40; 4]),
            },
        };
        let config = QualityControlConfig {
            penalty_window_size: 2,
            max_relative_quality_score: 60,
            trim_context_margin: 1,
        };
        let mut callability = ReadCallability::in_phase(4);
        let result = analyze(&trace, &calls, &callability, &config, 1)?;
        assert_eq!(
            (result.trim_start_0based, result.trim_end_0based_exclusive),
            (0, 4)
        );
        assert_eq!(result.retained_sequence, "AAAA");
        callability.callable_start_0based = 2;
        callability.callable_end_0based_exclusive = 3;
        let narrowed = analyze(&trace, &calls, &callability, &config, 1)?;
        assert_eq!(
            (
                narrowed.trim_start_0based,
                narrowed.trim_end_0based_exclusive
            ),
            (2, 3)
        );
        assert_eq!(narrowed.retained_sequence, "A");
        // A dephased neighbour anchors the span; a mixed one is trimmed.
        callability.segments = [
            (0, 2, PhaseState::Dephased),
            (2, 3, PhaseState::InPhase),
            (3, 4, PhaseState::Mixed),
        ]
        .map(
            |(start, end, state)| crate::model::callability::PhaseSegment {
                call_start_0based: start,
                call_end_0based_exclusive: end,
                state,
                after_repeat: false,
                shadow: None,
            },
        )
        .to_vec();
        let anchored = analyze(&trace, &calls, &callability, &config, 1)?;
        assert_eq!(
            (
                anchored.trim_start_0based,
                anchored.trim_end_0based_exclusive
            ),
            (1, 3)
        );
        callability.callable_end_0based_exclusive = 2;
        assert!(analyze(&trace, &calls, &callability, &config, 1).is_err());
        assert!(result.per_call.iter().all(|quality| {
            quality.relative_quality_score == 60 && quality.vendor_quality_applies
        }));
        Ok(())
    }
    fn one_call_fixture(
        vendor_agrees: Option<bool>,
        qualities: Option<Vec<u8>>,
    ) -> (Chromatogram, BaseCalls) {
        let calls = BaseCalls {
            calls: vec![BaseCall {
                index_0based: 0,
                locus_position_0based: 2,
                window_start_0based: 1,
                window_end_0based_exclusive: 4,
                peaks: std::array::from_fn(|channel| ChannelPeak {
                    base: Nucleotide::ALL[channel],
                    height: 1,
                    position_0based: 2,
                    source: PeakSource::LocusFallback,
                }),
                primary_peak_evidence: None,
                primary: 'A',
                ambiguity: 'A',
                qualifying_channels: vec![Nucleotide::A],
                vendor_agrees,
            }],
            primary_sequence: "A".into(),
        };
        let trace = Chromatogram {
            source_name: "synthetic.ab1".into(),
            source_sha256: String::new(),
            channels: std::array::from_fn(|_| vec![0; 8]),
            locus_positions: vec![2],
            duplicate_loci: 0,
            vendor: VendorEvidence {
                primary: Some("T".into()),
                qualities,
            },
        };
        (trace, calls)
    }

    fn one_call_config() -> QualityControlConfig {
        QualityControlConfig {
            penalty_window_size: 1,
            max_relative_quality_score: 60,
            trim_context_margin: 0,
        }
    }

    // SRS-QC-002
    #[test]
    fn vendor_quality_applies_only_when_the_vendor_call_agrees() -> Result<()> {
        let callability = ReadCallability::in_phase(1);
        for (agrees, qualities, applies) in [
            (Some(true), Some(vec![40]), true),
            (Some(false), Some(vec![40]), false),
            (None, Some(vec![40]), false),
            (Some(true), None, false),
        ] {
            let (trace, calls) = one_call_fixture(agrees, qualities);
            let result = analyze(&trace, &calls, &callability, &one_call_config(), 0)?;
            assert_eq!(result.per_call[0].vendor_quality_applies, applies);
            // The relative score never depends on the vendor quality.
            assert_eq!(result.per_call[0].relative_quality_score, 60);
        }
        Ok(())
    }

    // SRS-QC-005
    #[test]
    fn rejects_a_callable_span_that_does_not_fit_the_calls() {
        let (trace, calls) = one_call_fixture(Some(true), None);
        for (start, end) in [(0, 2), (1, 1), (1, 0)] {
            let mut callability = ReadCallability::in_phase(1);
            callability.callable_start_0based = start;
            callability.callable_end_0based_exclusive = end;
            let error = analyze(&trace, &calls, &callability, &one_call_config(), 0).err();
            assert!(
                matches!(
                    error,
                    Some(crate::error::SangerError::QualityControl(
                        QualityControlError::InvalidCallableSpan { calls: 1, .. }
                    ))
                ),
                "{error:?}"
            );
        }
    }
}
