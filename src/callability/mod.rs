//! Signal-derived read callability (`dna.read_callability/v1`).
//!
//! From the read's own signal, in trace order and without reference, profile,
//! or filename knowledge, the stage derives per-position features, segments the
//! read into phase states with a hysteresis state machine, classifies every
//! masked segment, and reports a typed mask with its callable span. Decoded
//! channels, loci, calls, and locus evidence are never modified.
//!
//! The core (`features`, `runs`, `phase`, `shadow`, `classify`, `mask`) depends only on
//! the plain records of `model::callability`; `sanger` is the one adapter that
//! reads Sanger types.

mod classify;
mod features;
mod mask;
mod phase;
mod runs;
mod sanger;
mod shadow;

use crate::config::Config;
use crate::error::Result;
use crate::model::basecalls::BaseCalls;
use crate::model::callability::{PositionEvidence, ReadCallability};
use crate::model::sanger::Chromatogram;
use crate::model::signal::SignalAnalysis;

/// Plain thresholds of the core, assembled from the configuration.
#[derive(Debug, Clone, Copy)]
struct Settings {
    weak_amplitude_fraction: f64,
    repeat_min_length: usize,
    shadow: shadow::Rules,
    thresholds: phase::Thresholds,
}

impl Settings {
    fn from_config(config: &Config) -> Self {
        Self {
            weak_amplitude_fraction: config.callability.weak_amplitude_fraction,
            repeat_min_length: config.callability.repeat_min_length,
            shadow: shadow::Rules {
                main_minimum: config.callability.minimum_main_share,
                far_maximum: config.callability.maximum_far_share,
                shadow_minimum: config.callability.minimum_shadow_share,
            },
            thresholds: phase::Thresholds {
                window: config.callability.window_calls,
                onset: config.callability.onset_defect_fraction,
                exit: config.callability.exit_defect_fraction,
                double_ratio: config.basecalling.secondary_peak_ratio,
            },
        }
    }
}

/// Derives the callability of one Sanger read.
pub(crate) fn analyze(
    trace: &Chromatogram,
    calls: &BaseCalls,
    signal: &SignalAnalysis,
    config: &Config,
) -> Result<ReadCallability> {
    let evidence = sanger::evidence(trace, calls, signal)?;
    derive(&evidence, &Settings::from_config(config))
}

/// The modality-generic core: evidence records in call order to callability.
fn derive(evidence: &[PositionEvidence], settings: &Settings) -> Result<ReadCallability> {
    let features = features::calculate(evidence, settings.weak_amplitude_fraction)?;
    let primary = evidence
        .iter()
        .map(|position| position.primary)
        .collect::<Vec<_>>();
    let repeats = runs::find(&primary, settings.repeat_min_length);
    let prior = phase::prior_windows(features.len(), &repeats, settings.thresholds.window);
    let defects = phase::defects(&features, &settings.thresholds);
    let callable = phase::segment(&defects, &prior, &repeats, &settings.thresholds);
    let segments = classify::segments(
        &callable,
        &classify::Context {
            evidence,
            features: &features,
            defects: &defects,
            repeats: &repeats,
            window: settings.thresholds.window,
            rules: settings.shadow,
        },
    );
    mask::build(repeats, segments, features.len())
}

#[cfg(test)]
mod tests {
    use crate::model::callability::{PhaseState, RepeatUnit};

    use super::*;

    const SETTINGS: Settings = Settings {
        weak_amplitude_fraction: 0.1,
        repeat_min_length: 8,
        shadow: shadow::Rules {
            main_minimum: 0.35,
            far_maximum: 0.12,
            shadow_minimum: 0.1,
        },
        thresholds: phase::Thresholds {
            window: 8,
            onset: 0.375,
            exit: 0.125,
            double_ratio: 0.33,
        },
    };

    fn channel(base: char) -> usize {
        match base {
            'A' => 0,
            'C' => 1,
            'G' => 2,
            'T' => 3,
            other => panic!("unsupported synthetic base {other}"),
        }
    }

    /// A clean ladder spelling `sequence`; from `shadow_from` on, every
    /// position also carries each `(offset, fraction)` shadow: that fraction of
    /// the base `offset` calls away, added to its channel.
    fn ladder(
        sequence: &str,
        shadow_from: usize,
        shadows: &[(isize, f64)],
    ) -> Vec<PositionEvidence> {
        let bases = sequence.chars().collect::<Vec<_>>();
        bases
            .iter()
            .enumerate()
            .map(|(index, &base)| {
                let mut amplitudes = [0.0; 4];
                amplitudes[channel(base)] = 1_000.0;
                if index >= shadow_from {
                    for &(offset, fraction) in shadows {
                        if let Some(&shadow) = index
                            .checked_add_signed(offset)
                            .and_then(|neighbour| bases.get(neighbour))
                        {
                            amplitudes[channel(shadow)] += 1_000.0 * fraction;
                        }
                    }
                }
                PositionEvidence {
                    amplitudes,
                    coordinate: 2 + index * 12,
                    primary: Some(channel(base)),
                }
            })
            .collect()
    }

    #[test]
    fn keeps_a_clean_read_entirely_in_phase() -> Result<()> {
        let evidence = ladder("ACGTCAGTACGATCGTACCTGAGTACGA", usize::MAX, &[]);
        let callability = derive(&evidence, &SETTINGS)?;
        assert_eq!(callability.segments.len(), 1);
        assert_eq!(callability.segments[0].state, PhaseState::InPhase);
        assert_eq!(callability.masked_count(), 0);
        assert_eq!(
            (
                callability.callable_start_0based,
                callability.callable_end_0based_exclusive
            ),
            (0, 28)
        );
        assert!(callability.repeats.is_empty());
        Ok(())
    }

    #[test]
    fn masks_a_shadow_ladder_after_a_long_homopolymer_as_dephased() -> Result<()> {
        let sequence = format!("ACGTAGTCAGTACG{}TAGCTAGCATGCATGACTGACTAG", "C".repeat(9));
        let evidence = ladder(&sequence, 23, &[(-1, 0.4)]);
        let callability = derive(&evidence, &SETTINGS)?;
        assert_eq!(
            callability.repeats,
            [crate::model::callability::RepeatRun {
                call_start_0based: 14,
                call_end_0based_exclusive: 23,
                unit: RepeatUnit::Homopolymer(1),
            }]
        );
        assert_eq!(callability.segments.len(), 2, "{:?}", callability.segments);
        assert_eq!(callability.segments[0].state, PhaseState::InPhase);
        assert_eq!(callability.segments[0].call_end_0based_exclusive, 23);
        assert_eq!(callability.segments[1].state, PhaseState::Dephased);
        assert!(callability.segments[1].after_repeat);
        assert_eq!(
            callability.segments[1]
                .shadow
                .map(|shadow| shadow.offsets().collect::<Vec<_>>()),
            Some(vec![-1])
        );
        assert_eq!(callability.callable_end_0based_exclusive, 23);
        Ok(())
    }

    #[test]
    fn masks_a_two_sided_shadow_ladder_after_a_long_homopolymer_as_dephased() -> Result<()> {
        let sequence = format!("ACGTAGTCAGTACG{}TAGCTAGCATGCATGACTGACTAG", "C".repeat(9));
        let evidence = ladder(&sequence, 23, &[(-1, 0.4), (1, 0.4)]);
        let callability = derive(&evidence, &SETTINGS)?;
        let last = callability.segments.last().copied();
        assert_eq!(
            last.map(|segment| segment.state),
            Some(PhaseState::Dephased)
        );
        assert_eq!(last.map(|segment| segment.after_repeat), Some(true));
        assert_eq!(
            last.and_then(|segment| segment.shadow)
                .map(|shadow| shadow.offsets().collect::<Vec<_>>()),
            Some(vec![-1, 1])
        );
        Ok(())
    }

    #[test]
    fn relabels_without_moving_segment_boundaries() -> Result<()> {
        let sequence = format!("ACGTAGTCAGTACG{}TAGCTAGCATGCATGACTGACTAG", "C".repeat(9));
        let evidence = ladder(&sequence, 23, &[(-1, 0.4), (1, 0.4)]);
        let strict = Settings {
            shadow: shadow::Rules {
                main_minimum: 1.0,
                far_maximum: 0.0,
                shadow_minimum: 1.0,
            },
            ..SETTINGS
        };
        let relaxed = derive(&evidence, &SETTINGS)?;
        let strict = derive(&evidence, &strict)?;
        let bounds = |callability: &ReadCallability| {
            callability
                .segments
                .iter()
                .map(|segment| (segment.call_start_0based, segment.call_end_0based_exclusive))
                .collect::<Vec<_>>()
        };
        assert_eq!(bounds(&relaxed), bounds(&strict));
        assert_eq!(relaxed.mask.len(), strict.mask.len());
        assert_eq!(relaxed.masked_count(), strict.masked_count());
        assert_eq!(
            strict.segments.last().map(|segment| segment.state),
            Some(PhaseState::Mixed)
        );
        Ok(())
    }

    #[test]
    fn masks_scattered_double_peaks_as_mixed_without_a_repeat() -> Result<()> {
        let mut evidence = ladder("ACGTCAGTACGATCGTACCTGAGTACGATCGATCGTAGCT", usize::MAX, &[]);
        for (offset, position) in evidence.iter_mut().enumerate().skip(20) {
            // A secondary channel unrelated to the neighbours' primaries.
            let primary = position.primary.unwrap_or(0);
            let secondary = (primary + 2 + offset % 2) % 4;
            position.amplitudes[secondary] = 500.0;
        }
        let callability = derive(&evidence, &SETTINGS)?;
        let masked = callability
            .segments
            .iter()
            .filter(|segment| segment.state != PhaseState::InPhase)
            .collect::<Vec<_>>();
        assert_eq!(masked.len(), 1, "{:?}", callability.segments);
        assert_eq!(masked[0].state, PhaseState::Mixed);
        assert!(!masked[0].after_repeat);
        assert_eq!(callability.callable_end_0based_exclusive, 20);
        Ok(())
    }

    #[test]
    fn masks_a_weak_tail() -> Result<()> {
        let mut evidence = ladder("ACGTCAGTACGATCGTACCTGAGTACGATCGA", usize::MAX, &[]);
        for position in evidence.iter_mut().skip(24) {
            position.amplitudes = position.amplitudes.map(|amplitude| amplitude / 100.0);
        }
        let callability = derive(&evidence, &SETTINGS)?;
        assert_eq!(
            callability.segments.last().map(|segment| segment.state),
            Some(PhaseState::Weak)
        );
        assert_eq!(callability.callable_end_0based_exclusive, 24);
        Ok(())
    }
}
