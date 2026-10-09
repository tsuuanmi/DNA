//! Segment classification: the phase state of every masked run of positions
//! and whether it follows a repeat run.

use crate::model::callability::{PhaseSegment, PhaseState, PositionFeatures, RepeatRun};
use crate::signal_processing::round_metric;

use super::features::SHIFT_OFFSETS;
use super::phase::Defect;

/// Turns per-position callable flags into ordered segments, classifying every
/// masked run by its dominant defect class.
pub(super) fn segments(
    callable: &[bool],
    defects: &[Defect],
    features: &[PositionFeatures],
    repeats: &[RepeatRun],
    window: usize,
    shift_coherence: f64,
) -> Vec<PhaseSegment> {
    let mut segments = Vec::new();
    let mut start = 0;
    while start < callable.len() {
        let flag = callable[start];
        let end = start
            + callable[start..]
                .iter()
                .take_while(|&&next| next == flag)
                .count();
        segments.push(if flag {
            PhaseSegment {
                call_start_0based: start,
                call_end_0based_exclusive: end,
                state: PhaseState::InPhase,
                after_repeat: false,
                coherence: None,
                modal_offset: None,
            }
        } else {
            classify(
                start,
                end,
                defects,
                features,
                repeats,
                window,
                shift_coherence,
            )
        });
        start = end;
    }
    segments
}

/// Classifies one masked run.
///
/// The dominant defect class wins, with ties resolved in the order weak,
/// double, spacing. A double-peak run is dephased when at least
/// `shift_coherence` of its double peaks share the modal shift offset, and
/// mixed otherwise.
fn classify(
    start: usize,
    end: usize,
    defects: &[Defect],
    features: &[PositionFeatures],
    repeats: &[RepeatRun],
    window: usize,
    shift_coherence: f64,
) -> PhaseSegment {
    let mut weak = 0_usize;
    let mut double = 0_usize;
    let mut spacing = 0_usize;
    let mut offsets = [0_usize; SHIFT_OFFSETS.len()];
    for index in start..end {
        match defects[index] {
            Defect::Weak => weak += 1,
            Defect::Double => {
                double += 1;
                if let Some(slot) = features[index]
                    .shift_offset
                    .and_then(|offset| SHIFT_OFFSETS.iter().position(|&known| known == offset))
                {
                    offsets[slot] += 1;
                }
            }
            Defect::Spacing => spacing += 1,
            Defect::None => {}
        }
    }
    let (coherence, modal_offset) = if double == 0 {
        (None, None)
    } else {
        // Ties favour the nearest offset, which comes first in the fixed order.
        let (slot, count) =
            offsets
                .iter()
                .copied()
                .enumerate()
                .fold(
                    (0, 0),
                    |best, (slot, count)| {
                        if count > best.1 { (slot, count) } else { best }
                    },
                );
        (
            Some(round_metric(count as f64 / double as f64)),
            (count > 0).then_some(SHIFT_OFFSETS[slot]),
        )
    };
    let state = if weak >= double && weak >= spacing {
        PhaseState::Weak
    } else if double >= spacing {
        if coherence.is_some_and(|coherence| coherence >= shift_coherence) {
            PhaseState::Dephased
        } else {
            PhaseState::Mixed
        }
    } else {
        PhaseState::Irregular
    };
    let after_repeat = repeats.iter().any(|run| {
        let first = run.call_end_0based_exclusive;
        start >= first && start < first.saturating_add(window)
    });
    PhaseSegment {
        call_start_0based: start,
        call_end_0based_exclusive: end,
        state,
        after_repeat,
        coherence,
        modal_offset,
    }
}

#[cfg(test)]
mod tests {
    use crate::model::callability::RepeatUnit;

    use super::*;

    fn feature(shift_offset: Option<i8>) -> PositionFeatures {
        PositionFeatures {
            primary_channel: Some(0),
            secondary_channel: Some(1),
            dominance: 0.2,
            secondary_ratio: 0.6,
            shift_offset,
            spacing_deviation: 0.0,
            weak: false,
        }
    }

    #[test]
    fn labels_coherent_doubles_dephased_and_scattered_doubles_mixed() {
        let callable = [
            true, true, true, false, false, false, false, true, true, true,
        ];
        let defects = [
            Defect::None,
            Defect::None,
            Defect::None,
            Defect::Double,
            Defect::Double,
            Defect::Double,
            Defect::Double,
            Defect::None,
            Defect::None,
            Defect::None,
        ];
        let coherent = [
            feature(None),
            feature(None),
            feature(None),
            feature(Some(1)),
            feature(Some(1)),
            feature(Some(1)),
            feature(Some(-1)),
            feature(None),
            feature(None),
            feature(None),
        ];
        let coherent = segments(&callable, &defects, &coherent, &[], 4, 0.75);
        assert_eq!(coherent.len(), 3);
        assert_eq!(coherent[0].state, PhaseState::InPhase);
        assert_eq!(coherent[1].state, PhaseState::Dephased);
        assert_eq!(coherent[1].coherence, Some(0.75));
        assert_eq!(coherent[1].modal_offset, Some(1));
        assert!(!coherent[1].after_repeat);
        assert_eq!(
            (
                coherent[2].call_start_0based,
                coherent[2].call_end_0based_exclusive
            ),
            (7, 10)
        );

        let scattered = [
            feature(None),
            feature(None),
            feature(None),
            feature(Some(1)),
            feature(Some(-2)),
            feature(None),
            feature(Some(3)),
            feature(None),
            feature(None),
            feature(None),
        ];
        let scattered = segments(&callable, &defects, &scattered, &[], 4, 0.75);
        assert_eq!(scattered[1].state, PhaseState::Mixed);
        assert_eq!(scattered[1].coherence, Some(0.25));
        assert_eq!(scattered[1].modal_offset, Some(1));
    }

    #[test]
    fn labels_weak_and_irregular_runs_and_attributes_runs_after_a_repeat() {
        let callable = [false, false, false, true, true, true, true, false, false];
        let defects = [
            Defect::Weak,
            Defect::Double,
            Defect::Weak,
            Defect::None,
            Defect::None,
            Defect::None,
            Defect::None,
            Defect::Spacing,
            Defect::None,
        ];
        let features = vec![feature(None); callable.len()];
        let repeats = [RepeatRun {
            call_start_0based: 3,
            call_end_0based_exclusive: 7,
            unit: RepeatUnit::Homopolymer(2),
        }];
        let labelled = segments(&callable, &defects, &features, &repeats, 4, 0.75);
        assert_eq!(labelled[0].state, PhaseState::Weak);
        assert_eq!(labelled[0].coherence, Some(0.0));
        assert_eq!(labelled[0].modal_offset, None);
        assert!(!labelled[0].after_repeat);
        assert_eq!(labelled[2].state, PhaseState::Irregular);
        assert!(labelled[2].after_repeat);
    }
}
