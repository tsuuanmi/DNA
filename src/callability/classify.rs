//! Segment classification: the phase state of every masked run of positions
//! and whether it follows a repeat run.

use std::ops::Range;

use crate::model::callability::{
    PhaseSegment, PhaseState, PositionEvidence, PositionFeatures, RepeatRun, ShadowSummary,
};

use super::phase::Defect;
use super::shadow::{self, Rules, ShadowFit};

/// Everything classification reads besides the callable flags.
pub(super) struct Context<'a> {
    pub(super) evidence: &'a [PositionEvidence],
    pub(super) features: &'a [PositionFeatures],
    pub(super) defects: &'a [Defect],
    pub(super) repeats: &'a [RepeatRun],
    pub(super) window: usize,
    pub(super) rules: Rules,
}

/// Turns per-position callable flags into ordered segments, classifying every
/// masked run by its dominant defect class.
pub(super) fn segments(callable: &[bool], context: &Context<'_>) -> Vec<PhaseSegment> {
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
                shadow: None,
            }
        } else {
            classify(start..end, context)
        });
        start = end;
    }
    segments
}

/// Classifies one masked run.
///
/// The dominant defect class wins, with ties resolved in the order weak,
/// double, spacing. A double-peak run is dephased when its shadow fit explains
/// it as the main ladder plus near shadows, and mixed otherwise.
fn classify(calls: Range<usize>, context: &Context<'_>) -> PhaseSegment {
    let mut weak = 0_usize;
    let mut double = 0_usize;
    let mut spacing = 0_usize;
    for defect in &context.defects[calls.clone()] {
        match defect {
            Defect::Weak => weak += 1,
            Defect::Double => double += 1,
            Defect::Spacing => spacing += 1,
            Defect::None => {}
        }
    }
    let (state, shadow) = if weak >= double && weak >= spacing {
        (PhaseState::Weak, None)
    } else if double >= spacing {
        double_state(
            shadow::fit(context.evidence, context.features, calls.clone()),
            &context.rules,
        )
    } else {
        (PhaseState::Irregular, None)
    };
    let after_repeat = context.repeats.iter().any(|run| {
        let first = run.call_end_0based_exclusive;
        calls.start >= first && calls.start < first.saturating_add(context.window)
    });
    PhaseSegment {
        call_start_0based: calls.start,
        call_end_0based_exclusive: calls.end,
        state,
        after_repeat,
        shadow,
    }
}

/// State of a double-peak run from its shadow fit; a run without a fit is
/// mixed.
fn double_state(fit: Option<ShadowFit>, rules: &Rules) -> (PhaseState, Option<ShadowSummary>) {
    match fit {
        Some(fit) => {
            let (dephased, summary) = shadow::summarize(&fit, rules);
            let state = if dephased {
                PhaseState::Dephased
            } else {
                PhaseState::Mixed
            };
            (state, Some(summary))
        }
        None => (PhaseState::Mixed, None),
    }
}

#[cfg(test)]
mod tests {
    use crate::model::callability::RepeatUnit;

    use super::*;

    const RULES: Rules = Rules {
        main_minimum: 0.35,
        far_maximum: 0.12,
        shadow_minimum: 0.1,
    };

    fn feature() -> PositionFeatures {
        PositionFeatures {
            primary_channel: Some(0),
            secondary_channel: Some(1),
            dominance: 0.2,
            secondary_ratio: 0.6,
            spacing_deviation: 0.0,
            weak: false,
        }
    }

    fn evidence(positions: usize) -> Vec<PositionEvidence> {
        (0..positions)
            .map(|index| PositionEvidence {
                amplitudes: [1_000.0, 0.0, 0.0, 0.0],
                coordinate: index * 12,
                primary: Some(0),
            })
            .collect()
    }

    #[test]
    fn labels_double_peak_runs_from_their_shadow_fit() {
        let dephased = ShadowFit {
            weights: [0.0, 0.0, 0.2, 0.6, 0.2, 0.0, 0.0],
        };
        let (state, summary) = double_state(Some(dephased), &RULES);
        assert_eq!(state, PhaseState::Dephased);
        assert_eq!(
            summary.map(|summary| summary.offsets().collect::<Vec<_>>()),
            Some(vec![-1, 1])
        );
        let mixed = ShadowFit {
            weights: [0.04, 0.08, 0.11, 0.39, 0.08, 0.05, 0.07],
        };
        let (state, summary) = double_state(Some(mixed), &RULES);
        assert_eq!(state, PhaseState::Mixed);
        assert!(summary.is_some());
        assert_eq!(double_state(None, &RULES), (PhaseState::Mixed, None));
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
        let features = vec![feature(); callable.len()];
        let evidence = evidence(callable.len());
        let repeats = [RepeatRun {
            call_start_0based: 3,
            call_end_0based_exclusive: 7,
            unit: RepeatUnit::Homopolymer(2),
        }];
        let context = Context {
            evidence: &evidence,
            features: &features,
            defects: &defects,
            repeats: &repeats,
            window: 4,
            rules: RULES,
        };
        let labelled = segments(&callable, &context);
        assert_eq!(labelled.len(), 3);
        assert_eq!(labelled[0].state, PhaseState::Weak);
        assert_eq!(labelled[0].shadow, None);
        assert!(!labelled[0].after_repeat);
        assert_eq!(labelled[1].state, PhaseState::InPhase);
        assert_eq!(labelled[2].state, PhaseState::Irregular);
        assert!(labelled[2].after_repeat);
    }

    #[test]
    fn labels_a_short_double_peak_run_mixed() {
        let callable = [true, true, false, false, false, true, true];
        let defects = [
            Defect::None,
            Defect::None,
            Defect::Double,
            Defect::Double,
            Defect::Double,
            Defect::None,
            Defect::None,
        ];
        let features = vec![feature(); callable.len()];
        let evidence = evidence(callable.len());
        let context = Context {
            evidence: &evidence,
            features: &features,
            defects: &defects,
            repeats: &[],
            window: 4,
            rules: RULES,
        };
        let labelled = segments(&callable, &context);
        assert_eq!(labelled[1].state, PhaseState::Mixed);
        assert_eq!(labelled[1].shadow, None);
    }
}
