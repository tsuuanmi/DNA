//! Phase segmentation: rolling defect fractions, repeat priors, and the
//! hysteresis state machine that decides where the read is callable.

use crate::model::callability::{PositionFeatures, RepeatRun};

/// Relative spacing deviation at or above which a position is a spacing defect.
pub(super) const SPACING_DEFECT: f64 = 0.5;

/// Thresholds the segmentation applies.
#[derive(Debug, Clone, Copy)]
pub(super) struct Thresholds {
    /// Calls per rolling window.
    pub(super) window: usize,
    /// Defect fraction at or above which an in-phase stretch ends.
    pub(super) onset: f64,
    /// Defect fraction at or below which a masked stretch ends.
    pub(super) exit: f64,
    /// Secondary ratio at or above which a position carries a double peak.
    pub(super) double_ratio: f64,
}

/// The one defect class a position belongs to, by precedence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Defect {
    None,
    Weak,
    Double,
    Spacing,
}

/// Marks the positions inside the window that follows each repeat run; the
/// calls of any run are never inside a prior window, so a run that follows
/// another closely keeps its full-threshold onset.
pub(super) fn prior_windows(positions: usize, repeats: &[RepeatRun], window: usize) -> Vec<bool> {
    let mut prior = vec![false; positions];
    for run in repeats {
        let start = run.call_end_0based_exclusive.min(positions);
        let end = start.saturating_add(window).min(positions);
        prior[start..end].fill(true);
    }
    for run in repeats {
        let start = run.call_start_0based.min(positions);
        let end = run.call_end_0based_exclusive.min(positions);
        prior[start..end].fill(false);
    }
    prior
}

/// Classifies every position into one defect class.
pub(super) fn defects(features: &[PositionFeatures], thresholds: &Thresholds) -> Vec<Defect> {
    features
        .iter()
        .map(|feature| {
            if feature.weak {
                Defect::Weak
            } else if feature.secondary_ratio >= thresholds.double_ratio {
                Defect::Double
            } else if feature.spacing_deviation >= SPACING_DEFECT {
                Defect::Spacing
            } else {
                Defect::None
            }
        })
        .collect()
}

/// Decides per position whether the read is callable there.
///
/// The forward statistic at `i` is the defect fraction over `[i, i + window)`.
/// An in-phase stretch ends where that fraction reaches the onset threshold
/// (or the lower `exit + 1/window` inside a repeat prior window), localized at
/// the first defective position of the window; a masked stretch ends where the
/// fraction drops to the exit threshold. Decisions freeze once fewer than half
/// a window remains, and a callable island shorter than one window is absorbed
/// by its masked neighbours.
pub(super) fn segment(defects: &[Defect], prior: &[bool], thresholds: &Thresholds) -> Vec<bool> {
    let positions = defects.len();
    let window = thresholds.window.max(1);
    let mut cumulative = vec![0_usize; positions + 1];
    for (index, defect) in defects.iter().enumerate() {
        cumulative[index + 1] = cumulative[index] + usize::from(*defect != Defect::None);
    }
    let fraction = |index: usize| {
        let end = index.saturating_add(window).min(positions);
        (cumulative[end] - cumulative[index]) as f64 / (end - index) as f64
    };
    let prior_onset = thresholds.exit + 1.0 / window as f64;
    let freeze_below = window.div_ceil(2);

    let mut callable = vec![false; positions];
    let mut in_phase = positions > 0 && fraction(0) < thresholds.onset;
    let mut index = 0;
    while index < positions {
        if positions - index < freeze_below {
            callable[index..].fill(in_phase);
            break;
        }
        if in_phase {
            let threshold = if prior[index] {
                prior_onset
            } else {
                thresholds.onset
            };
            if fraction(index) >= threshold {
                let end = index.saturating_add(window).min(positions);
                let onset = (index..end)
                    .find(|&position| defects[position] != Defect::None)
                    .unwrap_or(index);
                callable[index..onset].fill(true);
                in_phase = false;
                index = onset;
                continue;
            }
            callable[index] = true;
        } else if fraction(index) <= thresholds.exit {
            in_phase = true;
            callable[index] = true;
        }
        index += 1;
    }
    absorb_islands(&mut callable, window);
    callable
}

/// Masks every callable run shorter than one window.
fn absorb_islands(callable: &mut [bool], window: usize) {
    let mut start = 0;
    while start < callable.len() {
        if !callable[start] {
            start += 1;
            continue;
        }
        let length = callable[start..].iter().take_while(|&&flag| flag).count();
        if length < window {
            callable[start..start + length].fill(false);
        }
        start += length;
    }
}

#[cfg(test)]
mod tests {
    use crate::model::callability::RepeatUnit;

    use super::*;

    const THRESHOLDS: Thresholds = Thresholds {
        window: 4,
        onset: 0.5,
        exit: 0.0,
        double_ratio: 0.33,
    };

    fn pattern(text: &str) -> Vec<Defect> {
        text.chars()
            .map(|symbol| match symbol {
                '.' => Defect::None,
                'd' => Defect::Double,
                'w' => Defect::Weak,
                's' => Defect::Spacing,
                other => panic!("unsupported defect symbol {other}"),
            })
            .collect()
    }

    fn render(callable: &[bool]) -> String {
        callable
            .iter()
            .map(|&flag| if flag { 'c' } else { 'm' })
            .collect()
    }

    fn run(text: &str, thresholds: &Thresholds) -> String {
        let defects = pattern(text);
        render(&segment(&defects, &vec![false; defects.len()], thresholds))
    }

    #[test]
    fn keeps_a_clean_read_callable_and_masks_a_uniformly_bad_one() {
        assert_eq!(run("............", &THRESHOLDS), "cccccccccccc");
        assert_eq!(run("dd.dd.dd.dd.", &THRESHOLDS), "mmmmmmmmmmmm");
        assert_eq!(run("ww.ss.dd.ws.", &THRESHOLDS), "mmmmmmmmmmmm");
    }

    #[test]
    fn localizes_the_onset_at_the_first_defect_and_exits_on_a_clean_window() {
        assert_eq!(
            run("......dd.dd.............", &THRESHOLDS),
            "ccccccmmmmmccccccccccccc"
        );
    }

    #[test]
    fn masks_a_defective_tail_and_recovers_after_a_transient_stretch() {
        assert_eq!(run("..........dd", &THRESHOLDS), "ccccccccccmm");
        assert_eq!(
            run("........dddd....................", &THRESHOLDS),
            "ccccccccmmmmcccccccccccccccccccc"
        );
    }

    #[test]
    fn lowers_the_onset_threshold_inside_a_repeat_prior_window() {
        let defects = pattern("........d...d...........");
        assert_eq!(
            render(&segment(&defects, &vec![false; defects.len()], &THRESHOLDS)),
            "cccccccccccccccccccccccc"
        );
        let mut prior = vec![false; defects.len()];
        prior[7..11].fill(true);
        assert_eq!(
            render(&segment(&defects, &prior, &THRESHOLDS)),
            "ccccccccmmmmmccccccccccc"
        );
    }

    #[test]
    fn freezes_the_last_half_window_in_the_current_state() {
        let wide = Thresholds {
            window: 6,
            ..THRESHOLDS
        };
        assert_eq!(run("...........d", &wide), "cccccccccccc");
        assert_eq!(run("dddddddddd.d", &wide), "mmmmmmmmmmmm");
    }

    #[test]
    fn absorbs_callable_islands_shorter_than_one_window() {
        let mut callable = [
            false, true, true, true, false, true, true, true, true, false,
        ];
        absorb_islands(&mut callable, 4);
        assert_eq!(
            callable,
            [
                false, false, false, false, false, true, true, true, true, false
            ]
        );
        let mut leading = [true, true, false, false];
        absorb_islands(&mut leading, 3);
        assert_eq!(leading, [false; 4]);
    }

    #[test]
    fn classifies_defects_by_precedence() {
        let feature = |weak, ratio, spacing| PositionFeatures {
            primary_channel: Some(0),
            secondary_channel: Some(1),
            dominance: 0.5,
            secondary_ratio: ratio,
            spacing_deviation: spacing,
            weak,
        };
        let features = [
            feature(true, 0.9, 0.9),
            feature(false, 0.33, 0.9),
            feature(false, 0.2, 0.5),
            feature(false, 0.2, 0.4),
            feature(false, 0.1, 0.0),
        ];
        assert_eq!(
            defects(&features, &THRESHOLDS),
            [
                Defect::Weak,
                Defect::Double,
                Defect::Spacing,
                Defect::None,
                Defect::None
            ]
        );
    }

    #[test]
    fn prior_windows_start_after_each_run_and_never_cover_a_run() {
        let runs = [RepeatRun {
            call_start_0based: 2,
            call_end_0based_exclusive: 10,
            unit: RepeatUnit::Homopolymer(1),
        }];
        let prior = prior_windows(14, &runs, 4);
        let expected = (0..14)
            .map(|index| (10..14).contains(&index))
            .collect::<Vec<_>>();
        assert_eq!(prior, expected);
        assert_eq!(prior_windows(10, &runs, 4), [false; 10]);
        let adjacent = [
            runs[0],
            RepeatRun {
                call_start_0based: 11,
                call_end_0based_exclusive: 13,
                unit: RepeatUnit::Homopolymer(2),
            },
        ];
        let prior = prior_windows(16, &adjacent, 4);
        let expected = (0..16)
            .map(|index| index == 10 || (13..16).contains(&index))
            .collect::<Vec<_>>();
        assert_eq!(prior, expected);
    }
}
