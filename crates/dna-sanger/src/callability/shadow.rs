//! Shadow model of the callability core.
//!
//! Over a stretch of calls, the normalized amplitudes are fitted as a
//! non-negative mix of one-hot templates of the read's own primary calls
//! shifted by `-3..=+3` calls: offset zero is the main ladder, the others are
//! slippage shadows. Superimposed ladders offset by one call are explained by
//! the main ladder and near shadows; mixed signal is not. The fitted weights
//! are model coefficients, never mixture or heteroplasmy fractions.

use std::ops::Range;

use crate::model::callability::{
    PositionEvidence, PositionFeatures, SHADOW_OFFSETS, ShadowSummary,
};
use crate::signal_processing::round_metric;

/// Number of model offsets; slot `k + 3` holds offset `k`.
const SLOTS: usize = 7;
/// Model offsets per slot.
const OFFSETS: [isize; SLOTS] = [-3, -2, -1, 0, 1, 2, 3];
/// Slot of the main ladder.
const MAIN: usize = 3;
/// Shadow slots in the fixed search order `+1, -1, +2, -2, +3, -3`.
const SEARCH_ORDER: [usize; SLOTS - 1] = [4, 2, 5, 1, 6, 0];
/// Fewest usable positions a fit needs: one more than the model's parameters.
pub(super) const MINIMUM_FIT_POSITIONS: usize = SLOTS + 1;
/// Relative pivot below which a subset of templates is treated as singular.
const SINGULAR_PIVOT: f64 = 1e-9;
/// Relative error-sum improvement a later subset needs to replace the best.
const MINIMUM_IMPROVEMENT: f64 = 1e-9;

/// Classification thresholds of the shadow model.
#[derive(Debug, Clone, Copy)]
pub(super) struct Rules {
    pub(super) main_minimum: f64,
    pub(super) far_maximum: f64,
    pub(super) shadow_minimum: f64,
}

/// Non-negative template weights of one stretch, rounded to six decimals;
/// `weights[k + 3]` belongs to offset `k`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct ShadowFit {
    pub(super) weights: [f64; SLOTS],
}

/// Fits the shadow model over `calls`.
///
/// A position contributes when it is not weak and its primary call is
/// resolved; a neighbouring template outside the read or on an unresolved call
/// contributes nothing. Returns `None` with fewer than
/// [`MINIMUM_FIT_POSITIONS`] usable positions or when every weight is zero.
pub(super) fn fit(
    evidence: &[PositionEvidence],
    features: &[PositionFeatures],
    calls: Range<usize>,
) -> Option<ShadowFit> {
    let mut gram = [[0.0; SLOTS]; SLOTS];
    let mut cross = [0.0; SLOTS];
    let mut squares = 0.0;
    let mut positions = 0_usize;
    for index in calls {
        let position = &evidence[index];
        let total = position.amplitudes.iter().sum::<f64>();
        if features[index].weak || position.primary.is_none() || total <= 0.0 {
            continue;
        }
        let profile = position.amplitudes.map(|amplitude| amplitude / total);
        let templates = OFFSETS.map(|offset| {
            index
                .checked_add_signed(offset)
                .and_then(|neighbour| evidence.get(neighbour))
                .and_then(|neighbour| neighbour.primary)
        });
        for (slot, template) in templates.iter().enumerate() {
            let Some(channel) = *template else {
                continue;
            };
            cross[slot] += profile[channel];
            for (other, other_template) in templates.iter().enumerate() {
                if *other_template == Some(channel) {
                    gram[slot][other] += 1.0;
                }
            }
        }
        squares += profile.iter().map(|value| value * value).sum::<f64>();
        positions += 1;
    }
    if positions < MINIMUM_FIT_POSITIONS {
        return None;
    }
    let weights = best_weights(&gram, &cross, squares)?.map(round_metric);
    (weights.iter().sum::<f64>() > 0.0).then_some(ShadowFit { weights })
}

/// Exact non-negative least squares by enumerating every subset of shadow
/// templates next to the main ladder, fewest templates first in the fixed
/// search order; a later subset wins only with a strictly smaller error sum.
fn best_weights(
    gram: &[[f64; SLOTS]; SLOTS],
    cross: &[f64; SLOTS],
    squares: f64,
) -> Option<[f64; SLOTS]> {
    let mut subsets = (0_u8..1 << SEARCH_ORDER.len()).collect::<Vec<_>>();
    subsets.sort_by_key(|subset| (subset.count_ones(), *subset));
    let tolerance = MINIMUM_IMPROVEMENT * squares.max(1.0);
    let mut best: Option<(f64, [f64; SLOTS])> = None;
    for subset in subsets {
        let mut active = vec![MAIN];
        active.extend(
            SEARCH_ORDER
                .iter()
                .enumerate()
                .filter(|(bit, _)| subset & (1 << bit) != 0)
                .map(|(_, &slot)| slot),
        );
        let Some(weights) = solve(gram, cross, &active) else {
            continue;
        };
        if weights.iter().any(|&weight| weight < 0.0) {
            continue;
        }
        let error = error_sum(gram, cross, squares, &weights);
        if best.is_none_or(|(best_error, _)| error < best_error - tolerance) {
            best = Some((error, weights));
        }
    }
    best.map(|(_, weights)| weights)
}

/// Least-squares weights of the `active` templates by Cholesky factorization
/// of their Gram submatrix, or `None` when the submatrix is singular.
fn solve(
    gram: &[[f64; SLOTS]; SLOTS],
    cross: &[f64; SLOTS],
    active: &[usize],
) -> Option<[f64; SLOTS]> {
    let size = active.len();
    let scale = active
        .iter()
        .map(|&slot| gram[slot][slot])
        .fold(0.0, f64::max);
    let mut lower = vec![vec![0.0; size]; size];
    for row in 0..size {
        for column in 0..=row {
            let value = gram[active[row]][active[column]]
                - dot(&lower[row][..column], &lower[column][..column]);
            if row == column {
                if value <= SINGULAR_PIVOT * scale {
                    return None;
                }
                lower[row][row] = value.sqrt();
            } else {
                lower[row][column] = value / lower[column][column];
            }
        }
    }
    let mut forward = vec![0.0; size];
    for row in 0..size {
        forward[row] =
            (cross[active[row]] - dot(&lower[row][..row], &forward[..row])) / lower[row][row];
    }
    let mut solution = vec![0.0; size];
    for row in (0..size).rev() {
        let later = (row + 1..size)
            .map(|inner| lower[inner][row] * solution[inner])
            .sum::<f64>();
        solution[row] = (forward[row] - later) / lower[row][row];
    }
    let mut weights = [0.0; SLOTS];
    for (&slot, value) in active.iter().zip(solution) {
        weights[slot] = value;
    }
    Some(weights)
}

/// Dot product of two equally long slices.
fn dot(left: &[f64], right: &[f64]) -> f64 {
    left.iter()
        .zip(right)
        .map(|(left, right)| left * right)
        .sum()
}

/// Residual sum of squares `y·y − 2 w·b + wᵀ G w`, clamped at zero.
fn error_sum(
    gram: &[[f64; SLOTS]; SLOTS],
    cross: &[f64; SLOTS],
    squares: f64,
    weights: &[f64; SLOTS],
) -> f64 {
    let mut error = squares;
    for slot in 0..SLOTS {
        error -= 2.0 * weights[slot] * cross[slot];
        for other in 0..SLOTS {
            error += weights[slot] * gram[slot][other] * weights[other];
        }
    }
    error.max(0.0)
}

/// Summarizes a fit and decides whether it explains the stretch as dephased:
/// the main share reaches `minimum_main_share`, the far share stays at or
/// below `maximum_far_share`, and a one-call shadow reaches
/// `minimum_shadow_share`.
pub(super) fn summarize(fit: &ShadowFit, rules: &Rules) -> (bool, ShadowSummary) {
    let total = fit.weights.iter().sum::<f64>();
    let shares = fit.weights.map(|weight| weight / total);
    let main_share = round_metric(shares[MAIN]);
    let far_share = round_metric(shares[0] + shares[1] + shares[5] + shares[6]);
    let offsets =
        SHADOW_OFFSETS.map(|offset| round_metric(shares[slot(offset)]) >= rules.shadow_minimum);
    let near = SHADOW_OFFSETS
        .iter()
        .zip(offsets)
        .any(|(&offset, reported)| offset.abs() == 1 && reported);
    let dephased = near && main_share >= rules.main_minimum && far_share <= rules.far_maximum;
    (
        dephased,
        ShadowSummary {
            main_share,
            far_share,
            offsets,
        },
    )
}

/// Slot of a shadow offset.
fn slot(offset: i8) -> usize {
    OFFSETS
        .iter()
        .position(|&candidate| candidate == isize::from(offset))
        .unwrap_or(MAIN)
}

#[cfg(test)]
mod tests {
    use super::*;

    const RULES: Rules = Rules {
        main_minimum: 0.35,
        far_maximum: 0.12,
        shadow_minimum: 0.1,
    };

    /// A non-repetitive primary sequence over channels A, C, G, T.
    const SEQUENCE: &str = "ACGTTGCAAGCTCATGGACTTACGGATCCAGTGCATAGCTGACCTAGTCGAT";

    fn channel(base: char) -> usize {
        match base {
            'A' => 0,
            'C' => 1,
            'G' => 2,
            'T' => 3,
            other => panic!("unsupported synthetic base {other}"),
        }
    }

    /// Evidence whose amplitudes mix the primary ladder with `shadows`, each an
    /// offset and the weight of the shifted copy.
    fn mixture(sequence: &str, main: f64, shadows: &[(isize, f64)]) -> Vec<PositionEvidence> {
        let channels = sequence.chars().map(channel).collect::<Vec<_>>();
        (0..channels.len())
            .map(|index| {
                let mut amplitudes = [0.0; 4];
                amplitudes[channels[index]] += 1_000.0 * main;
                for &(offset, weight) in shadows {
                    if let Some(&shadow) = index
                        .checked_add_signed(offset)
                        .and_then(|neighbour| channels.get(neighbour))
                    {
                        amplitudes[shadow] += 1_000.0 * weight;
                    }
                }
                PositionEvidence {
                    amplitudes,
                    coordinate: index * 12,
                    primary: Some(channels[index]),
                }
            })
            .collect()
    }

    fn features(evidence: &[PositionEvidence]) -> Vec<PositionFeatures> {
        evidence
            .iter()
            .map(|_| PositionFeatures {
                primary_channel: Some(0),
                secondary_channel: Some(1),
                dominance: 0.5,
                secondary_ratio: 0.5,
                spacing_deviation: 0.0,
                weak: false,
            })
            .collect()
    }

    fn fitted(evidence: &[PositionEvidence], calls: Range<usize>) -> Option<ShadowFit> {
        fit(evidence, &features(evidence), calls)
    }

    #[test]
    fn recovers_a_two_sided_shadow_ladder() {
        let evidence = mixture(SEQUENCE, 0.5, &[(-1, 0.25), (1, 0.25)]);
        // Interior positions only, so every template is defined.
        let Some(fit) = fitted(&evidence, 3..evidence.len() - 3) else {
            panic!("a two-sided ladder must be fitted");
        };
        assert_eq!(fit.weights, [0.0, 0.0, 0.25, 0.5, 0.25, 0.0, 0.0]);
        let (dephased, summary) = summarize(&fit, &RULES);
        assert!(dephased);
        assert_eq!(summary.main_share, 0.5);
        assert_eq!(summary.far_share, 0.0);
        assert_eq!(summary.offsets().collect::<Vec<_>>(), [-1, 1]);
    }

    #[test]
    fn recovers_a_one_sided_shadow_ladder() {
        let evidence = mixture(SEQUENCE, 0.7, &[(-1, 0.3)]);
        let fit = fitted(&evidence, 3..evidence.len() - 3);
        assert_eq!(
            fit.map(|fit| fit.weights),
            Some([0.0, 0.0, 0.3, 0.7, 0.0, 0.0, 0.0])
        );
    }

    #[test]
    fn never_reports_negative_weights() {
        let mut evidence = mixture(SEQUENCE, 1.0, &[]);
        // Signal on channels unrelated to any neighbour pulls no template up.
        for position in &mut evidence {
            let unrelated = (position.primary.unwrap_or(0) + 2) % 4;
            position.amplitudes[unrelated] += 300.0;
        }
        let fit = fitted(&evidence, 0..evidence.len());
        assert!(
            fit.is_some_and(|fit| fit.weights.iter().all(|&weight| weight >= 0.0)),
            "{fit:?}"
        );
    }

    #[test]
    fn treats_coincident_templates_as_one() {
        // Inside a dinucleotide repeat the offset-two template equals the main
        // ladder; the smaller subset wins and the duplicate stays at zero.
        let evidence = mixture(&"AC".repeat(20), 1.0, &[]);
        let fit = fitted(&evidence, 3..evidence.len() - 3);
        assert_eq!(
            fit.map(|fit| fit.weights),
            Some([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0])
        );
    }

    #[test]
    fn is_deterministic() {
        let evidence = mixture(SEQUENCE, 0.5, &[(-1, 0.2), (1, 0.2), (2, 0.1)]);
        assert_eq!(
            fitted(&evidence, 0..evidence.len()),
            fitted(&evidence, 0..evidence.len())
        );
    }

    #[test]
    fn needs_enough_usable_positions() {
        let mut evidence = mixture(SEQUENCE, 1.0, &[]);
        assert_eq!(fitted(&evidence, 0..MINIMUM_FIT_POSITIONS - 1), None);
        assert!(fitted(&evidence, 0..MINIMUM_FIT_POSITIONS).is_some());
        evidence[3].primary = None;
        assert_eq!(fitted(&evidence, 0..MINIMUM_FIT_POSITIONS), None);
        let mut weak = features(&evidence);
        weak[4].weak = true;
        assert_eq!(fit(&evidence, &weak, 4..MINIMUM_FIT_POSITIONS + 4), None);
        assert!(fit(&evidence, &weak, 4..MINIMUM_FIT_POSITIONS + 5).is_some());
    }

    #[test]
    fn fits_at_the_read_edges() {
        let evidence = mixture(SEQUENCE, 0.6, &[(-1, 0.4)]);
        let fit = fitted(&evidence, 0..evidence.len());
        assert!(
            fit.is_some_and(|fit| fit.weights[MAIN] > 0.5 && fit.weights[2] > 0.3),
            "{fit:?}"
        );
    }

    fn summary(weights: [f64; 7]) -> (bool, Vec<i8>) {
        let (dephased, summary) = summarize(&ShadowFit { weights }, &RULES);
        (dephased, summary.offsets().collect())
    }

    #[test]
    fn separates_dephased_from_mixed_stretches() {
        // Moderate two-sided slippage.
        assert_eq!(
            summary([0.0, 0.03, 0.17, 0.5, 0.17, 0.04, 0.0]),
            (true, vec![-1, 1])
        );
        // Several length populations of similar size after a long run.
        assert!(!summary([0.02, 0.04, 0.16, 0.24, 0.18, 0.03, 0.04]).0);
        // Globally mixed signal with shadows at every offset.
        assert!(!summary([0.04, 0.08, 0.11, 0.39, 0.08, 0.05, 0.07]).0);
        // Double peaks no neighbour explains.
        assert_eq!(
            summary([0.0, 0.0, 0.0, 0.7, 0.0, 0.0, 0.0]),
            (false, vec![])
        );
        // Dinucleotide slippage gives shadows two calls away only.
        assert_eq!(
            summary([0.0, 0.3, 0.0, 0.4, 0.0, 0.3, 0.0]),
            (false, vec![-2, 2])
        );
    }

    #[test]
    fn applies_the_thresholds_inclusively() {
        // main 0.35, far 0.12, near shadow 0.1 exactly.
        let (dephased, summary) = summarize(
            &ShadowFit {
                weights: [0.06, 0.0, 0.1, 0.35, 0.43, 0.0, 0.06],
            },
            &RULES,
        );
        assert_eq!(summary.main_share, 0.35);
        assert_eq!(summary.far_share, 0.12);
        assert!(dephased);
    }
}
