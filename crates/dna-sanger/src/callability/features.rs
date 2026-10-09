//! Per-position features of the callability core.
//!
//! Everything here is derived from plain [`PositionEvidence`] records in call
//! order; nothing refers to the chromatogram, the reference, or a profile.

use crate::model::callability::{PositionEvidence, PositionFeatures};
use crate::signal_processing::round_metric;
use dna_kernel::error::{CallabilityError, Result};

/// Half-width, in calls, of the window whose median spacing anchors the
/// spacing deviation of one position.
const SPACING_HALF_WINDOW: usize = 8;
/// Smallest position count the spacing feature handles.
pub(super) const MINIMUM_POSITIONS: usize = 2;

/// Derives one feature record per position.
///
/// A position is weak when it carries no positive signal or its primary
/// amplitude is below `weak_amplitude_fraction` of the read's median primary
/// amplitude; the read-wide median keeps homopolymers, whose local noise
/// estimate is inflated, from looking weak.
pub(super) fn calculate(
    evidence: &[PositionEvidence],
    weak_amplitude_fraction: f64,
) -> Result<Vec<PositionFeatures>> {
    if evidence.len() < MINIMUM_POSITIONS {
        return Err(CallabilityError::TooFewPositions {
            positions: evidence.len(),
            minimum: MINIMUM_POSITIONS,
        }
        .into());
    }
    let ranked = evidence
        .iter()
        .map(|position| rank(position.amplitudes))
        .collect::<Vec<_>>();
    let spacings = evidence
        .windows(2)
        .map(|pair| pair[1].coordinate.saturating_sub(pair[0].coordinate) as f64)
        .collect::<Vec<_>>();
    let weak_below = weak_amplitude_fraction * median_primary_amplitude(evidence, &ranked);
    Ok(evidence
        .iter()
        .enumerate()
        .map(|(index, position)| {
            let (primary_channel, secondary_channel) = ranked[index];
            let total = position.amplitudes.iter().sum::<f64>();
            let (dominance, secondary_ratio, weak) = match (primary_channel, secondary_channel) {
                (Some(primary), Some(secondary)) => (
                    round_metric(
                        (position.amplitudes[primary] - position.amplitudes[secondary]) / total,
                    ),
                    round_metric(position.amplitudes[secondary] / position.amplitudes[primary]),
                    position.amplitudes[primary] < weak_below,
                ),
                _ => (0.0, 0.0, true),
            };
            PositionFeatures {
                primary_channel,
                secondary_channel,
                dominance,
                secondary_ratio,
                spacing_deviation: spacing_deviation(index, &spacings),
                weak,
            }
        })
        .collect())
}

/// Median primary amplitude over the positions with positive signal, or zero
/// when there is none.
fn median_primary_amplitude(
    evidence: &[PositionEvidence],
    ranked: &[(Option<usize>, Option<usize>)],
) -> f64 {
    let mut amplitudes = evidence
        .iter()
        .zip(ranked)
        .filter_map(|(position, &(primary, _))| primary.map(|channel| position.amplitudes[channel]))
        .collect::<Vec<_>>();
    if amplitudes.is_empty() {
        return 0.0;
    }
    amplitudes.sort_by(f64::total_cmp);
    let middle = amplitudes.len() / 2;
    if amplitudes.len().is_multiple_of(2) {
        f64::midpoint(amplitudes[middle - 1], amplitudes[middle])
    } else {
        amplitudes[middle]
    }
}

/// Ranks the two strongest channels; ties favour the lower channel index and
/// a position without positive signal has no primary channel.
fn rank(amplitudes: [f64; 4]) -> (Option<usize>, Option<usize>) {
    let mut order = [0_usize, 1, 2, 3];
    order.sort_by(|left, right| {
        amplitudes[*right]
            .total_cmp(&amplitudes[*left])
            .then_with(|| left.cmp(right))
    });
    if amplitudes[order[0]] > 0.0 {
        (Some(order[0]), Some(order[1]))
    } else {
        (None, None)
    }
}

/// Relative deviation of the position's neighbouring spacings from the median
/// spacing of the surrounding window; zero when the median carries no
/// information.
fn spacing_deviation(index: usize, spacings: &[f64]) -> f64 {
    let window_start = index.saturating_sub(SPACING_HALF_WINDOW);
    let window_end = index
        .saturating_add(SPACING_HALF_WINDOW)
        .min(spacings.len());
    let mut window = spacings[window_start..window_end].to_vec();
    window.sort_by(f64::total_cmp);
    let middle = window.len() / 2;
    let median = if window.len().is_multiple_of(2) {
        f64::midpoint(window[middle - 1], window[middle])
    } else {
        window[middle]
    };
    if median <= 0.0 {
        return 0.0;
    }
    let deviation = [index.checked_sub(1), Some(index)]
        .into_iter()
        .flatten()
        .filter_map(|spacing| spacings.get(spacing))
        .map(|spacing| (spacing - median).abs() / median)
        .fold(0.0_f64, f64::max);
    round_metric(deviation)
}

#[cfg(test)]
mod tests {
    use dna_kernel::error::Error;

    use super::*;

    fn position(amplitudes: [f64; 4], coordinate: usize) -> PositionEvidence {
        PositionEvidence {
            amplitudes,
            coordinate,
            primary: None,
        }
    }

    #[test]
    fn ranks_channels_with_lower_index_winning_ties() {
        assert_eq!(rank([5.0, 9.0, 9.0, 1.0]), (Some(1), Some(2)));
        assert_eq!(rank([0.0, 0.0, 0.0, 0.0]), (None, None));
        assert_eq!(rank([0.0, 3.0, 0.0, 0.0]), (Some(1), Some(0)));
    }

    #[test]
    fn computes_dominance_ratio_and_weakness() -> Result<()> {
        let features = calculate(
            &[
                position([1_000.0, 30.0, 20.0, 10.0], 10),
                position([0.0, 500.0, 500.0, 0.0], 22),
                position([0.0, 0.0, 0.0, 0.0], 34),
                position([40.0, 0.0, 0.0, 0.0], 46),
            ],
            0.1,
        )?;
        assert_eq!(features[0].dominance, 0.915_094);
        assert_eq!(features[0].secondary_ratio, 0.03);
        assert!(!features[0].weak);
        assert_eq!(features[1].dominance, 0.0);
        assert_eq!(features[1].secondary_ratio, 1.0);
        assert!(!features[1].weak);
        assert_eq!(features[2].primary_channel, None);
        assert!(features[2].weak);
        // Median primary amplitude is 500; 40 is below one tenth of it.
        assert!(features[3].weak);
        Ok(())
    }

    #[test]
    fn measures_spacing_deviation_against_the_local_median() -> Result<()> {
        let coordinates = [0, 12, 24, 36, 48, 60, 90, 102, 114];
        let evidence = coordinates
            .iter()
            .map(|&coordinate| position([1_000.0, 0.0, 0.0, 0.0], coordinate))
            .collect::<Vec<_>>();
        let features = calculate(&evidence, 0.1)?;
        assert_eq!(features[0].spacing_deviation, 0.0);
        assert_eq!(features[5].spacing_deviation, 1.5);
        assert_eq!(features[6].spacing_deviation, 1.5);
        assert_eq!(features[7].spacing_deviation, 0.0);
        Ok(())
    }

    #[test]
    fn rejects_fewer_than_two_positions() {
        assert!(matches!(
            calculate(&[position([1.0, 0.0, 0.0, 0.0], 0)], 0.1),
            Err(Error::Callability(CallabilityError::TooFewPositions {
                positions: 1,
                minimum: 2
            }))
        ));
    }
}
