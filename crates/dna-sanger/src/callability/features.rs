//! Per-position features of the callability core.
//!
//! Everything here is derived from plain [`PositionEvidence`] records in call
//! order; nothing refers to the chromatogram, the reference, or a profile.

use crate::error::{CallabilityError, Result};
use crate::model::callability::{PositionEvidence, PositionFeatures};
use crate::signal_processing::round_metric;

/// Half-width, in calls, of the window whose median spacing anchors the
/// spacing deviation of one position.
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
    spacing_half_window: usize,
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
    let spacings = spacings(evidence);
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
                spacing_deviation: spacing_deviation(index, &spacings, spacing_half_window),
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

/// The spacings between consecutive positions, in call order.
fn spacings(evidence: &[PositionEvidence]) -> Vec<f64> {
    evidence
        .windows(2)
        .map(|pair| pair[1].coordinate.saturating_sub(pair[0].coordinate) as f64)
        .collect()
}

/// The median spacing of the window around position `index`: up to
/// `half_window` spacings on each side; zero for an empty window.
fn local_median(index: usize, spacings: &[f64], half_window: usize) -> f64 {
    let window_start = index.saturating_sub(half_window);
    let window_end = index.saturating_add(half_window).min(spacings.len());
    let mut window = spacings
        .get(window_start..window_end)
        .unwrap_or_default()
        .to_vec();
    if window.is_empty() {
        return 0.0;
    }
    window.sort_by(f64::total_cmp);
    let middle = window.len() / 2;
    if window.len().is_multiple_of(2) {
        f64::midpoint(window[middle - 1], window[middle])
    } else {
        window[middle]
    }
}

/// The junctions whose spacing marks a missing or an extra peak, each by the
/// position before it (SRS-CALL-026). With `M` the median spacing of that
/// position's window, a junction marks a missing peak when its spacing
/// exceeds `M`, and its span together with either neighbouring junction
/// exceeds `2M`, each by at least `spacing_defect` of `M`: the measure and
/// level of a spacing defect (SRS-CALL-002). It marks an extra peak in the
/// mirror case: its spacing falls short of `M`, and its span together with
/// either neighbouring junction falls short of `2M`, each by at least
/// `spacing_defect` of `M`. A peak the basecaller missed leaves a gap about
/// twice the read's spacing, and the ladder stays a call behind; a peak it
/// called twice leaves two calls a fraction of the spacing apart, and the
/// ladder a call ahead. A peak merely placed late or early is followed by a
/// spacing that makes up for it.
pub(super) fn gaps(
    evidence: &[PositionEvidence],
    half_window: usize,
    spacing_defect: f64,
) -> Vec<usize> {
    let spacings = spacings(evidence);
    (0..spacings.len())
        .filter(|&index| {
            let median = local_median(index, &spacings, half_window);
            if median <= 0.0 {
                return false;
            }
            let spacing = spacings[index];
            let neighbours = [index.checked_sub(1), index.checked_add(1)]
                .into_iter()
                .flatten()
                .filter_map(|neighbour| spacings.get(neighbour))
                .copied()
                .collect::<Vec<_>>();
            // `sign` 1 asks for a span wider than expected, -1 for narrower.
            let defect = |sign: f64| {
                let off = |span: f64, expected: f64| {
                    round_metric(sign * (span - expected * median) / median) >= spacing_defect
                };
                off(spacing, 1.0)
                    && neighbours
                        .iter()
                        .all(|&neighbour| off(spacing + neighbour, 2.0))
            };
            defect(1.0) || defect(-1.0)
        })
        .collect()
}

/// Relative deviation of the position's neighbouring spacings from the median
/// spacing of the surrounding window; zero when the median carries no
/// information.
fn spacing_deviation(index: usize, spacings: &[f64], half_window: usize) -> f64 {
    let median = local_median(index, spacings, half_window);
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
    use crate::error::SangerError as Error;

    use super::*;

    fn position(amplitudes: [f64; 4], coordinate: usize) -> PositionEvidence {
        PositionEvidence {
            amplitudes,
            coordinate,
            primary: None,
        }
    }

    // SRS-CALL-002
    #[test]
    fn ranks_channels_with_lower_index_winning_ties() {
        assert_eq!(rank([5.0, 9.0, 9.0, 1.0]), (Some(1), Some(2)));
        assert_eq!(rank([0.0, 0.0, 0.0, 0.0]), (None, None));
        assert_eq!(rank([0.0, 3.0, 0.0, 0.0]), (Some(1), Some(0)));
    }

    // SRS-CALL-002
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
            8,
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

    // SRS-CALL-002
    #[test]
    fn measures_spacing_deviation_against_the_local_median() -> Result<()> {
        let coordinates = [0, 12, 24, 36, 48, 60, 90, 102, 114];
        let evidence = coordinates
            .iter()
            .map(|&coordinate| position([1_000.0, 0.0, 0.0, 0.0], coordinate))
            .collect::<Vec<_>>();
        let features = calculate(&evidence, 0.1, 8)?;
        assert_eq!(features[0].spacing_deviation, 0.0);
        assert_eq!(features[5].spacing_deviation, 1.5);
        assert_eq!(features[6].spacing_deviation, 1.5);
        assert_eq!(features[7].spacing_deviation, 0.0);
        Ok(())
    }

    // SRS-CALL-026
    #[test]
    fn a_gap_is_a_spacing_at_least_the_defect_off_the_local_median() {
        let evidence = |coordinates: &[usize]| {
            coordinates
                .iter()
                .map(|&coordinate| position([1_000.0, 0.0, 0.0, 0.0], coordinate))
                .collect::<Vec<_>>()
        };
        // Spacings 12, 12, 12, 12, 12, 24, 12, 12, 12: the junction after
        // position 5 is twice the median, a missing peak.
        let missing = evidence(&[0, 12, 24, 36, 48, 60, 84, 96, 108, 120]);
        assert_eq!(gaps(&missing, 8, 0.5), [5]);
        assert_eq!(gaps(&missing, 8, 1.0), [5]);
        assert!(gaps(&missing, 8, 1.5).is_empty());
        // Spacings 12, 12, 4, 12, 12, 12, 12, 12, 12: two calls 4 scans apart
        // on one peak, an extra peak, the ladder a call ahead after it.
        let extra = evidence(&[0, 12, 24, 28, 40, 52, 64, 76, 88, 100]);
        assert_eq!(gaps(&extra, 8, 0.5), [2]);
        assert!(gaps(&extra, 8, 1.0).is_empty());
        // Spacings 12, 12, 12, 12, 12, 18, 6, 12, 12: a peak placed late,
        // whose next spacing makes up for it, marks neither.
        let late = evidence(&[0, 12, 24, 36, 48, 60, 78, 84, 96, 108]);
        assert!(gaps(&late, 8, 0.5).is_empty());
        // Spacings 12, 12, 12, 12, 12, 6, 18, 12, 12: nor does a peak placed
        // early.
        let early = evidence(&[0, 12, 24, 36, 48, 60, 66, 84, 96, 108]);
        assert!(gaps(&early, 8, 0.5).is_empty());
        // A constant spacing has none.
        let even = (0..10)
            .map(|index| position([1_000.0, 0.0, 0.0, 0.0], index * 12))
            .collect::<Vec<_>>();
        assert!(gaps(&even, 8, 0.5).is_empty());
    }

    #[test]
    fn rejects_fewer_than_two_positions() {
        assert!(matches!(
            calculate(&[position([1.0, 0.0, 0.0, 0.0], 0)], 0.1, 8),
            Err(Error::Callability(CallabilityError::TooFewPositions {
                positions: 1,
                minimum: 2
            }))
        ));
    }
}
