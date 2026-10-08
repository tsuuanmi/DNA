//! Per-call ambiguity and peak-spacing penalties.

use crate::error::{QualityControlError, Result};
use crate::model::basecalls::BaseCalls;

/// Penalty vector and best contiguous section.
#[derive(Debug)]
pub(crate) struct PenaltyResult {
    pub(crate) penalties: Vec<i32>,
    pub(crate) best_start: usize,
    pub(crate) best_end: usize,
    pub(crate) best_average: f64,
}

/// Computes safe local penalties and the minimum-sum section.
pub(crate) fn calculate(
    calls: &BaseCalls,
    window_size: usize,
    best_fraction: f64,
) -> Result<PenaltyResult> {
    if calls.is_empty() || window_size == 0 {
        return Err(QualityControlError::EmptyPenaltyInput.into());
    }
    let count = calls.len();
    let mean_spacing = if count > 1 {
        calls
            .calls
            .windows(2)
            .map(|pair| pair[1].locus_position_0based - pair[0].locus_position_0based)
            .sum::<usize>() as f64
            / (count - 1) as f64
    } else {
        0.0
    };
    let half = window_size / 2;
    let mut penalties = Vec::with_capacity(count);
    for index in 0..count {
        let start = index.saturating_sub(half);
        let end = start.saturating_add(window_size).min(count);
        let ambiguity = calls.calls[start..end]
            .iter()
            .filter(|call| !matches!(call.ambiguity, 'A' | 'C' | 'G' | 'T'))
            .count();
        let mut distances = calls.calls[start..end]
            .windows(2)
            .map(|pair| pair[1].locus_position_0based - pair[0].locus_position_0based);
        let first = distances.next();
        let spacing_penalty = if let Some(first) = first {
            let (minimum, maximum) = distances.fold((first, first), |(minimum, maximum), value| {
                (minimum.min(value), maximum.max(value))
            });
            spacing_penalty(minimum, maximum, mean_spacing)
        } else {
            0
        };
        let ambiguity = i32::try_from(ambiguity)
            .map_err(|_| QualityControlError::Overflow("ambiguity penalty overflow"))?;
        penalties.push(ambiguity.saturating_add(spacing_penalty));
    }

    let best_length = best_section_length(count, best_fraction);
    let mut current: i64 = penalties[..best_length]
        .iter()
        .map(|value| i64::from(*value))
        .sum();
    let mut best_sum = current;
    let mut best_start = 0;
    for start in 1..=count - best_length {
        current -= i64::from(penalties[start - 1]);
        current += i64::from(penalties[start + best_length - 1]);
        if current < best_sum {
            best_sum = current;
            best_start = start;
        }
    }
    Ok(PenaltyResult {
        penalties,
        best_start,
        best_end: best_start + best_length,
        best_average: best_sum as f64 / best_length as f64,
    })
}

/// Floors the mean deviation of the extreme local peak spacings from the trace mean.
#[expect(
    clippy::cast_possible_truncation,
    reason = "spacings are bounded by the ABIF sample count, far below i32::MAX"
)]
fn spacing_penalty(minimum: usize, maximum: usize, mean_spacing: f64) -> i32 {
    f64::midpoint(
        (maximum as f64 - mean_spacing).abs(),
        (minimum as f64 - mean_spacing).abs(),
    )
    .floor() as i32
}

/// Length of the best section: the validated fraction of all calls, at least one.
#[expect(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "the configured fraction lies in (0, 1], so the floored product lies in [0, count]"
)]
fn best_section_length(count: usize, best_fraction: f64) -> usize {
    ((count as f64 * best_fraction).floor() as usize)
        .max(1)
        .min(count)
}

#[cfg(test)]
mod tests {
    use crate::model::basecalls::{BaseCall, ChannelPeak, PeakSource};
    use crate::model::nucleotide::Nucleotide;

    use super::*;

    fn calls(locus_positions: &[usize], ambiguities: &[char]) -> BaseCalls {
        let calls = locus_positions
            .iter()
            .zip(ambiguities)
            .enumerate()
            .map(|(index, (&locus_position, &ambiguity))| BaseCall {
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
                ambiguity,
                qualifying_channels: vec![Nucleotide::A],
                vendor_agrees: None,
            })
            .collect();
        BaseCalls {
            calls,
            primary_sequence: "A".repeat(locus_positions.len()),
        }
    }

    #[test]
    fn calculates_deterministic_ambiguity_and_spacing_penalties() -> Result<()> {
        let result = calculate(&calls(&[0, 4, 8, 20], &['A', 'A', 'N', 'A']), 3, 0.5)?;
        assert_eq!(result.penalties, vec![3, 3, 5, 6]);
        assert_eq!((result.best_start, result.best_end), (0, 2));
        assert_eq!(result.best_average, 3.0);
        Ok(())
    }
}
