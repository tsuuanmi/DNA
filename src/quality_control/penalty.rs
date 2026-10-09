//! Per-call ambiguity and peak-spacing penalties.

use crate::error::{QualityControlError, Result};
use crate::model::basecalls::BaseCalls;
use crate::model::nucleotide::is_canonical;

/// Computes safe local ambiguity and spacing penalties per call.
pub(crate) fn calculate(calls: &BaseCalls, window_size: usize) -> Result<Vec<i32>> {
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
            .filter(|call| !is_canonical(call.ambiguity))
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

    Ok(penalties)
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
        let penalties = calculate(&calls(&[0, 4, 8, 20], &['A', 'A', 'N', 'A']), 3)?;
        assert_eq!(penalties, vec![3, 3, 5, 6]);
        Ok(())
    }
}
