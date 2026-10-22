//! Electrophoretic spikes (SRS-SIG-012): a narrow event that rises in all
//! four analysed channels at once, far above the read's local signal.
//!
//! A spike is not DNA. Every dye shows it, each at the scan its mobility
//! correction moves it to, so the four apexes lie within about one call
//! spacing of each other. A call whose event falls inside it reads the
//! spike's channel, not the template's base. This module only finds the
//! spikes; the evidence adapter withholds the calls (SRS-CALL-027).

use crate::median;
use crate::model::basecalls::{BaseCall, BaseCalls};
use crate::model::sanger::Chromatogram;
use crate::model::signal::SignalSpike;

use super::SpikeParams;

/// The read's spikes, in sample order: within `span` median spacings, every
/// channel has a local maximum at least `height` times the local level, the
/// median primary peak of the `level_calls` calls on each side of the call
/// whose window holds it (that call left out). A spike's samples run from
/// the earliest apex's half-height start to the latest apex's half-height
/// end, and its calls are those whose primary peak, the highest of their
/// selected peaks, lies among them. None without a median spacing.
pub(super) fn detect(
    trace: &Chromatogram,
    calls: &BaseCalls,
    median_spacing: Option<f64>,
    parameters: &SpikeParams,
) -> Vec<SignalSpike> {
    let calls = &calls.calls;
    let Some(spacing) = median_spacing.filter(|spacing| *spacing > 0.0) else {
        return Vec::new();
    };
    if calls.is_empty() {
        return Vec::new();
    }
    let levels = levels(calls, parameters.level_calls);
    // The call whose window holds a sample: the last whose window starts at
    // or before it.
    let nearest = |sample: usize| {
        calls
            .partition_point(|call| call.window_start_0based <= sample)
            .saturating_sub(1)
    };
    let mut apexes = Vec::new();
    for (channel, values) in trace.channels.iter().enumerate() {
        for sample in 1..values.len().saturating_sub(1) {
            let value = values[sample];
            if values[sample - 1] < value
                && value >= values[sample + 1]
                && f64::from(value) >= parameters.height * levels[nearest(sample)]
            {
                apexes.push((sample, channel));
            }
        }
    }
    apexes.sort_unstable();
    // The span is a whole number of samples: a spacing of a few samples at
    // most a few times over.
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let span = (parameters.span * spacing).ceil() as usize;
    let mut spikes = Vec::new();
    let mut consumed = None;
    for &(start, _) in &apexes {
        if consumed.is_some_and(|last| start <= last) {
            continue;
        }
        let end = start + span;
        // Each channel's highest qualifying apex within the span.
        let mut highest: [Option<usize>; 4] = [None; 4];
        for &(sample, channel) in apexes
            .iter()
            .skip_while(|(sample, _)| *sample < start)
            .take_while(|(sample, _)| *sample <= end)
        {
            let values = &trace.channels[channel];
            if highest[channel].is_none_or(|best| values[sample] > values[best]) {
                highest[channel] = Some(sample);
            }
        }
        let Some(apex) = highest.iter().copied().collect::<Option<Vec<_>>>() else {
            continue;
        };
        let (first, last) =
            apex.iter()
                .enumerate()
                .fold((usize::MAX, 0), |(first, last), (channel, &sample)| {
                    let (left, right) = half_height(&trace.channels[channel], sample);
                    (first.min(left), last.max(right))
                });
        let samples = first..last + 1;
        let inside = calls
            .iter()
            .enumerate()
            .filter(|(_, call)| primary(call).is_some_and(|peak| samples.contains(&peak.0)))
            .map(|(index, _)| index);
        let calls_start = inside.clone().min().unwrap_or(0);
        let calls_end = inside.max().map_or(calls_start, |last| last + 1);
        spikes.push(SignalSpike {
            samples,
            calls: calls_start..calls_end,
        });
        consumed = apex.iter().copied().max();
    }
    spikes
}

/// A call's primary peak: the position and height of the highest of its
/// selected peaks.
fn primary(call: &BaseCall) -> Option<(usize, i32)> {
    call.peaks
        .iter()
        .max_by_key(|peak| peak.height)
        .map(|peak| (peak.position_0based, peak.height))
}

/// Each call's local level: the median primary peak of the `side` calls on
/// each side, the call left out, and at least one.
fn levels(calls: &[BaseCall], side: usize) -> Vec<f64> {
    let primaries = calls
        .iter()
        .map(|call| f64::from(primary(call).map_or(0, |peak| peak.1).max(0)))
        .collect::<Vec<_>>();
    (0..primaries.len())
        .map(|index| {
            let mut around = primaries[index.saturating_sub(side)..index]
                .iter()
                .chain(primaries.iter().skip(index + 1).take(side))
                .copied()
                .collect::<Vec<_>>();
            median::midpoint(&mut around).map_or(1.0, |median| median.max(1.0))
        })
        .collect()
}

/// The samples around an apex where the channel stays at or above half its
/// height, as the first and last of them.
fn half_height(values: &[i32], apex: usize) -> (usize, usize) {
    let half = f64::from(values[apex]) / 2.0;
    let at_least = |sample: usize| f64::from(values[sample]) >= half;
    let mut left = apex;
    while left > 0 && at_least(left - 1) {
        left -= 1;
    }
    let mut right = apex;
    while right + 1 < values.len() && at_least(right + 1) {
        right += 1;
    }
    (left, right)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::basecalls::{ChannelPeak, PeakSource};
    use crate::model::sanger::VendorEvidence;
    use dna_kernel::model::nucleotide::Nucleotide;

    const SPACING: usize = 12;

    fn parameters() -> SpikeParams {
        SpikeParams {
            height: 3.0,
            span: 1.5,
            level_calls: 6,
        }
    }

    /// A read of `count` calls, one peak of 500 per call in a rotating
    /// channel at the call's centre.
    fn read(count: usize) -> (Chromatogram, BaseCalls) {
        let samples = count * SPACING;
        let mut channels: [Vec<i32>; 4] = std::array::from_fn(|_| vec![0; samples]);
        let mut calls = Vec::new();
        for index in 0..count {
            let center = index * SPACING + SPACING / 2;
            let primary = index % 4;
            for (offset, value) in [500, 380, 260, 140, 20].into_iter().enumerate() {
                channels[primary][center + offset] = value;
                channels[primary][center - offset] = value;
            }
            let base = Nucleotide::ALL[primary];
            calls.push(BaseCall {
                index_0based: index,
                locus_position_0based: center,
                window_start_0based: index * SPACING,
                window_end_0based_exclusive: (index + 1) * SPACING,
                peaks: std::array::from_fn(|channel| ChannelPeak {
                    base: Nucleotide::ALL[channel],
                    height: if channel == primary { 500 } else { 0 },
                    position_0based: center,
                    source: PeakSource::LocalMaximum,
                }),
                primary_peak_evidence: None,
                primary: base.as_char(),
                ambiguity: base.as_char(),
                qualifying_channels: vec![base],
                vendor_agrees: None,
            });
        }
        let trace = Chromatogram {
            source_name: "synthetic.ab1".into(),
            source_sha256: String::new(),
            channels,
            locus_positions: calls
                .iter()
                .map(|call| call.locus_position_0based)
                .collect(),
            duplicate_loci: 0,
            vendor: VendorEvidence::default(),
        };
        let primary_sequence = calls.iter().map(|call| call.primary).collect();
        (
            trace,
            BaseCalls {
                calls,
                primary_sequence,
            },
        )
    }

    /// Adds a narrow peak of `height` to `channel` at `apex`.
    fn add(trace: &mut Chromatogram, channel: usize, apex: usize, height: i32) {
        trace.channels[channel][apex] += height;
        for offset in 1..4 {
            let value = height / (1 << offset);
            trace.channels[channel][apex + offset] += value;
            trace.channels[channel][apex - offset] += value;
        }
    }

    // SRS-SIG-012
    /// Four channels peaking within one spacing, each far above the local
    /// level, are one spike: its calls are those whose primary peak it
    /// covers. The apexes are staggered as mobility correction staggers them.
    #[test]
    fn finds_a_spike_in_all_four_channels() {
        let (mut trace, calls) = read(40);
        let center = 20 * SPACING + SPACING / 2;
        for (channel, offset) in [(3, 0), (0, 5), (1, 6), (2, 12)] {
            add(&mut trace, channel, center - 6 + offset, 4000);
        }
        let spikes = detect(&trace, &calls, Some(12.0), &parameters());
        assert_eq!(spikes.len(), 1);
        assert_eq!(spikes[0].calls, 20..21);
        assert!(spikes[0].samples.contains(&center));
    }

    // SRS-SIG-012
    /// Three tall channels, or four spread wider than the span, or four not
    /// tall enough against the local level, make no spike; nor does a read
    /// without a median spacing.
    #[test]
    fn ignores_three_channels_a_wide_spread_and_a_low_event() {
        let center = 20 * SPACING + SPACING / 2;
        let (mut three, calls) = read(40);
        for channel in 0..3 {
            add(&mut three, channel, center, 4000);
        }
        assert!(detect(&three, &calls, Some(12.0), &parameters()).is_empty());
        let (mut wide, _) = read(40);
        for (channel, offset) in [(0, 0), (1, 6), (2, 12), (3, 30)] {
            add(&mut wide, channel, center - 12 + offset, 4000);
        }
        assert!(detect(&wide, &calls, Some(12.0), &parameters()).is_empty());
        let (mut low, _) = read(40);
        for channel in 0..4 {
            add(&mut low, channel, center, 1000);
        }
        assert!(detect(&low, &calls, Some(12.0), &parameters()).is_empty());
        let (mut tall, _) = read(40);
        for channel in 0..4 {
            add(&mut tall, channel, center, 4000);
        }
        assert_eq!(detect(&tall, &calls, Some(12.0), &parameters()).len(), 1);
        assert!(detect(&tall, &calls, None, &parameters()).is_empty());
    }
}
