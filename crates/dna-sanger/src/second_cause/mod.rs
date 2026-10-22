//! Why a call shows a second peak (SRS-CALL-021).
//!
//! A second peak has several causes a reviewer tells apart by its shape: the
//! tail of a neighbouring call of that base, a wide dye blob, a slippage
//! shadow of a neighbour's base, a peak barely above the noise, or a second
//! template. [`classify`] labels every call whose second base's share
//! exceeds `signal_processing.second_floor`, and every call masked as mixed,
//! from the read's own calls only: the primary channel's selected peak is the
//! call's centre, half the distance between its neighbours' primary peaks its
//! spacing, and the peak widths signal processing measures at the call's
//! primary event its widths. That width is the second channel's around the
//! primary event, not around its own apex; it is close for a centred peak,
//! the only one the width rule judges.
//!
//! The rules apply in order, the first that holds giving the cause:
//!
//! 1. `noise`: the second peak is below `second_cause_noise_ratio` times the
//!    call's lowest channel;
//! 2. `neighbour_tail`: the second channel has no local maximum in the call's
//!    window, or its peak lies at least `second_cause_tail_offset` spacings
//!    off the centre towards a neighbouring call of the second base;
//! 3. `off_centre`: the peak lies that far off the centre otherwise;
//! 4. `dye_blob`: the second channel is at least `second_cause_blob_width`
//!    times as wide as the median primary peak of the
//!    `second_cause_width_calls` calls on each side;
//! 5. `phase_shift`: the second base is a neighbouring call's base;
//! 6. `second_template` otherwise, which a run signal can show to be a shared
//!    `context` shadow (`run_signal::refine`).
//!
//! A call masked as mixed whose second share does not exceed the floor is a
//! `mixed_stretch`. The labels decide nothing: no verdict reads them.

use dna_kernel::model::nucleotide::Nucleotide;
use dna_kernel::read_evidence::SecondBase;

use crate::median;
use crate::model::basecalls::{BaseCall, BaseCalls, PeakSource};
use crate::model::callability::{PhaseState, ReadCallability};
use crate::model::locus_evidence::channel_width;
use crate::model::second_cause::{SecondCause, SecondPeak};
use crate::model::signal::SignalAnalysis;
use crate::signal_processing::SignalProcessingConfig;

/// The second base of a call's signal (SRS-CALL-015): the highest peak of the
/// other three channels, with its share of the primary and that peak
/// together. Whether it stands out is [`classify`]'s judgement against the
/// floor, set by the evidence adapter.
#[must_use]
pub fn second_base(call: &BaseCall) -> Option<SecondBase> {
    let primary = Nucleotide::from_char(call.primary)?.channel_index();
    let first = call.peaks.get(primary)?.height;
    let (channel, second) = call
        .peaks
        .iter()
        .enumerate()
        .filter(|(channel, _)| *channel != primary)
        .map(|(channel, peak)| (channel, peak.height))
        .max_by_key(|(channel, height)| (*height, std::cmp::Reverse(*channel)))?;
    let total = f64::from(first) + f64::from(second);
    (second > 0 && total > 0.0).then(|| SecondBase {
        base: Nucleotide::ALL[channel].as_char(),
        share: f64::from(second) / total,
        stands_out: false,
    })
}

/// The second peak of every call, by call index: `None` for a call whose
/// second share does not exceed the floor and that is not masked as mixed.
#[must_use]
pub fn classify(
    calls: &BaseCalls,
    signal: &SignalAnalysis,
    callability: &ReadCallability,
    config: &SignalProcessingConfig,
) -> Vec<Option<SecondPeak>> {
    (0..calls.calls.len())
        .map(|index| {
            let second = second_base(&calls.calls[index]);
            let mixed = callability.mask.get(index).copied().flatten() == Some(PhaseState::Mixed);
            match second {
                Some(second) if second.share > config.second_floor => Some(SecondPeak {
                    second: SecondBase {
                        base: second.base,
                        share: second.share,
                        stands_out: true,
                    },
                    cause: shape(calls, signal, index, second.base, config),
                }),
                _ if mixed => Some(SecondPeak {
                    second: SecondBase {
                        base: second.map_or('N', |second| second.base),
                        share: second.map_or(0.0, |second| second.share),
                        stands_out: false,
                    },
                    cause: SecondCause::MixedStretch,
                }),
                _ => None,
            }
        })
        .collect()
}

/// The cause of call `index`'s second peak of base `second` from its shape.
fn shape(
    calls: &BaseCalls,
    signal: &SignalAnalysis,
    index: usize,
    second: char,
    config: &SignalProcessingConfig,
) -> SecondCause {
    let parameters = &config.second_cause;
    let call = &calls.calls[index];
    let Some(other) = Nucleotide::from_char(second) else {
        return SecondCause::SecondTemplate;
    };
    let peak = call.peaks[other.channel_index()];
    let lowest = call
        .peaks
        .iter()
        .map(|peak| peak.height)
        .min()
        .unwrap_or(0)
        .max(1);
    if f64::from(peak.height) < parameters.noise_ratio * f64::from(lowest) {
        return SecondCause::Noise;
    }
    let neighbor = |offset: isize| {
        index
            .checked_add_signed(offset)
            .and_then(|neighbor| calls.calls.get(neighbor))
    };
    let (previous, next) = (neighbor(-1), neighbor(1));
    let center = |call: &BaseCall| {
        Nucleotide::from_char(call.primary).map_or(call.locus_position_0based, |base| {
            call.peaks[base.channel_index()].position_0based
        })
    };
    // Signed offset of the second peak from the call's centre, in spacings.
    let offset = previous.zip(next).and_then(|(previous, next)| {
        let spacing = (center(next) as f64 - center(previous) as f64) / 2.0;
        (spacing > 0.0).then(|| (peak.position_0based as f64 - center(call) as f64) / spacing)
    });
    let towards = offset.is_some_and(|offset| {
        (offset > 0.0 && next.is_some_and(|next| next.primary == second))
            || (offset < 0.0 && previous.is_some_and(|previous| previous.primary == second))
    });
    let shifted = offset.is_some_and(|offset| offset.abs() >= parameters.tail_offset);
    if peak.source != PeakSource::LocalMaximum || (shifted && towards) {
        return SecondCause::NeighborTail;
    }
    if shifted {
        return SecondCause::OffCenter;
    }
    if wide(calls, signal, index, other, config) {
        return SecondCause::DyeBlob;
    }
    if previous.is_some_and(|previous| previous.primary == second)
        || next.is_some_and(|next| next.primary == second)
    {
        return SecondCause::PhaseShift;
    }
    SecondCause::SecondTemplate
}

/// Whether call `index`'s `second` channel is at least the blob width ratio
/// times as wide as the median primary peak of its neighbouring calls.
fn wide(
    calls: &BaseCalls,
    signal: &SignalAnalysis,
    index: usize,
    second: Nucleotide,
    config: &SignalProcessingConfig,
) -> bool {
    let parameters = &config.second_cause;
    let width = |call: usize, channel: Nucleotide| {
        channel_width(&signal.loci, call, channel.channel_index())
    };
    let Some(own) = width(index, second) else {
        return false;
    };
    let first = index.saturating_sub(parameters.width_calls);
    let last = (index + parameters.width_calls).min(calls.calls.len().saturating_sub(1));
    let mut widths = (first..=last)
        .filter(|&call| call != index)
        .filter_map(|call| width(call, Nucleotide::from_char(calls.calls[call].primary)?))
        .collect::<Vec<_>>();
    median::upper(&mut widths, Ord::cmp)
        .is_some_and(|median| own as f64 >= parameters.blob_width * median as f64)
}

// SRS-CALL-021
#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::basecalls::ChannelPeak;
    use crate::model::locus_evidence::LocusEvidence;
    use crate::model::signal::SangerIntegrity;

    /// Ten samples between calls.
    const SPACING: usize = 10;

    fn config() -> SignalProcessingConfig {
        SignalProcessingConfig {
            second_floor: 0.20,
            ..crate::signal_processing::fixture(10)
        }
    }

    /// A read of `bases`, each call's primary peak 1000 at its centre and its
    /// other channels `others` there, every channel four samples wide.
    fn read(bases: &str, others: i32) -> (BaseCalls, SignalAnalysis, ReadCallability) {
        let calls = bases
            .chars()
            .enumerate()
            .map(|(index, base)| BaseCall {
                index_0based: index,
                locus_position_0based: index * SPACING,
                window_start_0based: index * SPACING,
                window_end_0based_exclusive: index * SPACING + SPACING / 2,
                peaks: std::array::from_fn(|channel| ChannelPeak {
                    base: Nucleotide::ALL[channel],
                    height: if Nucleotide::ALL[channel].as_char() == base {
                        1000
                    } else {
                        others
                    },
                    position_0based: index * SPACING,
                    source: PeakSource::LocalMaximum,
                }),
                primary_peak_evidence: None,
                primary: base,
                ambiguity: base,
                qualifying_channels: Vec::new(),
                vendor_agrees: None,
            })
            .collect::<Vec<_>>();
        let count = calls.len();
        let loci = (0..count)
            .map(|index| LocusEvidence {
                call_index_0based: index,
                locus_position_0based: index * SPACING,
                window_start_0based: index * SPACING,
                window_end_0based_exclusive: index * SPACING + 1,
                context_call_start_0based: index,
                context_call_end_0based_exclusive: index + 1,
                context_sample_start_0based: index * SPACING,
                context_sample_end_0based_exclusive: index * SPACING + 1,
                event_position_0based: index * SPACING,
                channel_heights: [0; 4],
                channel_baselines: [0.0; 4],
                channel_noise_sigmas: [1.0; 4],
                corrected_amplitudes: [0.0; 4],
                snrs: [0.0; 4],
                profile: None,
                channel_widths: [4; 4],
            })
            .collect();
        (
            BaseCalls {
                primary_sequence: bases.into(),
                calls,
            },
            SignalAnalysis {
                integrity: SangerIntegrity {
                    locus_count: count,
                    duplicate_locus_count: 0,
                    vendor_primary_count: None,
                    vendor_quality_count: None,
                    minimum_locus_spacing: None,
                    median_locus_spacing: None,
                    maximum_locus_spacing: None,
                    clipped_channel_samples: 0,
                    maximum_to_median_event_signal_ratio: None,
                },
                loci,
                windows: Vec::new(),
                noisy_regions: Vec::new(),
                spikes: Vec::new(),
            },
            ReadCallability::in_phase(count),
        )
    }

    /// The cause at call 2 of `ACAGT` (an A between a C and a G) once
    /// `change` set its second peak.
    fn cause_at_two(
        others: i32,
        change: impl FnOnce(&mut BaseCalls, &mut SignalAnalysis, &mut ReadCallability),
    ) -> Option<SecondCause> {
        let (mut calls, mut signal, mut callability) = read("ACAGT", others);
        change(&mut calls, &mut signal, &mut callability);
        classify(&calls, &signal, &callability, &config())[2].map(|peak| peak.cause)
    }

    /// Sets call 2's `base` channel to `height` at `offset` samples from its
    /// centre.
    fn second(calls: &mut BaseCalls, base: Nucleotide, height: i32, offset: isize) {
        let peak = &mut calls.calls[2].peaks[base.channel_index()];
        peak.height = height;
        peak.position_0based = (2 * SPACING).saturating_add_signed(offset);
    }

    #[test]
    fn a_clean_call_has_no_second_peak() {
        assert_eq!(cause_at_two(10, |_, _, _| {}), None);
    }

    #[test]
    fn a_second_peak_barely_above_the_lowest_channel_is_noise() {
        // 280 against 1000 is a share of 0.22, but under three times 100.
        let cause = cause_at_two(100, |calls, _, _| second(calls, Nucleotide::T, 280, 0));
        assert_eq!(cause, Some(SecondCause::Noise));
    }

    #[test]
    fn a_second_peak_without_a_maximum_or_towards_its_base_is_a_neighbor_tail() {
        let rising = cause_at_two(10, |calls, _, _| {
            second(calls, Nucleotide::T, 400, 0);
            calls.calls[2].peaks[Nucleotide::T.channel_index()].source = PeakSource::LocusFallback;
        });
        assert_eq!(rising, Some(SecondCause::NeighborTail));
        // Half a spacing towards the C before it.
        let shifted = cause_at_two(10, |calls, _, _| second(calls, Nucleotide::C, 400, -5));
        assert_eq!(shifted, Some(SecondCause::NeighborTail));
    }

    #[test]
    fn a_shifted_second_peak_towards_no_call_of_its_base_is_off_center() {
        // Half a spacing towards the G after it, a T.
        let cause = cause_at_two(10, |calls, _, _| second(calls, Nucleotide::T, 400, 5));
        assert_eq!(cause, Some(SecondCause::OffCenter));
    }

    #[test]
    fn a_wide_centerd_second_peak_is_a_dye_blob() {
        // Eight samples wide against the neighbours' four.
        let cause = cause_at_two(10, |calls, signal, _| {
            second(calls, Nucleotide::T, 400, 0);
            signal.loci[2].channel_widths[Nucleotide::T.channel_index()] = 8;
        });
        assert_eq!(cause, Some(SecondCause::DyeBlob));
    }

    #[test]
    fn a_centerd_second_peak_of_a_neighbors_base_is_a_phase_shift() {
        let cause = cause_at_two(10, |calls, _, _| second(calls, Nucleotide::G, 400, 0));
        assert_eq!(cause, Some(SecondCause::PhaseShift));
    }

    #[test]
    fn an_unexplained_centerd_second_peak_is_a_second_template() {
        let cause = cause_at_two(10, |calls, _, _| second(calls, Nucleotide::T, 400, 0));
        assert_eq!(cause, Some(SecondCause::SecondTemplate));
        // At the floor's share, 250 against 1000, there is no label.
        let floor = cause_at_two(10, |calls, _, _| second(calls, Nucleotide::T, 250, 0));
        assert_eq!(floor, None);
    }

    #[test]
    fn a_mixed_call_with_a_small_second_peak_is_a_mixed_stretch() {
        let cause = cause_at_two(10, |_, _, callability| {
            callability.mask[2] = Some(PhaseState::Mixed);
        });
        assert_eq!(cause, Some(SecondCause::MixedStretch));
    }
}
