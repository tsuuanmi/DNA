//! The completion event of each reference-free stage: what the stage made
//! of the read, and the settings it ran with. Logging only; nothing reads
//! these counts.

use std::time::Instant;

use super::SangerConfig;
use crate::model::basecalls::{BaseCalls, PeakSource};
use crate::model::callability::{PhaseState, ReadCallability};
use crate::model::quality::QualityControlResult;
use crate::model::signal::SignalAnalysis;

/// Logs `basecalling_completed`.
pub(super) fn basecalling(calls: &BaseCalls, config: &SangerConfig, stage_started: Instant) {
    let canonical_primary = calls
        .calls
        .iter()
        .filter(|call| call.primary != 'N')
        .count();
    let unresolved_primary = calls.len() - canonical_primary;
    let two_channel_iupac = calls
        .calls
        .iter()
        .filter(|call| call.qualifying_channels.len() == 2 && call.ambiguity != 'N')
        .count();
    let multi_channel_unresolved = calls
        .calls
        .iter()
        .filter(|call| call.qualifying_channels.len() > 2 && call.ambiguity == 'N')
        .count();
    let calls_with_fallback = calls
        .calls
        .iter()
        .filter(|call| {
            call.peaks
                .iter()
                .any(|peak| peak.source == PeakSource::LocusFallback)
        })
        .count();
    let vendor_compared = calls
        .calls
        .iter()
        .filter(|call| call.vendor_agrees.is_some())
        .count();
    let vendor_disagreements = calls
        .calls
        .iter()
        .filter(|call| call.vendor_agrees == Some(false))
        .count();
    tracing::info!(
        event = "basecalling_completed",
        elapsed_ms = stage_started.elapsed().as_millis(),
        calls = calls.len(),
        canonical_primary,
        unresolved_primary,
        two_channel_iupac,
        multi_channel_unresolved,
        calls_with_ploc_fallback = calls_with_fallback,
        vendor_compared,
        vendor_disagreements,
        secondary_peak_ratio = %format_args!("{:.4}", config.basecalling.secondary_peak_ratio),
    );
}

/// Logs `signal_processing_completed`.
pub(super) fn signal_processing(
    signal: &SignalAnalysis,
    config: &SangerConfig,
    stage_started: Instant,
) {
    let maximum_secondary_snr = signal
        .windows
        .iter()
        .map(|window| window.maximum_secondary_snr)
        .fold(0.0_f64, f64::max);
    let profiled_loci = signal
        .loci
        .iter()
        .filter(|locus| locus.profile.is_some())
        .count();
    let locus_vendor_length_mismatches = signal.integrity.vendor_length_mismatch_count();
    let clipped_channel_samples = signal.integrity.clipped_channel_samples;
    tracing::info!(
        event = "signal_processing_completed",
        elapsed_ms = stage_started.elapsed().as_millis(),
        loci = signal.loci.len(),
        profiled_loci,
        windows = signal.windows.len(),
        noisy_windows = signal.noisy_window_count(),
        noisy_regions = signal.noisy_regions.len(),
        noisy_calls = signal.noisy_call_count(),
        spikes = signal.spikes.len(),
        spike_calls = signal.spike_call_count(),
        window_size_bases = config.signal_processing.window_size_bases,
        minimum_noisy_windows = config.signal_processing.minimum_noisy_windows,
        minimum_primary_snr = %format_args!("{:.4}", config.signal_processing.minimum_primary_snr),
        maximum_secondary_snr = %format_args!("{maximum_secondary_snr:.4}"),
        ploc_vendor_length_mismatches = locus_vendor_length_mismatches,
        clipped_channel_samples,
        maximum_to_median_event_signal_ratio = ?signal.integrity.maximum_to_median_event_signal_ratio,
    );
}

/// Logs `callability_completed`.
pub(super) fn callability(
    calls: &BaseCalls,
    callability: &ReadCallability,
    config: &SangerConfig,
    stage_started: Instant,
) {
    let segment_map = callability
        .segments
        .iter()
        .map(|segment| {
            let shadow = segment.shadow.map_or_else(String::new, |shadow| {
                let offsets = if segment.state == PhaseState::Dephased {
                    let offsets = shadow
                        .offsets()
                        .map(|offset| format!("{offset:+}"))
                        .collect::<Vec<_>>();
                    format!("({})", offsets.join(","))
                } else {
                    String::new()
                };
                format!(
                    "{offsets}[{:.3},{:.3}]",
                    shadow.main_share, shadow.far_share
                )
            });
            format!(
                "{}..{}:{}{shadow}{}",
                segment.call_start_0based,
                segment.call_end_0based_exclusive,
                segment.state.label(),
                if segment.after_repeat { "+repeat" } else { "" },
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    let callable_fraction = callability.callable_count() as f64 / calls.len() as f64;
    tracing::info!(
        event = "callability_completed",
        elapsed_ms = stage_started.elapsed().as_millis(),
        calls = calls.len(),
        repeats = callability.repeats.len(),
        segments = callability.segments.len(),
        in_phase_segments = callability.segment_count(PhaseState::InPhase),
        dephased_segments = callability.segment_count(PhaseState::Dephased),
        mixed_segments = callability.segment_count(PhaseState::Mixed),
        weak_segments = callability.segment_count(PhaseState::Weak),
        irregular_segments = callability.segment_count(PhaseState::Irregular),
        masked_calls = callability.masked_count(),
        callable = %format_args!(
            "{}..{}",
            callability.callable_start_0based, callability.callable_end_0based_exclusive
        ),
        callable_fraction = %format_args!("{callable_fraction:.4}"),
        segment_map = %segment_map,
        window_calls = config.callability.window_calls,
        onset_defect_fraction = %format_args!("{:.4}", config.callability.onset_defect_fraction),
        exit_defect_fraction = %format_args!("{:.4}", config.callability.exit_defect_fraction),
        minimum_main_share = %format_args!("{:.4}", config.callability.minimum_main_share),
        maximum_far_share = %format_args!("{:.4}", config.callability.maximum_far_share),
        minimum_shadow_share = %format_args!("{:.4}", config.callability.minimum_shadow_share),
        weak_amplitude_fraction = %format_args!("{:.4}", config.callability.weak_amplitude_fraction),
        repeat_min_length = config.callability.repeat_min_length,
        minimum_callable_calls = config.callability.minimum_callable_calls,
    );
}

/// Logs `quality_control_completed`.
pub(super) fn quality_control(
    calls: &BaseCalls,
    quality: &QualityControlResult,
    stage_started: Instant,
) {
    let score_min = quality
        .per_call
        .iter()
        .map(|quality| quality.relative_quality_score)
        .min()
        .unwrap_or(0);
    let score_max = quality
        .per_call
        .iter()
        .map(|quality| quality.relative_quality_score)
        .max()
        .unwrap_or(0);
    let score_mean = if quality.per_call.is_empty() {
        0.0
    } else {
        quality
            .per_call
            .iter()
            .map(|quality| u64::from(quality.relative_quality_score))
            .sum::<u64>() as f64
            / quality.per_call.len() as f64
    };
    let max_penalty = quality
        .per_call
        .iter()
        .map(|quality| quality.penalty)
        .max()
        .unwrap_or(0);
    let vendor_quality_applicable = quality
        .per_call
        .iter()
        .filter(|quality| quality.vendor_quality_applies)
        .count();
    let trimmed_left = quality.trim_start_0based;
    let trimmed_right = calls
        .len()
        .saturating_sub(quality.trim_end_0based_exclusive);
    let retained_fraction = quality.retained_sequence.len() as f64 / calls.len() as f64;
    tracing::info!(
        event = "quality_control_completed",
        elapsed_ms = stage_started.elapsed().as_millis(),
        trim = %format_args!(
            "{}..{}",
            quality.trim_start_0based, quality.trim_end_0based_exclusive
        ),
        retained = quality.retained_sequence.len(),
        trimmed_left,
        trimmed_right,
        retained_fraction = %format_args!("{retained_fraction:.4}"),
        relative_score_min = score_min,
        relative_score_mean = %format_args!("{score_mean:.2}"),
        relative_score_max = score_max,
        max_penalty,
        vendor_quality_applicable,
    );
}
