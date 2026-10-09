//! Shared reference-free read processing for reusable scientific capabilities.

use std::time::Instant;

use crate::basecalling;
use crate::callability;
use crate::config::Config;
use crate::error::{CallabilityError, Result};
use crate::model::basecalls::{BaseCalls, PeakSource};
use crate::model::callability::{PhaseState, ReadCallability, ReadRejection};
use crate::model::quality::QualityControlResult;
use crate::model::sanger::Chromatogram;
use crate::model::signal::SignalAnalysis;
use crate::quality_control;
use crate::signal_processing;

/// Warning counts shared by command-level operational summaries.
pub(crate) struct ReadWarnings {
    pub(crate) unresolved_primary_calls: usize,
    pub(crate) multi_channel_unresolved_calls: usize,
    pub(crate) vendor_disagreements: usize,
    pub(crate) locus_vendor_length_mismatches: usize,
    pub(crate) clipped_channel_samples: usize,
}

/// Scientific read products shared by reference-free and reference-guided paths.
pub(crate) struct ProcessedRead {
    pub(crate) calls: BaseCalls,
    pub(crate) signal: SignalAnalysis,
    pub(crate) callability: ReadCallability,
    pub(crate) quality: QualityControlResult,
    pub(crate) warnings: ReadWarnings,
}

/// A read after the stages that decide whether it can be analyzed at all.
pub(crate) struct PreparedRead {
    pub(crate) calls: BaseCalls,
    pub(crate) signal: SignalAnalysis,
    pub(crate) callability: ReadCallability,
    pub(crate) warnings: ReadWarnings,
}

impl PreparedRead {
    /// Why the read cannot be analyzed, or `None` when it has enough callable
    /// calls.
    pub(crate) fn rejection(&self, config: &Config) -> Option<ReadRejection> {
        let callable_calls = self.callability.callable_count();
        let minimum_callable_calls = config.callability.minimum_callable_calls;
        (callable_calls < minimum_callable_calls).then_some(ReadRejection {
            callable_calls,
            minimum_callable_calls,
        })
    }
}

/// Runs the scientific stages that require no reference, emitting one
/// `tracing` stage span and completion event per stage; a read with too few
/// callable calls fails typed.
pub(crate) fn process(trace: &Chromatogram, config: &Config) -> Result<ProcessedRead> {
    let prepared = prepare(trace, config)?;
    finish(trace, prepared, config)
}

/// Runs basecalling, signal processing, and callability.
pub(crate) fn prepare(trace: &Chromatogram, config: &Config) -> Result<PreparedRead> {
    let stage = tracing::info_span!("basecalling").entered();
    let stage_started = Instant::now();
    let calls = basecalling::call(trace, &config.basecalling)?;
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

    drop(stage);
    let stage = tracing::info_span!("signal_processing").entered();
    let stage_started = Instant::now();
    let signal = signal_processing::analyze(trace, &calls, &config.signal_processing)?;
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
        window_size_bases = config.signal_processing.window_size_bases,
        minimum_noisy_windows = config.signal_processing.minimum_noisy_windows,
        minimum_primary_snr = %format_args!("{:.4}", config.signal_processing.minimum_primary_snr),
        maximum_secondary_snr = %format_args!("{maximum_secondary_snr:.4}"),
        ploc_vendor_length_mismatches = locus_vendor_length_mismatches,
        clipped_channel_samples,
        maximum_to_median_event_signal_ratio = ?signal.integrity.maximum_to_median_event_signal_ratio,
    );

    drop(stage);
    let stage = tracing::info_span!("callability").entered();
    let stage_started = Instant::now();
    let callability = callability::analyze(trace, &calls, &signal, config)?;
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

    drop(stage);
    Ok(PreparedRead {
        calls,
        signal,
        callability,
        warnings: ReadWarnings {
            unresolved_primary_calls: unresolved_primary,
            multi_channel_unresolved_calls: multi_channel_unresolved,
            vendor_disagreements,
            locus_vendor_length_mismatches,
            clipped_channel_samples,
        },
    })
}

/// Runs quality control on a prepared read; a read with too few callable
/// calls fails typed.
pub(crate) fn finish(
    trace: &Chromatogram,
    prepared: PreparedRead,
    config: &Config,
) -> Result<ProcessedRead> {
    if let Some(rejection) = prepared.rejection(config) {
        return Err(CallabilityError::TooFewCallableCalls {
            callable: rejection.callable_calls,
            minimum: rejection.minimum_callable_calls,
        }
        .into());
    }
    let PreparedRead {
        calls,
        signal,
        callability,
        warnings,
    } = prepared;
    let _stage = tracing::info_span!("quality_control").entered();
    let stage_started = Instant::now();
    let quality = quality_control::analyze(
        trace,
        &calls,
        &callability,
        &config.quality_control,
        config.variant_calling.read_end_margin,
    )?;
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

    Ok(ProcessedRead {
        calls,
        signal,
        callability,
        quality,
        warnings,
    })
}
