//! Shared basecalling, signal analysis, quality control, and stage logging.

use std::time::Instant;

use crate::basecalling;
use crate::config::Config;
use crate::error::Result;
use crate::logger::Logger;
use crate::model::basecalls::{BaseCalls, PeakSource};
use crate::model::quality::QualityControlResult;
use crate::model::signal::DNAAnalysis;
use crate::model::trace::Chromatogram;
use crate::quality_control;
use crate::signal_processing;

/// Warning counts shared by command-level operational summaries.
pub(crate) struct ReadWarnings {
    pub(crate) unresolved_primary_calls: usize,
    pub(crate) multi_channel_unresolved_calls: usize,
    pub(crate) vendor_disagreements: usize,
    pub(crate) ploc_vendor_length_mismatches: usize,
    pub(crate) clipped_channel_samples: usize,
}

/// Scientific read products shared by reference-free and reference-guided paths.
pub(crate) struct ProcessedRead {
    pub(crate) calls: BaseCalls,
    pub(crate) signal: DNAAnalysis,
    pub(crate) quality: QualityControlResult,
    pub(crate) warnings: ReadWarnings,
}

/// Runs and logs the scientific stages that require no reference.
pub(crate) fn process(
    trace: &Chromatogram,
    config: &Config,
    logger: &mut Logger,
    stage: &mut &'static str,
) -> Result<ProcessedRead> {
    *stage = "basecalling";
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
                .any(|peak| peak.source == PeakSource::PlocFallback)
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
    logger.info(
        module_path!(),
        line!(),
        format_args!(
            concat!(
                "event=basecalling_completed elapsed_ms={} calls={} canonical_primary={} ",
                "unresolved_primary={} two_channel_iupac={} multi_channel_unresolved={} ",
                "calls_with_ploc_fallback={} vendor_compared={} vendor_disagreements={} ",
                "secondary_peak_ratio={:.4}"
            ),
            stage_started.elapsed().as_millis(),
            calls.len(),
            canonical_primary,
            unresolved_primary,
            two_channel_iupac,
            multi_channel_unresolved,
            calls_with_fallback,
            vendor_compared,
            vendor_disagreements,
            config.basecalling.secondary_peak_ratio
        ),
    )?;

    *stage = "signal_processing";
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
    let ploc_vendor_length_mismatches = signal.integrity.vendor_length_mismatch_count();
    let clipped_channel_samples = signal.integrity.clipped_channel_samples;
    logger.info(
        module_path!(),
        line!(),
        format_args!(
            concat!(
                "event=signal_processing_completed elapsed_ms={} loci={} profiled_loci={} ",
                "windows={} noisy_windows={} noisy_regions={} noisy_calls={} ",
                "window_size_bases={} minimum_noisy_windows={} minimum_primary_snr={:.4} ",
                "maximum_secondary_snr={:.4} ploc_vendor_length_mismatches={} ",
                "clipped_channel_samples={} maximum_to_median_event_signal_ratio={:?}"
            ),
            stage_started.elapsed().as_millis(),
            signal.loci.len(),
            profiled_loci,
            signal.windows.len(),
            signal.noisy_window_count(),
            signal.noisy_regions.len(),
            signal.noisy_call_count(),
            config.signal_processing.window_size_bases,
            config.signal_processing.minimum_noisy_windows,
            config.signal_processing.minimum_primary_snr,
            maximum_secondary_snr,
            ploc_vendor_length_mismatches,
            clipped_channel_samples,
            signal.integrity.maximum_to_median_event_signal_ratio
        ),
    )?;

    *stage = "quality_control";
    let stage_started = Instant::now();
    let quality = quality_control::analyze(trace, &calls, &config.quality_control)?;
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
    logger.info(
        module_path!(),
        line!(),
        format_args!(
            concat!(
                "event=quality_control_completed elapsed_ms={} trim={}..{} retained={} ",
                "trimmed_left={} trimmed_right={} retained_fraction={:.4} ",
                "relative_score_min={} relative_score_mean={:.2} relative_score_max={} ",
                "max_penalty={} vendor_quality_applicable={}"
            ),
            stage_started.elapsed().as_millis(),
            quality.trim_start_0based,
            quality.trim_end_0based_exclusive,
            quality.retained_sequence.len(),
            trimmed_left,
            trimmed_right,
            retained_fraction,
            score_min,
            score_mean,
            score_max,
            max_penalty,
            vendor_quality_applicable
        ),
    )?;

    Ok(ProcessedRead {
        calls,
        signal,
        quality,
        warnings: ReadWarnings {
            unresolved_primary_calls: unresolved_primary,
            multi_channel_unresolved_calls: multi_channel_unresolved,
            vendor_disagreements,
            ploc_vendor_length_mismatches,
            clipped_channel_samples,
        },
    })
}
