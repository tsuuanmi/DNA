//! Shared projection of merged observational signal-quality regions.

use crate::model::result::{
    DNAQualityResult, IntervalResult, NoisyRegionResult, TraceIntegrityResult,
};
use crate::model::signal::{DNAAnalysis, TraceIntegrity};

/// Projects merged noisy regions while omitting internal rolling windows.
pub(super) fn project(signal: DNAAnalysis) -> DNAQualityResult {
    let integrity = project_integrity(&signal.integrity);
    let noisy_regions = signal
        .noisy_regions
        .into_iter()
        .map(|region| NoisyRegionResult {
            calls: IntervalResult {
                start: region.call_start_0based,
                end: region.call_end_0based_exclusive,
            },
            samples: IntervalResult {
                start: region.sample_start_0based,
                end: region.sample_end_0based_exclusive,
            },
            minimum_primary_snr: region.minimum_primary_snr,
        })
        .collect();
    DNAQualityResult {
        integrity,
        noisy_regions,
    }
}

/// Projects immutable trace-integrity evidence for reuse by sample read summaries.
pub(super) fn project_integrity(integrity: &TraceIntegrity) -> TraceIntegrityResult {
    TraceIntegrityResult {
        ploc_count: integrity.ploc_count,
        vendor_primary_count: integrity.vendor_primary_count,
        vendor_quality_count: integrity.vendor_quality_count,
        minimum_ploc_spacing: integrity.minimum_ploc_spacing,
        median_ploc_spacing: integrity.median_ploc_spacing,
        maximum_ploc_spacing: integrity.maximum_ploc_spacing,
        clipped_channel_samples: integrity.clipped_channel_samples,
        maximum_to_median_event_signal_ratio: integrity.maximum_to_median_event_signal_ratio,
    }
}
