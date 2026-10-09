//! Observation-only rolling signal-quality analysis.

mod features;
mod integrity;
mod locus_evidence;
mod regions;
mod statistics;

pub(crate) use statistics::round_metric;

use crate::config::SignalProcessingConfig;
use crate::error::Result;
use crate::model::basecalls::BaseCalls;
use crate::model::sanger::Chromatogram;
use crate::model::signal::SignalAnalysis;

/// Calculates rolling SNR features, basecall-independent locus evidence, and merged noisy regions.
pub(crate) fn analyze(
    trace: &Chromatogram,
    calls: &BaseCalls,
    config: &SignalProcessingConfig,
) -> Result<SignalAnalysis> {
    let windows = features::calculate(trace, calls, config)?;
    let loci = locus_evidence::calculate(trace, config)?;
    let noisy_regions = regions::merge(&windows, config.minimum_noisy_windows);
    let integrity = integrity::assess(trace, &loci)?;
    Ok(SignalAnalysis {
        integrity,
        loci,
        windows,
        noisy_regions,
    })
}
