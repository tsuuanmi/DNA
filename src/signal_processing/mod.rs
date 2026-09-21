//! Observation-only rolling signal-quality analysis.

mod features;
mod integrity;
mod locus_evidence;
mod regions;
mod statistics;

use crate::config::DNAProcessingConfig;
use crate::error::Result;
use crate::model::basecalls::BaseCalls;
use crate::model::signal::DNAAnalysis;
use crate::model::trace::Chromatogram;

/// Calculates rolling SNR features, basecall-independent locus evidence, and merged noisy regions.
pub(crate) fn analyze(
    trace: &Chromatogram,
    calls: &BaseCalls,
    config: &DNAProcessingConfig,
) -> Result<DNAAnalysis> {
    let windows = features::calculate(trace, calls, config)?;
    let loci = locus_evidence::calculate(trace, config)?;
    let noisy_regions = regions::merge(&windows, config.minimum_noisy_windows);
    let integrity = integrity::assess(trace, &loci)?;
    Ok(DNAAnalysis {
        integrity,
        loci,
        windows,
        noisy_regions,
    })
}
