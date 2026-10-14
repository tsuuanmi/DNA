//! Observation-only rolling signal-quality analysis.

mod config;

pub use config::{RawSignalProcessingConfig, SignalProcessingConfig};

mod features;
mod integrity;
mod locus_evidence;
mod regions;
mod statistics;

pub(crate) use statistics::round_metric;

use crate::model::basecalls::BaseCalls;
use crate::model::sanger::Chromatogram;
use crate::model::signal::SignalAnalysis;
use dna_kernel::error::Result;

/// Calculates rolling SNR features, basecall-independent locus evidence, and merged noisy regions.
///
/// Both feature windows and locus contexts estimate channel statistics over
/// runs of whole locus windows; each distinct span is estimated once.
pub(crate) fn analyze(
    trace: &Chromatogram,
    calls: &BaseCalls,
    config: &SignalProcessingConfig,
) -> Result<SignalAnalysis> {
    let mut estimates = statistics::SpanEstimates::new(trace);
    let windows = features::calculate(calls, config, &mut estimates)?;
    let loci = locus_evidence::calculate(config, &mut estimates)?;
    let noisy_regions = regions::merge(&windows, config.minimum_noisy_windows);
    let integrity = integrity::assess(trace, &loci)?;
    Ok(SignalAnalysis {
        integrity,
        loci,
        windows,
        noisy_regions,
    })
}
