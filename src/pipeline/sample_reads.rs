//! Shared multi-trace read processing for sample operations.

use crate::config::Config;
use crate::error::Result;
use crate::model::read_observation::ReadObservation;
use crate::model::reference::Reference;
use crate::model::sanger::Chromatogram;
use crate::profile::Profile;
use crate::variant_analysis;

pub(crate) struct CompletedSampleReads {
    pub(crate) reads: Vec<ReadObservation>,
    pub(crate) warning_total: usize,
}

pub(crate) fn build(
    traces: &[Chromatogram],
    reference: &Reference,
    config: &Config,
    profile: &Profile,
) -> Result<CompletedSampleReads> {
    let mut reads = Vec::with_capacity(traces.len());
    let mut warning_total = 0usize;

    for (index, trace) in traces.iter().enumerate() {
        tracing::info!(
            event = "sample_read_started",
            read_index = index,
            trace_name = ?trace.source_name,
            trace_sha256 = %trace.source_sha256,
        );
        let completed = variant_analysis::observation::build(trace, reference, config, profile)?;
        warning_total += completed.warning_total;
        tracing::info!(
            event = "sample_read_completed",
            read_index = index,
            trace_sha256 = %completed.read.input_sha256,
            orientation = ?completed.read.alignment.orientation,
            segments = completed.read.alignment.reference_segments.len(),
            variants = completed.read.variants.reported.len(),
        );
        reads.push(completed.read);
    }

    Ok(CompletedSampleReads {
        reads,
        warning_total,
    })
}
