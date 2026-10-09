//! Shared multi-trace read processing for sample operations.

use crate::config::Config;
use crate::error::Result;
use crate::model::read_observation::ReadObservation;
use crate::model::reference::Reference;
use crate::model::sample_evidence::RejectedSampleRead;
use crate::model::sanger::Chromatogram;
use crate::profile::Profile;
use crate::read_processing;
use crate::variant_analysis;

pub(crate) struct CompletedSampleReads {
    pub(crate) reads: Vec<ReadObservation>,
    pub(crate) rejected: Vec<RejectedSampleRead>,
    pub(crate) warning_total: usize,
}

/// Processes every trace independently; a read with too few callable calls is
/// recorded as rejected and the remaining reads continue.
pub(crate) fn build(
    traces: &[Chromatogram],
    reference: &Reference,
    config: &Config,
    profile: &Profile,
) -> Result<CompletedSampleReads> {
    let mut reads = Vec::with_capacity(traces.len());
    let mut rejected = Vec::new();
    let mut warning_total = 0usize;

    for (index, trace) in traces.iter().enumerate() {
        tracing::info!(
            event = "sample_read_started",
            read_index = index,
            trace_name = ?trace.source_name,
            trace_sha256 = %trace.source_sha256,
        );
        let prepared = read_processing::prepare(trace, config)?;
        if let Some(rejection) = prepared.rejection(config) {
            tracing::warn!(
                event = "sample_read_rejected",
                read_index = index,
                trace_sha256 = %trace.source_sha256,
                reason = "callable_calls_below_minimum",
                callable_calls = rejection.callable_calls,
                minimum_callable_calls = rejection.minimum_callable_calls,
            );
            rejected.push(RejectedSampleRead {
                input_name: trace.source_name.clone(),
                input_sha256: trace.source_sha256.clone(),
                integrity: prepared.signal.integrity,
                callability: prepared.callability,
                rejection,
            });
            continue;
        }
        let processed = read_processing::finish(trace, prepared, config)?;
        let completed =
            variant_analysis::observation::observe(trace, processed, reference, config, profile)?;
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
        rejected,
        warning_total,
    })
}
