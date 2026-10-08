//! Reference-free base re-calling, quality trimming, and result publication.

use std::path::Path;
use std::time::Instant;

use crate::error::Result;
use crate::input::sanger;
use crate::operation_log::OperationLog;
use crate::pipeline::path;
use crate::read_processing::{self, ProcessedRead};
use crate::report::{self, CompletedBasecall};

use super::Operation;

/// Runs one complete AB1-to-basecalls JSON operation.
pub(crate) fn run(trace: &Path, config_path: &Path) -> Result<()> {
    let trace_stem = path::trace_stem(trace)?;
    Operation::begin(trace_stem, "basecall_failed", || {
        tracing::info!(
            event = "basecall_started",
            version = env!("CARGO_PKG_VERSION"),
            trace_path = ?trace.display().to_string(),
        );
    })?
    .run(|log, started| basecall(trace, config_path, log, started))
}

fn basecall(trace: &Path, config_path: &Path, log: &OperationLog, started: Instant) -> Result<()> {
    let stage = tracing::info_span!("input_loading").entered();
    let stage_started = Instant::now();
    let prepared = sanger::prepare_basecall(trace, config_path)?;
    let output = path::basecall_output(trace)?;
    let inputs = sanger::load_basecall(prepared)?;
    tracing::info!(
        event = "basecall_inputs_loaded",
        elapsed_ms = stage_started.elapsed().as_millis(),
        trace_name = ?inputs.trace.source_name,
        trace_sha256 = %inputs.trace.source_sha256,
        samples = inputs.trace.sample_count(),
        call_loci = inputs.trace.call_count(),
        vendor_primary = inputs.trace.vendor.primary.is_some(),
        vendor_quality = inputs.trace.vendor.qualities.is_some(),
        config_path = ?inputs.config.source_path.display().to_string(),
        config_sha256 = %inputs.config.source_sha256,
        output_path = ?output.display().to_string(),
    );

    drop(stage);
    let ProcessedRead {
        calls,
        signal,
        quality,
        warnings,
    } = read_processing::process(&inputs.trace, &inputs.config)?;
    let warning_total = warnings.unresolved_primary_calls
        + warnings.multi_channel_unresolved_calls
        + warnings.vendor_disagreements
        + warnings.locus_vendor_length_mismatches
        + warnings.clipped_channel_samples;
    if warning_total > 0 {
        tracing::warn!(
            event = "basecall_warning_summary",
            total = warning_total,
            unresolved_primary_calls = warnings.unresolved_primary_calls,
            multi_channel_unresolved_calls = warnings.multi_channel_unresolved_calls,
            vendor_disagreements = warnings.vendor_disagreements,
            ploc_vendor_length_mismatches = warnings.locus_vendor_length_mismatches,
            clipped_channel_samples = warnings.clipped_channel_samples,
        );
    }

    let stage = tracing::info_span!("reporting").entered();
    let stage_started = Instant::now();
    let result = report::build_basecall(CompletedBasecall {
        config: inputs.config,
        trace: inputs.trace,
        calls,
        signal,
        quality,
    })?;
    let schema_version = result.schema_version;
    let call_count = result.read.call_count;
    let retained = result.read.retained.len();
    let bytes = report::serialize(&result)?;

    drop(stage);
    let _stage = tracing::info_span!("result_publication").entered();
    tracing::info!(
        event = "basecall_ready_for_publication",
        elapsed_ms = stage_started.elapsed().as_millis(),
        total_elapsed_ms = started.elapsed().as_millis(),
        schema = schema_version,
        calls = call_count,
        retained,
        warnings = warning_total,
        output_path = ?output.display().to_string(),
        bytes = bytes.len(),
    );
    log.sync()?;
    report::publish(&output, &bytes)
}
