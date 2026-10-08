//! Single-read analysis orchestration, reporting, and publication.

use std::path::Path;
use std::time::Instant;

use crate::error::Result;
use crate::input::sanger;
use crate::operation_log::OperationLog;
use crate::pipeline::path;
use crate::report::{self, CompletedAnalysis};
use crate::variant_analysis;

use super::Operation;

/// Runs one complete AB1-to-JSON analysis with one per-trace append-only log.
pub(crate) fn run(trace: &Path, reference: &Path, config_path: &Path) -> Result<()> {
    let trace_stem = path::trace_stem(trace)?;
    Operation::begin(trace_stem, "analysis_failed", || {
        tracing::info!(
            event = "analysis_started",
            version = env!("CARGO_PKG_VERSION"),
            trace_path = ?trace.display().to_string(),
            reference_path = ?reference.display().to_string(),
        );
    })?
    .run(|log, started| analyze(trace, reference, config_path, log, started))
}

fn analyze(
    trace: &Path,
    reference: &Path,
    config_path: &Path,
    log: &OperationLog,
    analysis_started: Instant,
) -> Result<()> {
    let stage = tracing::info_span!("input_loading").entered();
    let stage_started = Instant::now();
    let prepared = sanger::prepare_analysis(trace, reference, config_path)?;
    let output = path::analysis_output(trace)?;
    let inputs = sanger::load_analysis(prepared)?;
    tracing::info!(
        event = "inputs_loaded",
        elapsed_ms = stage_started.elapsed().as_millis(),
        trace_name = ?inputs.trace.source_name,
        trace_sha256 = %inputs.trace.source_sha256,
        samples = inputs.trace.sample_count(),
        call_loci = inputs.trace.call_count(),
        vendor_primary = inputs.trace.vendor.primary.is_some(),
        vendor_quality = inputs.trace.vendor.qualities.is_some(),
        reference_name = ?inputs.reference.name,
        reference_sha256 = %inputs.reference.sequence_sha256,
        topology = ?inputs.reference.topology,
        reference_bases = inputs.reference.len(),
        config_path = ?inputs.config.source_path.display().to_string(),
        config_sha256 = %inputs.config.source_sha256,
        output_path = ?output.display().to_string(),
    );

    drop(stage);
    let completed =
        variant_analysis::observation::build(&inputs.trace, &inputs.reference, &inputs.config)?;

    let stage = tracing::info_span!("reporting").entered();
    let stage_started = Instant::now();
    let warning_total = completed.warning_total;
    let result = report::build_analysis(CompletedAnalysis {
        reference: inputs.reference,
        read: completed.read,
    })?;
    let result_variants = result.variants.len();
    let schema_version = result.schema_version;
    let bytes = report::serialize(&result)?;

    drop(stage);
    let _stage = tracing::info_span!("result_publication").entered();
    tracing::info!(
        event = "result_ready_for_publication",
        elapsed_ms = stage_started.elapsed().as_millis(),
        total_elapsed_ms = analysis_started.elapsed().as_millis(),
        schema = schema_version,
        variants = result_variants,
        warnings = warning_total,
        output_path = ?output.display().to_string(),
        bytes = bytes.len(),
    );
    log.sync()?;
    report::publish(&output, &bytes)
}
