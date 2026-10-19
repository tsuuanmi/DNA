//! Reference-free read evidence and result publication (PROP-0011): the
//! modality's per-call evidence exactly as the core receives it.

use std::path::Path;
use std::time::Instant;

use crate::error::Result;
use crate::input::sanger;
use crate::operation_log::OperationLog;
use crate::pipeline::path;
use crate::pipeline::plugins as plugin;
use crate::report::{self, CompletedEvidence};
use dna_sanger::read_processing::{self, ProcessedRead};

use super::Operation;

/// Runs one complete AB1-to-evidence JSON operation.
pub(crate) fn run(trace: &Path, config_path: &Path) -> Result<()> {
    let trace_stem = path::trace_stem(trace)?;
    Operation::begin(trace_stem, "evidence_failed", || {
        tracing::info!(
            event = "evidence_started",
            version = env!("CARGO_PKG_VERSION"),
            trace_path = ?trace.display().to_string(),
        );
    })?
    .run(|log, started| evidence(trace, config_path, log, started))
}

fn evidence(trace: &Path, config_path: &Path, log: &OperationLog, started: Instant) -> Result<()> {
    let stage = tracing::info_span!("input_loading").entered();
    let prepared = sanger::prepare_basecall(trace, config_path)?;
    let output = path::evidence_output(trace)?;
    let inputs = sanger::load_basecall(prepared)?;
    tracing::info!(
        event = "evidence_inputs_loaded",
        trace_sha256 = %inputs.trace.source_sha256,
        call_loci = inputs.trace.call_count(),
        config_sha256 = %inputs.config.source_sha256,
        output_path = ?output.display().to_string(),
    );
    drop(stage);

    let ProcessedRead {
        calls,
        signal,
        callability,
        quality,
        second,
        mixtures,
        ..
    } = read_processing::process(&inputs.trace, &inputs.config.sanger)?;
    let stage = tracing::info_span!("read_evidence").entered();
    let evidence = read_processing::read_evidence(
        &calls,
        &signal,
        &callability,
        &quality,
        &second,
        &mixtures,
        &inputs.config.sanger.sanger_evidence,
    )?;
    drop(stage);

    let _stage = tracing::info_span!("reporting").entered();
    let stage_started = Instant::now();
    let result = report::build_evidence(CompletedEvidence {
        config: inputs.config,
        trace: inputs.trace,
        calls,
        callability,
        quality,
        second,
        evidence,
        plugins: plugin::EVIDENCE,
    })?;
    let bytes = report::serialize(&result)?;
    tracing::info!(
        event = "evidence_ready_for_publication",
        elapsed_ms = stage_started.elapsed().as_millis(),
        total_elapsed_ms = started.elapsed().as_millis(),
        schema = result.schema_version,
        calls = result.read.call_count,
        output_path = ?output.display().to_string(),
        bytes = bytes.len(),
    );
    log.sync()?;
    report::publish(&output, &bytes)
}
