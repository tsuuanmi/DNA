//! Run signal construction and publication (PROP-0009).

use std::collections::BTreeSet;
use std::path::PathBuf;
use std::time::Instant;

use crate::error::Result;
use crate::input::sanger;
use crate::operation_log::OperationLog;
use crate::pipeline::path;
use crate::pipeline::plugins as plugin;
use crate::report::{self, CompletedRunSignal};
use dna_sanger::run_signal::{self, RunSignal};

use super::Operation;
use super::sample::validate_sample_id;

/// Runs one run signal operation with one run-level append-only log.
pub(crate) fn run(run_id: &str, traces: &[PathBuf], config_path: &std::path::Path) -> Result<()> {
    validate_sample_id(run_id)?;
    Operation::begin(run_id, "signal_failed", || {
        tracing::info!(
            event = "signal_started",
            version = env!("CARGO_PKG_VERSION"),
            run_id = ?run_id,
            traces = traces.len(),
        );
    })?
    .run(|log, started| signal(run_id, traces, config_path, log, started))
}

fn signal(
    run_id: &str,
    traces: &[PathBuf],
    config_path: &std::path::Path,
    log: &OperationLog,
    started: Instant,
) -> Result<()> {
    let stage = tracing::info_span!("input_loading").entered();
    let inputs = sanger::load_run(traces, config_path)?;
    let output = path::signal_output(run_id)?;
    tracing::info!(
        event = "signal_inputs_loaded",
        traces = inputs.traces.len(),
        config_sha256 = %inputs.config.source_sha256,
        output_path = ?output.display().to_string(),
    );
    drop(stage);

    let stage = tracing::info_span!("run_signal").entered();
    let stage_started = Instant::now();
    let edge_margin = inputs.config.core.variant_calling.read_end_margin;
    let mut contributions = Vec::with_capacity(inputs.traces.len());
    let mut rejected = 0;
    for trace in &inputs.traces {
        if let Some(read) = run_signal::read(trace, &inputs.config.sanger, edge_margin)? {
            contributions.push(read);
        } else {
            rejected += 1;
            tracing::warn!(
                event = "signal_read_rejected",
                trace_sha256 = %trace.source_sha256,
            );
        }
    }
    let identities = inputs
        .traces
        .iter()
        .map(|trace| trace.source_sha256.clone())
        .collect::<BTreeSet<_>>();
    let signal = RunSignal::build(run_id.to_owned(), identities, contributions);
    tracing::info!(
        event = "run_signal_completed",
        elapsed_ms = stage_started.elapsed().as_millis(),
        reads = signal.reads,
        rejected_reads = rejected,
        calls = signal.calls,
        keys = signal.keys.len(),
        selection = signal.selection,
    );
    drop(stage);

    let stage = tracing::info_span!("reporting").entered();
    let result = report::build_run_signal(CompletedRunSignal {
        signal,
        rejected_reads: rejected,
        configuration_sha256: inputs.config.source_sha256.clone(),
        plugins: plugin::SIGNAL,
    });
    let bytes = report::serialize(&result)?;
    drop(stage);
    let _stage = tracing::info_span!("result_publication").entered();
    tracing::info!(
        event = "signal_ready_for_publication",
        total_elapsed_ms = started.elapsed().as_millis(),
        schema = result.schema_version,
        output_path = ?output.display().to_string(),
        bytes = bytes.len(),
    );
    log.sync()?;
    report::publish(&output, &bytes)
}
