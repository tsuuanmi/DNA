//! Sample review orchestration and publication (PROP-0004).

use std::path::{Path, PathBuf};
use std::time::Instant;

use crate::input::sanger;
use crate::model::review_result::ItemKind;
use crate::operation_log::OperationLog;
use crate::pipeline::path;
use crate::pipeline::plugins as plugin;
use crate::report::{self, CompletedReview, Reviewed};
use crate::sample_analysis::{self, SampleOutcome};
use dna_kernel::error::Result;

use super::Operation;
use super::sample::validate_sample_id;

/// Runs one sample review operation with one sample-level append-only log.
pub(crate) fn run(
    sample_id: &str,
    traces: &[PathBuf],
    reference: &Path,
    config_path: &Path,
) -> Result<()> {
    validate_sample_id(sample_id)?;
    Operation::begin(sample_id, "review_failed", || {
        tracing::info!(
            event = "review_started",
            version = env!("CARGO_PKG_VERSION"),
            sample_id = ?sample_id,
            traces = traces.len(),
            reference_path = ?reference.display().to_string(),
        );
    })?
    .run(|log, started| review(sample_id, traces, reference, config_path, log, started))
}

fn review(
    sample_id: &str,
    traces: &[PathBuf],
    reference: &Path,
    config_path: &Path,
    log: &OperationLog,
    started: Instant,
) -> Result<()> {
    let stage = tracing::info_span!("input_loading").entered();
    let inputs = sanger::load_sample(traces, reference, config_path)?;
    let output = path::review_output(sample_id)?;
    tracing::info!(
        event = "review_inputs_loaded",
        traces = inputs.traces.len(),
        reference_sha256 = %inputs.reference.sequence_sha256,
        config_sha256 = %inputs.config.source_sha256,
        profile_id = ?inputs.profile.identity().id,
        output_path = ?output.display().to_string(),
    );
    let reference = inputs.reference.clone();
    let profile = inputs.profile.clone();
    let configuration_sha256 = inputs.config.source_sha256.clone();

    drop(stage);
    let outcome = sample_analysis::analyze(inputs, log, super::read_workers())?;

    let stage = tracing::info_span!("review").entered();
    let stage_started = Instant::now();
    let outcome = match outcome {
        SampleOutcome::Rejected(rejected) => Reviewed::Rejected(rejected),
        SampleOutcome::Analyzed(analysis) => {
            let call = analysis.call(sample_id).map_err(|error| error.to_string());
            if let Err(error) = &call {
                tracing::warn!(event = "review_call_failed", error = ?error);
            }
            Reviewed::Analyzed { analysis, call }
        }
    };
    let plugins = plugin::review(profile.notation.is_some(), profile.empop);
    let result = report::build_review(CompletedReview {
        sample_id: sample_id.to_owned(),
        reference,
        profile,
        configuration_sha256,
        outcome,
        plugins,
    })?;
    tracing::info!(
        event = "review_completed",
        elapsed_ms = stage_started.elapsed().as_millis(),
        status = ?result.status,
        reads = result.reads.len(),
        items = result.items.len(),
        checks = result
            .items
            .iter()
            .filter(|item| item.kind == ItemKind::Check)
            .count(),
        clear = result.summary.clear,
    );

    drop(stage);
    let _stage = tracing::info_span!("reporting").entered();
    let bytes = report::serialize(&result)?;
    tracing::info!(
        event = "review_ready_for_publication",
        total_elapsed_ms = started.elapsed().as_millis(),
        schema = result.schema_version,
        output_path = ?output.display().to_string(),
        bytes = bytes.len(),
    );
    log.sync()?;
    report::publish(&output, &bytes)
}
