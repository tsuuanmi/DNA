//! Sample review orchestration and publication (PROP-0004).

use std::path::{Path, PathBuf};
use std::time::Instant;

use crate::error::Result;
use crate::input::sanger;
use crate::input::signal as input_signal;
use crate::model::review_result::{ItemKind, LabelResult};
use crate::operation_log::OperationLog;
use crate::pipeline::path;
use crate::pipeline::plugins as plugin;
use crate::report::{self, CompletedReview, Reviewed};
use crate::sample_analysis::{self, ConsensusCall, ReviewCheck, SampleAnalysis, SampleOutcome};
use dna_sanger::run_signal::RunSignal;

use super::Operation;
use super::sample::validate_sample_id;

/// Runs one sample review operation with one sample-level append-only log.
pub(crate) fn run(
    sample_id: &str,
    traces: &[PathBuf],
    reference: &Path,
    signal: Option<&Path>,
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
    .run(|log, started| {
        review(
            sample_id,
            traces,
            reference,
            signal,
            config_path,
            log,
            started,
        )
    })
}

fn review(
    sample_id: &str,
    traces: &[PathBuf],
    reference: &Path,
    signal: Option<&Path>,
    config_path: &Path,
    log: &OperationLog,
    started: Instant,
) -> Result<()> {
    let stage = tracing::info_span!("input_loading").entered();
    let inputs = sanger::load_sample(traces, reference, config_path)?;
    let signal = signal.map(input_signal::load).transpose()?;
    if let Some((signal, sha256)) = &signal {
        tracing::info!(
            event = "review_signal_loaded",
            run_id = ?signal.run_id,
            sha256 = %sha256,
            reads = signal.reads,
            keys = signal.keys.len(),
        );
    }
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
            let call = analysis.call(sample_id);
            if let Err(error) = &call {
                tracing::warn!(event = "review_call_failed", error = ?error.to_string());
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
        signal,
    })?;
    let judged = result
        .items
        .iter()
        .flat_map(|item| &item.labels)
        .filter_map(|label| match label {
            LabelResult::PossibleHeteroplasmy { reads, .. } => Some(reads),
            _ => None,
        })
        .flatten()
        .map(|read| read.calls.is_some())
        .collect::<Vec<_>>();
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
        heteroplasmy_reads_run_level = judged.iter().filter(|&&level| level).count(),
        heteroplasmy_reads_floor = judged.iter().filter(|&&level| !level).count(),
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

impl SampleAnalysis {
    /// The items of this analysis's review list that a reviewer must look at
    /// (`check` items), given the call of its consensus and the run signal
    /// it is judged against, if any, by the rules `dna review` applies,
    /// without building or publishing the document.
    ///
    /// # Errors
    ///
    /// Returns [`ReportError`](crate::error::ReportError) when the analysis's
    /// records disagree, as `dna review` would.
    pub fn review_checks(
        &self,
        call: &ConsensusCall,
        signal: Option<&RunSignal>,
    ) -> Result<Vec<ReviewCheck>> {
        Ok(report::review_checks(self, call, signal)?
            .into_iter()
            .map(|(first, last, labels)| ReviewCheck {
                first,
                last,
                labels,
            })
            .collect())
    }
}
