//! Sample consensus orchestration and publication (PROP-0003).

use std::path::Path;
use std::time::Instant;

use crate::error::Result;
use crate::operation_log::OperationLog;
use crate::pipeline::path;
use crate::pipeline::plugins as plugin;
use crate::report::{self, CompletedConsensus};
use crate::sample_analysis::{self, reads};

use super::sample::validate_sample_id;
use super::{Operation, SampleInput};

/// Runs one sample consensus operation with one sample-level append-only log.
pub(crate) fn run(
    sample_id: &str,
    input: SampleInput<'_>,
    reference: &Path,
    run_context: Option<&Path>,
    config_path: &Path,
) -> Result<()> {
    validate_sample_id(sample_id)?;
    Operation::begin(sample_id, "consensus_failed", || {
        tracing::info!(
            event = "consensus_started",
            version = env!("CARGO_PKG_VERSION"),
            sample_id = ?sample_id,
            traces = input.len(),
            evidence = matches!(input, SampleInput::Evidence(_)),
            reference_path = ?reference.display().to_string(),
        );
    })?
    .run(|log, started| {
        build(
            sample_id,
            input,
            reference,
            run_context,
            config_path,
            log,
            started,
        )
    })
}

fn build(
    sample_id: &str,
    input: SampleInput<'_>,
    reference: &Path,
    run_context: Option<&Path>,
    config_path: &Path,
    log: &OperationLog,
    started: Instant,
) -> Result<()> {
    let stage = tracing::info_span!("input_loading").entered();
    let super::LoadedSample {
        loaded: mut inputs,
        reads,
    } = super::load_sample(input, reference, config_path)?;
    inputs.run_context = super::load_run_context(run_context)?;
    let (document, fasta) = path::consensus_outputs(sample_id)?;
    tracing::info!(
        event = "consensus_inputs_loaded",
        traces = reads.len(),
        reference_sha256 = %inputs.reference.sequence_sha256,
        config_sha256 = %inputs.config.source_sha256,
        profile_id = ?inputs.profile.identity().id,
        output_path = ?document.display().to_string(),
        fasta_path = ?fasta.display().to_string(),
    );

    drop(stage);
    let placed = reads::place(
        reads.reads(),
        &inputs.reference,
        &inputs.config,
        &inputs.profile,
        log,
        super::read_workers(),
    )?;

    let stage = tracing::info_span!("consensus").entered();
    let stage_started = Instant::now();
    let evidence = sample_analysis::aggregate(&placed, &inputs.config)?;
    let (consensus, run_context) =
        sample_analysis::build_consensus(&placed.called, &evidence, &inputs)?;
    tracing::info!(
        event = "consensus_completed",
        elapsed_ms = stage_started.elapsed().as_millis(),
        reads = evidence.reads.len(),
        rejected_reads = evidence.rejected_reads.len(),
        segments = consensus.segments.len(),
        sites = consensus.sites.len(),
        called_positions = consensus.summary.called_positions,
        contested_sites = consensus.summary.contested_sites,
        unresolved_positions = consensus.summary.unresolved_positions,
    );

    drop(stage);
    let _stage = tracing::info_span!("reporting").entered();
    let result = report::build_consensus(CompletedConsensus {
        sample_id: sample_id.to_owned(),
        reference: inputs.reference,
        profile: inputs.profile.identity().clone(),
        evidence,
        rejected: sample_analysis::by_identity(placed.rejected),
        consensus,
        plugins: plugin::CONSENSUS,
        run_context,
    })?;
    let records = report::consensus_fasta(&result);
    let bytes = report::serialize(&result)?;
    tracing::info!(
        event = "consensus_ready_for_publication",
        total_elapsed_ms = started.elapsed().as_millis(),
        schema = result.schema_version,
        segments = result.segments.len(),
        output_path = ?document.display().to_string(),
        bytes = bytes.len(),
    );
    log.sync()?;
    report::publish(&fasta, &records)?;
    report::publish(&document, &bytes)
}
