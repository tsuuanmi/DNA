//! Sample consensus orchestration and publication (PROP-0003).

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::time::Instant;

use crate::input::sanger;
use crate::operation_log::OperationLog;
use crate::pipeline::path;
use crate::pipeline::plugins as plugin;
use crate::report::{self, CompletedConsensus};
use dna_core::consensus;
use dna_core::model::sample_evidence::RejectedSampleRead;
use dna_core::sample as sample_science;
use dna_kernel::error::Result;

use super::sample::validate_sample_id;
use super::{Operation, sample_reads};

/// Runs one sample consensus operation with one sample-level append-only log.
pub(crate) fn run(
    sample_id: &str,
    traces: &[PathBuf],
    reference: &Path,
    config_path: &Path,
) -> Result<()> {
    validate_sample_id(sample_id)?;
    Operation::begin(sample_id, "consensus_failed", || {
        tracing::info!(
            event = "consensus_started",
            version = env!("CARGO_PKG_VERSION"),
            sample_id = ?sample_id,
            traces = traces.len(),
            reference_path = ?reference.display().to_string(),
        );
    })?
    .run(|log, started| build(sample_id, traces, reference, config_path, log, started))
}

fn build(
    sample_id: &str,
    traces: &[PathBuf],
    reference: &Path,
    config_path: &Path,
    log: &OperationLog,
    started: Instant,
) -> Result<()> {
    let stage = tracing::info_span!("input_loading").entered();
    let inputs = sanger::load_sample(traces, reference, config_path)?;
    let (document, fasta) = path::consensus_outputs(sample_id)?;
    tracing::info!(
        event = "consensus_inputs_loaded",
        traces = inputs.traces.len(),
        reference_sha256 = %inputs.reference.sequence_sha256,
        config_sha256 = %inputs.config.source_sha256,
        profile_id = ?inputs.profile.identity().id,
        output_path = ?document.display().to_string(),
        fasta_path = ?fasta.display().to_string(),
    );

    drop(stage);
    let completed_reads = sample_reads::build(
        &inputs.traces,
        &inputs.reference,
        &inputs.config,
        &inputs.profile,
    )?;
    let reads = completed_reads
        .reads
        .into_iter()
        .map(|read| read.called)
        .collect::<Vec<_>>();
    let rejected = completed_reads
        .rejected
        .iter()
        .map(|read| RejectedSampleRead {
            input_name: read.sanger.input_name.clone(),
            input_sha256: read.sanger.input_sha256.clone(),
        })
        .collect::<Vec<_>>();
    let causes = completed_reads
        .rejected
        .into_iter()
        .map(|read| (read.sanger.input_sha256.clone(), read))
        .collect::<BTreeMap<_, _>>();

    let stage = tracing::info_span!("consensus").entered();
    let stage_started = Instant::now();
    // Aggregation validates the read set and fixes the read order.
    let evidence = sample_science::aggregate(
        &reads.iter().collect::<Vec<_>>(),
        &rejected,
        &inputs.config.core.sample_reconciliation,
    )?;
    let by_identity = reads
        .iter()
        .map(|read| (read.input_sha256.as_str(), read))
        .collect::<BTreeMap<_, _>>();
    let ordered = evidence
        .reads
        .iter()
        .filter_map(|read| by_identity.get(read.input_sha256.as_str()).copied())
        .collect::<Vec<_>>();
    let consensus = consensus::build(&ordered, &inputs.reference, &inputs.config.core);
    tracing::info!(
        event = "consensus_completed",
        elapsed_ms = stage_started.elapsed().as_millis(),
        reads = ordered.len(),
        rejected_reads = rejected.len(),
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
        rejected: causes,
        consensus,
        plugins: plugin::CONSENSUS,
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
