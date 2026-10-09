//! Multi-read sample evidence orchestration and publication.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::time::Instant;

use crate::error::{Result, SampleError};
use crate::input::sanger;
use crate::model::sample_evidence::RejectedSampleRead;
use crate::operation_log::OperationLog;
use crate::pipeline::path;
use crate::plugin;
use crate::report::{self, CompletedSampleEvidence, SangerSampleEvidence};
use crate::sample as sample_science;

use super::{Operation, represent, sample_metrics, sample_reads};

/// Runs one sample-evidence operation with one sample-level append-only log.
pub(crate) fn run(
    sample_id: &str,
    traces: &[PathBuf],
    reference: &Path,
    config_path: &Path,
) -> Result<()> {
    validate_sample_id(sample_id)?;
    Operation::begin(sample_id, "sample_failed", || {
        tracing::info!(
            event = "sample_started",
            version = env!("CARGO_PKG_VERSION"),
            sample_id = ?sample_id,
            traces = traces.len(),
            reference_path = ?reference.display().to_string(),
        );
    })?
    .run(|log, started| sample(sample_id, traces, reference, config_path, log, started))
}

fn sample(
    sample_id: &str,
    traces: &[PathBuf],
    reference: &Path,
    config_path: &Path,
    log: &OperationLog,
    started: Instant,
) -> Result<()> {
    let stage = tracing::info_span!("input_loading").entered();
    let stage_started = Instant::now();
    let inputs = sanger::load_sample(traces, reference, config_path)?;
    let output = path::sample_output(sample_id)?;
    tracing::info!(
        event = "sample_inputs_loaded",
        elapsed_ms = stage_started.elapsed().as_millis(),
        sample_id = ?sample_id,
        traces = inputs.traces.len(),
        reference_name = ?inputs.reference.name,
        reference_sha256 = %inputs.reference.sequence_sha256,
        topology = ?inputs.reference.topology,
        reference_bases = inputs.reference.len(),
        config_path = ?inputs.config.source_path.display().to_string(),
        config_sha256 = %inputs.config.source_sha256,
        profile_id = ?inputs.profile.identity().id,
        profile_sha256 = %inputs.profile.identity().sha256,
        output_path = ?output.display().to_string(),
    );

    drop(stage);
    let completed_reads = sample_reads::build(
        &inputs.traces,
        &inputs.reference,
        &inputs.config,
        &inputs.profile,
    )?;
    let warning_total = completed_reads.warning_total;
    let mut reads = Vec::with_capacity(completed_reads.reads.len());
    let mut sanger_reads = BTreeMap::new();
    for read in completed_reads.reads {
        sanger_reads.insert(read.called.input_sha256.clone(), read.sanger);
        reads.push(read.called);
    }
    let rejected = completed_reads
        .rejected
        .iter()
        .map(|read| RejectedSampleRead {
            input_name: read.input_name.clone(),
            input_sha256: read.input_sha256.clone(),
        })
        .collect::<Vec<_>>();
    let sanger_rejected = completed_reads
        .rejected
        .into_iter()
        .map(|read| (read.input_sha256.clone(), read))
        .collect();

    let stage = tracing::info_span!("sample_aggregation").entered();
    let stage_started = Instant::now();
    let evidence = sample_science::aggregate(
        &reads.iter().collect::<Vec<_>>(),
        &rejected,
        &inputs.config.sample_reconciliation,
    )?;
    let metrics = sample_metrics::summarize(&evidence, &sanger_reads);
    tracing::info!(
        event = "sample_aggregation_completed",
        elapsed_ms = stage_started.elapsed().as_millis(),
        reads = evidence.reads.len(),
        rejected_reads = evidence.rejected_reads.len(),
        coverage_segments = evidence.coverage.len(),
        overlaps = evidence.overlaps.len(),
        eligible_overlaps = evidence
            .overlaps
            .iter()
            .filter(|overlap| overlap.eligible)
            .count(),
        locus_differences = evidence.locus_differences.len(),
        "{metrics}"
    );

    drop(stage);
    let stage = tracing::info_span!("nomenclature").entered();
    let notation = represent::represent(&reads, &inputs.reference, &inputs.profile)?;

    drop(stage);
    let stage = tracing::info_span!("reporting").entered();
    let stage_started = Instant::now();
    let plugins = if notation.is_some() {
        plugin::SAMPLE_WITH_NOTATION
    } else {
        plugin::SAMPLE
    };
    let result = report::build_sample(CompletedSampleEvidence {
        sample_id: sample_id.to_owned(),
        reference: inputs.reference,
        profile: inputs.profile.identity().clone(),
        evidence,
        notation,
        plugins,
        sanger: SangerSampleEvidence {
            reads: sanger_reads,
            rejected: sanger_rejected,
        },
    })?;
    let reads = result.reads.len();
    let rejected_reads = result.rejected_reads.len();
    let coverage_segments = result.coverage.len();
    let overlaps = result.overlaps.len();
    let locus_differences = result.locus_differences.len();
    let variants = result.variants.len();
    let schema_version = result.schema_version;
    let bytes = report::serialize(&result)?;

    drop(stage);
    let _stage = tracing::info_span!("result_publication").entered();
    tracing::info!(
        event = "sample_result_ready_for_publication",
        elapsed_ms = stage_started.elapsed().as_millis(),
        total_elapsed_ms = started.elapsed().as_millis(),
        schema = schema_version,
        reads,
        rejected_reads,
        coverage_segments,
        overlaps,
        locus_differences,
        variants,
        read_warnings = warning_total,
        output_path = ?output.display().to_string(),
        bytes = bytes.len(),
    );
    log.sync()?;
    report::publish(&output, &bytes)
}

pub(super) fn validate_sample_id(sample_id: &str) -> Result<()> {
    let mut characters = sample_id.chars();
    let valid_first = characters
        .next()
        .is_some_and(|value| value.is_ascii_alphanumeric());
    let valid_rest =
        characters.all(|value| value.is_ascii_alphanumeric() || matches!(value, '_' | '.' | '-'));
    if sample_id.len() > 128 || !valid_first || !valid_rest {
        return Err(SampleError::InvalidSampleId.into());
    }
    Ok(())
}
