//! Core-only call orchestration over reviewed consensus sequences (ADR-0069).

use std::path::{Path, PathBuf};
use std::time::Instant;

use crate::error::Result;
use crate::input::sequence;
use crate::operation_log::OperationLog;
use crate::pipeline::path;
use crate::plugin;
use crate::report::{self, CompletedCall};
use crate::variant_analysis::read_call::{self, ReadIdentity};

use super::sample::validate_sample_id;
use super::{Operation, represent};

/// Runs one core-only call with one sample-level append-only log.
pub(crate) fn run(
    sample_id: &str,
    sequences: &[PathBuf],
    reference: &Path,
    config_path: &Path,
) -> Result<()> {
    validate_sample_id(sample_id)?;
    Operation::begin(sample_id, "call_failed", || {
        tracing::info!(
            event = "call_started",
            version = env!("CARGO_PKG_VERSION"),
            sample_id = ?sample_id,
            sequence_files = sequences.len(),
            reference_path = ?reference.display().to_string(),
        );
    })?
    .run(|log, started| call(sample_id, sequences, reference, config_path, log, started))
}

fn call(
    sample_id: &str,
    sequences: &[PathBuf],
    reference: &Path,
    config_path: &Path,
    log: &OperationLog,
    started: Instant,
) -> Result<()> {
    let stage = tracing::info_span!("input_loading").entered();
    let inputs = sequence::load_call(sequences, reference, config_path)?;
    let output = path::call_output(sample_id)?;
    tracing::info!(
        event = "call_inputs_loaded",
        reads = inputs.reads.len(),
        reference_name = ?inputs.reference.name,
        reference_sha256 = %inputs.reference.sequence_sha256,
        config_sha256 = %inputs.config.source_sha256,
        profile_id = ?inputs.profile.identity().id,
        profile_sha256 = %inputs.profile.identity().sha256,
        output_path = ?output.display().to_string(),
    );

    drop(stage);
    let reads = inputs
        .reads
        .iter()
        .map(|read| {
            tracing::info!(
                event = "call_read_started",
                read = ?read.name,
                sha256 = %read.sha256,
                bases = read.sequence.len(),
            );
            read_call::call_read(
                ReadIdentity {
                    input_name: read.name.clone(),
                    input_sha256: read.sha256.clone(),
                },
                sequence::read_evidence(read)?,
                &inputs.reference,
                &inputs.config,
                &inputs.profile,
            )
        })
        .collect::<Result<Vec<_>>>()?;

    let stage = tracing::info_span!("nomenclature").entered();
    let notation = represent::represent(&reads, &inputs.reference, &inputs.profile)?;

    drop(stage);
    let _stage = tracing::info_span!("reporting").entered();
    let plugins = if notation.is_some() {
        plugin::CALL_WITH_NOTATION
    } else {
        plugin::CALL
    };
    let result = report::build_call(CompletedCall {
        sample_id: sample_id.to_owned(),
        reference: inputs.reference,
        profile: inputs.profile.identity().clone(),
        reads,
        notation,
        plugins,
    })?;
    let schema_version = result.schema_version;
    let variants = result
        .reads
        .iter()
        .map(|read| {
            read.variants
                .iter()
                .filter(|variant| variant.eligible)
                .count()
        })
        .sum::<usize>();
    let bytes = report::serialize(&result)?;
    tracing::info!(
        event = "call_ready_for_publication",
        total_elapsed_ms = started.elapsed().as_millis(),
        schema = schema_version,
        reads = result.reads.len(),
        eligible_variants = variants,
        output_path = ?output.display().to_string(),
        bytes = bytes.len(),
    );
    log.sync()?;
    report::publish(&output, &bytes)
}
