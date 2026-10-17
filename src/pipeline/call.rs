//! Core-only call orchestration over reviewed consensus sequences (ADR-0069).

use std::path::{Path, PathBuf};
use std::time::Instant;

use crate::config::Config;
use crate::input::sequence::{self, CallInputs, SequenceRead};
use crate::operation_log::OperationLog;
use crate::pipeline::path;
use crate::pipeline::plugins as plugin;
use crate::report::{self, CompletedCall, SampleNotation};
use dna_core::model::called_read::CalledRead;
use dna_core::read_call::{self, CoreRun, ReadIdentity};
use dna_kernel::error::Result;
use dna_kernel::model::reference::Reference;
use dna_kernel::profile::Profile;

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
    let completed = complete(sample_id, inputs)?;

    let _stage = tracing::info_span!("reporting").entered();
    let result = report::build_call(completed)?;
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

/// Calls every consensus sequence and derives the notation over loaded inputs,
/// without output paths, operation logs, or publication.
pub(super) fn complete(sample_id: &str, inputs: CallInputs) -> Result<CompletedCall> {
    let Called { reads, notation } = call_sequences(
        &inputs.reads,
        &inputs.config,
        &inputs.profile,
        &inputs.reference,
    )?;
    let plugins = if notation.is_some() {
        plugin::CALL_WITH_NOTATION
    } else {
        plugin::CALL
    };
    Ok(CompletedCall {
        sample_id: sample_id.to_owned(),
        reference: inputs.reference,
        profile: inputs.profile.identity().clone(),
        reads,
        notation,
        plugins,
    })
}

/// The called reads of a set of consensus sequences and their notation.
pub(super) struct Called {
    pub(super) reads: Vec<CalledRead>,
    /// `None` when the profile declares no notation.
    pub(super) notation: Option<SampleNotation>,
}

/// The reference interval a record name declares, as a consensus segment's
/// `<sample-id>_<first>-<last>` does (1-based, closed, within the reference);
/// `None` for any other name.
fn declared(name: &str, reference_length: usize) -> Option<[usize; 2]> {
    let (_, interval) = name.rsplit_once('_')?;
    let (first, last) = interval.split_once('-')?;
    let digits = |text: &str| !text.is_empty() && text.bytes().all(|byte| byte.is_ascii_digit());
    if !digits(first) || !digits(last) {
        return None;
    }
    let (first, last) = (first.parse::<usize>().ok()?, last.parse::<usize>().ok()?);
    (1 <= first && first <= last && last <= reference_length).then_some([first, last])
}

/// The scientific step of `call`: calls every consensus sequence with the
/// core, then derives the notation.
pub(super) fn call_sequences(
    sequences: &[SequenceRead],
    config: &Config,
    profile: &Profile,
    reference: &Reference,
) -> Result<Called> {
    let run = CoreRun {
        config: &config.core,
        configuration_sha256: &config.source_sha256,
        regions: &profile.regions,
    };
    let reads = sequences
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
                    declared: declared(&read.name, reference.sequence.len()),
                },
                sequence::read_evidence(read)?,
                reference,
                &run,
            )
            .map(CalledRead::without_excluded_positions)
        })
        .collect::<Result<Vec<_>>>()?;

    let _stage = tracing::info_span!("nomenclature").entered();
    let notation = represent::represent(&reads, reference, profile)?;
    Ok(Called { reads, notation })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_consensus_record_name_declares_its_interval() {
        assert_eq!(
            declared("2514526_15803-15920", 16_569),
            Some([15_803, 15_920])
        );
        assert_eq!(declared("a_b_1-70", 160), Some([1, 70]));
        for name in [
            "HV1",
            "sample_70-1",
            "sample_0-5",
            "sample_1-161",
            "sample_1-",
            "sample_+1-5",
            "sample-1-5",
        ] {
            assert_eq!(declared(name, 160), None, "{name}");
        }
    }
}
