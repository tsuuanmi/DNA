//! Post-calling notation of a variants document (ADR-0069 phase 4).

use std::path::Path;
use std::time::Instant;

use crate::error::{Result, VariantsError};
use crate::input::variants;
use crate::operation_log::OperationLog;
use crate::pipeline::path;
use crate::pipeline::plugins as plugin;
use crate::report::{self, CompletedNotation, NamedRead};
use dna_empop::conformance;
use dna_empop::export::{self, Record};

use super::Operation;
use super::represent;
use super::sample::validate_sample_id;

/// Derives the notation of one sample's variants document with one
/// sample-level append-only log.
pub(crate) fn run(
    sample_id: &str,
    document: &Path,
    reference: &Path,
    config_path: &Path,
) -> Result<()> {
    validate_sample_id(sample_id)?;
    Operation::begin(sample_id, "notation_failed", || {
        tracing::info!(
            event = "notation_started",
            version = env!("CARGO_PKG_VERSION"),
            sample_id = ?sample_id,
            document_path = ?document.display().to_string(),
            reference_path = ?reference.display().to_string(),
        );
    })?
    .run(|log, started| notation(sample_id, document, reference, config_path, log, started))
}

fn notation(
    sample_id: &str,
    document: &Path,
    reference: &Path,
    config_path: &Path,
    log: &OperationLog,
    started: Instant,
) -> Result<()> {
    let stage = tracing::info_span!("input_loading").entered();
    let inputs = variants::load_notation(sample_id, document, reference, config_path)?;
    let output = path::notation_output(sample_id)?;
    tracing::info!(
        event = "notation_inputs_loaded",
        reads = inputs.document.reads.len(),
        document_sha256 = %inputs.document.sha256,
        reference_sha256 = %inputs.reference.sequence_sha256,
        config_sha256 = %inputs.config.source_sha256,
        profile_id = ?inputs.profile.identity().id,
        profile_sha256 = %inputs.profile.identity().sha256,
        output_path = ?output.display().to_string(),
    );

    drop(stage);
    let stage = tracing::info_span!("nomenclature").entered();
    let reads = inputs
        .document
        .reads
        .iter()
        .map(|read| NamedRead {
            sha256: read.sha256.clone(),
            name: read.name.clone(),
        })
        .collect::<Vec<_>>();
    let records = inputs
        .document
        .reads
        .iter()
        .map(|read| Record {
            reference: read.reference.clone(),
            callable: read.callable.clone(),
            excluded: read.excluded.clone(),
        })
        .collect::<Vec<_>>();
    let notation = represent::represent_variants(
        inputs
            .document
            .reads
            .into_iter()
            .map(|read| (read.sha256, read.variants, read.callable)),
        &inputs.reference,
        &inputs.profile,
    )?
    .ok_or(VariantsError::NoNotation)?;

    drop(stage);
    let stage = tracing::info_span!("conformance").entered();
    let conformance = inputs.profile.conformance.as_ref().map(|rules| {
        let findings = notation
            .reads
            .iter()
            .map(|read| conformance::check(&inputs.reference, rules, &read.variants))
            .collect::<Vec<_>>();
        tracing::info!(
            event = "conformance_completed",
            rules = rules.rules.len(),
            findings = findings.iter().map(Vec::len).sum::<usize>(),
        );
        (rules.rules.clone(), findings)
    });

    drop(stage);
    let stage = tracing::info_span!("empop").entered();
    let empop = if inputs.profile.empop {
        let calls = report::notation_calls(&inputs.reference, &notation)?;
        let calls = calls
            .keys()
            .map(|call| export::Call {
                position: call.position,
                ordinal: call.ordinal,
                text: &call.text,
            })
            .collect::<Vec<_>>();
        let exported = export::export(&inputs.profile.regions, &records, &calls);
        tracing::info!(
            event = "empop_completed",
            ranges = exported.ranges.len(),
            calls = exported.mitotype.len(),
            cut = exported.cut.len(),
        );
        Some(exported)
    } else {
        None
    };

    drop(stage);
    let _stage = tracing::info_span!("reporting").entered();
    let plugins = plugin::notation(conformance.is_some(), empop.is_some());
    let result = report::build_notation(CompletedNotation {
        sample_id: sample_id.to_owned(),
        configuration_sha256: inputs.config.source_sha256.clone(),
        reference: inputs.reference,
        profile: inputs.profile.identity().clone(),
        source_sha256: inputs.document.sha256,
        reads,
        notation,
        conformance,
        empop,
        plugins,
    })?;
    let schema_version = result.schema_version;
    let calls = result.notation.calls.len();
    let bytes = report::serialize(&result)?;
    tracing::info!(
        event = "notation_ready_for_publication",
        total_elapsed_ms = started.elapsed().as_millis(),
        schema = schema_version,
        calls,
        output_path = ?output.display().to_string(),
        bytes = bytes.len(),
    );
    log.sync()?;
    report::publish(&output, &bytes)
}
