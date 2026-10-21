//! Single-read analysis orchestration, reporting, and publication: from a
//! trace, or from the read's recorded evidence (PROP-0011).

use std::path::{Path, PathBuf};
use std::time::Instant;

use crate::config::Config;
use crate::error::Result;
use crate::input::{evidence, sanger};
use crate::model::read_observation::ReadObservation;
use crate::model::sanger_read::SangerRead;
use crate::operation_log::OperationLog;
use crate::pipeline::path;
use crate::pipeline::plugins as plugin;
use crate::report::{self, CompletedAnalysis};
use crate::variant_analysis;
use dna_call::model::called_read::CalledRead;
use dna_call::read_call::ReadIdentity;
use dna_kernel::model::reference::Reference;
use dna_kernel::profile::Profile;

use super::Operation;

/// Runs one complete AB1-to-JSON analysis with one per-trace append-only log.
pub(crate) fn run(trace: &Path, reference: &Path, config_path: &Path) -> Result<()> {
    let trace_stem = path::trace_stem(trace)?;
    Operation::begin(trace_stem, "analysis_failed", || {
        tracing::info!(
            event = "analysis_started",
            version = env!("CARGO_PKG_VERSION"),
            trace_path = ?trace.display().to_string(),
            reference_path = ?reference.display().to_string(),
        );
    })?
    .run(|log, started| analyze(trace, trace_stem, reference, config_path, log, started))
}

/// Runs one analysis of a read from its `dna.read_evidence/v1` document,
/// logged and published under the stem of the trace it was written from.
pub(crate) fn run_evidence(document: &Path, reference: &Path, config_path: &Path) -> Result<()> {
    let stem = path::evidence_stem(document)?;
    Operation::begin(stem, "analysis_failed", || {
        tracing::info!(
            event = "analysis_started",
            version = env!("CARGO_PKG_VERSION"),
            evidence_path = ?document.display().to_string(),
            reference_path = ?reference.display().to_string(),
        );
    })?
    .run(|log, started| analyze_evidence(document, stem, reference, config_path, log, started))
}

fn analyze(
    trace: &Path,
    stem: &str,
    reference: &Path,
    config_path: &Path,
    log: &OperationLog,
    analysis_started: Instant,
) -> Result<()> {
    let stage = tracing::info_span!("input_loading").entered();
    let stage_started = Instant::now();
    let prepared = sanger::prepare_analysis(trace, reference, config_path)?;
    let output = path::analysis_output(stem)?;
    let inputs = sanger::load_analysis(prepared)?;
    tracing::info!(
        event = "inputs_loaded",
        elapsed_ms = stage_started.elapsed().as_millis(),
        trace_name = ?inputs.trace.source_name,
        trace_sha256 = %inputs.trace.source_sha256,
        samples = inputs.trace.sample_count(),
        call_loci = inputs.trace.call_count(),
        vendor_primary = inputs.trace.vendor.primary.is_some(),
        vendor_quality = inputs.trace.vendor.qualities.is_some(),
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
    let completed = variant_analysis::observation::build(
        &inputs.trace,
        &inputs.reference,
        &inputs.config,
        &inputs.profile,
    )?;
    let ReadObservation { called, sanger } = completed.read;
    let sanger = SangerRead::of(&sanger)?;
    publish(
        Analyzed {
            config: &inputs.config,
            profile: &inputs.profile,
            reference: inputs.reference,
            called,
            sanger,
            warning_total: completed.warning_total,
            output,
        },
        log,
        analysis_started,
    )
}

fn analyze_evidence(
    document: &Path,
    stem: &str,
    reference: &Path,
    config_path: &Path,
    log: &OperationLog,
    analysis_started: Instant,
) -> Result<()> {
    let stage = tracing::info_span!("input_loading").entered();
    let stage_started = Instant::now();
    let output = path::analysis_output(stem)?;
    let inputs = evidence::load_analysis(document, reference, config_path, plugin::EVIDENCE)?;
    let read = inputs.read;
    tracing::info!(
        event = "inputs_loaded",
        elapsed_ms = stage_started.elapsed().as_millis(),
        trace_name = ?read.input_name,
        trace_sha256 = %read.input_sha256,
        calls = read.evidence.calls().len(),
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
    let stage = tracing::info_span!("read_evidence").entered();
    variant_analysis::observation::record_evidence(&read.evidence, &inputs.config);
    drop(stage);
    let called = variant_analysis::observation::call(
        ReadIdentity {
            input_name: read.input_name,
            input_sha256: read.input_sha256,
            declared: None,
        },
        read.evidence,
        None,
        &inputs.reference,
        &inputs.config,
        &inputs.profile,
    )?;
    let warning_total =
        variant_analysis::observation::warning_total(&read.sanger.read_warnings(), &called);
    publish(
        Analyzed {
            config: &inputs.config,
            profile: &inputs.profile,
            reference: inputs.reference,
            called,
            sanger: read.sanger,
            warning_total,
            output,
        },
        log,
        analysis_started,
    )
}

/// One read placed and called, ready to report.
struct Analyzed<'a> {
    config: &'a Config,
    profile: &'a Profile,
    reference: Reference,
    called: CalledRead,
    sanger: SangerRead,
    warning_total: usize,
    output: PathBuf,
}

/// Builds the analysis document of one read and publishes it.
fn publish(analyzed: Analyzed<'_>, log: &OperationLog, analysis_started: Instant) -> Result<()> {
    let Analyzed {
        config,
        profile,
        reference,
        called,
        sanger,
        warning_total,
        output,
    } = analyzed;
    let stage = tracing::info_span!("reporting").entered();
    let stage_started = Instant::now();
    let trust = dna_call::read_trust(&called.evidence, &config.core.variant_calling);
    let result = report::build_analysis(CompletedAnalysis {
        reference,
        profile: profile.identity().clone(),
        called,
        sanger,
        trust,
        plugins: plugin::ANALYZE,
    })?;
    let result_variants = result.variants.len();
    let schema_version = result.schema_version;
    let bytes = report::serialize(&result)?;

    drop(stage);
    let _stage = tracing::info_span!("result_publication").entered();
    tracing::info!(
        event = "result_ready_for_publication",
        elapsed_ms = stage_started.elapsed().as_millis(),
        total_elapsed_ms = analysis_started.elapsed().as_millis(),
        schema = schema_version,
        variants = result_variants,
        warnings = warning_total,
        output_path = ?output.display().to_string(),
        bytes = bytes.len(),
    );
    log.sync()?;
    report::publish(&output, &bytes)
}
