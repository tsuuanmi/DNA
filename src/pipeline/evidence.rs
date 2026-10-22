//! Reference-free read evidence and result publication (PROP-0011): the
//! modality's per-call evidence exactly as the core receives it, or why the
//! modality set the read aside.

use std::path::Path;
use std::time::Instant;

use crate::error::Result;
use crate::input::sanger::{self, BasecallInputs};
use crate::model::read_observation::{RejectedRead, RejectionCause};
use crate::model::sanger_read::{Rejection, SangerRead};
use crate::operation_log::OperationLog;
use crate::pipeline::path;
use crate::pipeline::plugins as plugin;
use crate::report::{self, CompletedEvidence, Evidenced};
use dna_sanger::read_processing::{self, Admission};

use super::Operation;

/// Runs one complete AB1-to-evidence JSON operation.
pub(crate) fn run(trace: &Path, config_path: &Path) -> Result<()> {
    let trace_stem = path::trace_stem(trace)?;
    Operation::begin(trace_stem, "evidence_failed", || {
        tracing::info!(
            event = "evidence_started",
            version = env!("CARGO_PKG_VERSION"),
            trace_path = ?trace.display().to_string(),
        );
    })?
    .run(|log, started| evidence(trace, config_path, log, started))
}

fn evidence(trace: &Path, config_path: &Path, log: &OperationLog, started: Instant) -> Result<()> {
    let stage = tracing::info_span!("input_loading").entered();
    let prepared = sanger::prepare_basecall(trace, config_path)?;
    let output = path::evidence_output(trace)?;
    let inputs = sanger::load_basecall(prepared)?;
    tracing::info!(
        event = "evidence_inputs_loaded",
        trace_sha256 = %inputs.trace.source_sha256,
        call_loci = inputs.trace.call_count(),
        config_sha256 = %inputs.config.source_sha256,
        output_path = ?output.display().to_string(),
    );
    drop(stage);

    let outcome = outcome(&inputs)?;

    let _stage = tracing::info_span!("reporting").entered();
    let stage_started = Instant::now();
    let calls = match &outcome {
        Evidenced::Read { evidence, .. } => evidence.calls().len(),
        Evidenced::Rejected(_) => 0,
    };
    let result = report::build_evidence(CompletedEvidence {
        input: (inputs.trace.source_name, inputs.trace.source_sha256),
        config: inputs.config,
        outcome,
        plugins: plugin::EVIDENCE,
    })?;
    let bytes = report::serialize(&result)?;
    tracing::info!(
        event = "evidence_ready_for_publication",
        elapsed_ms = stage_started.elapsed().as_millis(),
        total_elapsed_ms = started.elapsed().as_millis(),
        schema = result.schema_version(),
        calls,
        output_path = ?output.display().to_string(),
        bytes = bytes.len(),
    );
    log.sync()?;
    report::publish(&output, &bytes)
}

/// What the Sanger plugin makes of one loaded trace: the read with its
/// evidence, or the read set aside as a sample command sets it aside (too
/// short for one signal window, or with too few callable calls).
fn outcome(inputs: &BasecallInputs) -> Result<Evidenced> {
    let rejected = |cause| {
        let read = RejectedRead {
            input_name: inputs.trace.source_name.clone(),
            input_sha256: inputs.trace.source_sha256.clone(),
            cause,
        };
        tracing::warn!(
            event = "evidence_read_rejected",
            trace_sha256 = %read.input_sha256,
        );
        Evidenced::Rejected(Box::new(Rejection::of(&read)))
    };
    Ok(
        match read_processing::admit(&inputs.trace, &inputs.config.sanger)? {
            Admission::TooShort(short) => rejected(RejectionCause::TooShort(short)),
            Admission::TooFewCallable { sanger, rejection } => {
                rejected(RejectionCause::Callability {
                    sanger: *sanger,
                    rejection,
                })
            }
            Admission::Processed(processed) => {
                let stage = tracing::info_span!("read_evidence").entered();
                let evidence = read_processing::read_evidence(
                    &processed,
                    &inputs.config.sanger.sanger_evidence,
                )?;
                drop(stage);
                let sanger = SangerRead::project(
                    &processed.calls,
                    &processed.signal,
                    &processed.callability,
                    &processed.quality,
                    &processed.second,
                )?;
                Evidenced::Read {
                    sanger: Box::new(sanger),
                    evidence: Box::new(evidence),
                }
            }
        },
    )
}

/// The read evidence document of `trace` under the configuration at
/// `config_path`, as `dna evidence` writes it, in memory (PROP-0011): what
/// an evaluation reads back in place of the trace.
///
/// # Errors
///
/// As `dna evidence`: an unreadable or invalid trace or configuration, or a
/// failure of the Sanger plugin other than setting the read aside.
pub(crate) fn document(trace: &Path, config_path: &Path) -> Result<Vec<u8>> {
    let inputs = sanger::load_basecall(sanger::prepare_basecall(trace, config_path)?)?;
    let outcome = outcome(&inputs)?;
    let result = report::build_evidence(CompletedEvidence {
        input: (inputs.trace.source_name, inputs.trace.source_sha256),
        config: inputs.config,
        outcome,
        plugins: plugin::EVIDENCE,
    })?;
    report::serialize(&result)
}
