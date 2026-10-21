//! Run context construction and publication over one laboratory batch
//! (PROP-0014).

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

use crate::error::{KernelError, Result, RunContextError};
use crate::operation_log::OperationLog;
use crate::pipeline::path;
use crate::pipeline::plugins as plugin;
use crate::report::{self, CompletedRunContext, Skipped};
use crate::sample_analysis::{self, Outcome, RunContext};
use dna_consensus::noise::LengthStatistics;

use super::sample::validate_sample_id;
use super::{Operation, SampleInput};

/// Largest accepted read list.
const MAX_LIST_BYTES: u64 = 16 * 1024 * 1024;

/// A batch's read list: `sample_id<TAB>path` lines of traces or of read
/// evidence documents.
#[derive(Debug, Clone, Copy)]
pub(crate) enum List<'a> {
    Traces(&'a Path),
    Evidence(&'a Path),
}

/// Runs one run context operation with one run-level append-only log.
pub(crate) fn run(
    run_id: &str,
    list: List<'_>,
    reference: &Path,
    config_path: &Path,
) -> Result<()> {
    validate_sample_id(run_id)?;
    Operation::begin(run_id, "run_context_failed", || {
        tracing::info!(
            event = "run_context_started",
            version = env!("CARGO_PKG_VERSION"),
            run_id = ?run_id,
            evidence = matches!(list, List::Evidence(_)),
            reference_path = ?reference.display().to_string(),
        );
    })?
    .run(|log, started| context(run_id, list, reference, config_path, log, started))
}

/// The samples of a read list, in first-seen order, each with its reads.
fn samples(list: &Path) -> Result<Vec<(String, Vec<PathBuf>)>> {
    let read_error = |source| KernelError::Read {
        kind: "run read list",
        path: list.to_path_buf(),
        source,
    };
    let bytes = fs::metadata(list).map_err(read_error)?.len();
    if bytes == 0 || bytes > MAX_LIST_BYTES {
        return Err(RunContextError::ReadList {
            line: 0,
            reason: "the list is empty or too large",
        }
        .into());
    }
    let text = fs::read_to_string(list).map_err(read_error)?;
    let mut samples: Vec<(String, Vec<PathBuf>)> = Vec::new();
    for (index, line) in text.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let Some((sample, path)) = line.split_once('\t') else {
            return Err(RunContextError::ReadList {
                line: index + 1,
                reason: "a line is not sample_id<TAB>path",
            }
            .into());
        };
        validate_sample_id(sample)?;
        match samples.iter_mut().find(|(known, _)| known == sample) {
            Some((_, paths)) => paths.push(PathBuf::from(path)),
            None => samples.push((sample.to_owned(), vec![PathBuf::from(path)])),
        }
    }
    if samples.is_empty() {
        return Err(RunContextError::ReadList {
            line: 0,
            reason: "the list names no read",
        }
        .into());
    }
    Ok(samples)
}

/// Why a listed input file is left out of the run context, if it is: a
/// file that is missing, not a regular file, empty, or unreadable. Such a
/// file fails every sample command, so its sample cannot be analysed with
/// it; its other reads still count, and the sample's own analysis, which
/// cannot use the file either, still leaves them out.
fn unusable(path: &Path) -> Option<&'static str> {
    match fs::metadata(path) {
        Err(_) => Some("missing"),
        Ok(metadata) if !metadata.is_file() => Some("not_a_file"),
        Ok(metadata) if metadata.len() == 0 => Some("empty"),
        Ok(_) => fs::File::open(path).err().map(|_| "unreadable"),
    }
}

/// A listed file's name, as a run context records it.
fn file_name(path: &Path) -> String {
    path.file_name().map_or_else(
        || path.display().to_string(),
        |name| name.to_string_lossy().into_owned(),
    )
}

fn context(
    run_id: &str,
    list: List<'_>,
    reference: &Path,
    config_path: &Path,
    log: &OperationLog,
    started: Instant,
) -> Result<()> {
    let samples = samples(match list {
        List::Traces(path) | List::Evidence(path) => path,
    })?;
    let output = path::run_context_output(run_id)?;
    let stage_started = Instant::now();
    let mut statistics: Option<LengthStatistics> = None;
    let mut reads = BTreeSet::new();
    let mut rejected = 0;
    let mut contributing = 0;
    let mut identity = None;
    let mut windows = Vec::new();
    let mut skipped = Vec::new();
    for (sample_id, listed) in &samples {
        let mut paths = Vec::new();
        for path in listed {
            match unusable(path) {
                Some(reason) => {
                    tracing::warn!(
                        event = "run_context_input_skipped",
                        sample_id = ?sample_id,
                        path = ?path.display().to_string(),
                        reason,
                    );
                    skipped.push(Skipped {
                        sample_id: sample_id.clone(),
                        name: file_name(path),
                        reason,
                    });
                }
                None => paths.push(path.clone()),
            }
        }
        if paths.is_empty() {
            continue;
        }
        let input = match list {
            List::Traces(_) => SampleInput::Traces(&paths),
            List::Evidence(_) => SampleInput::Evidence(&paths),
        };
        let super::LoadedSample {
            loaded,
            reads: sample_reads,
        } = super::load_sample(input, reference, config_path)?;
        if identity.is_none() {
            windows = sample_analysis::recorded_windows(&loaded.profile);
            identity = Some((
                loaded.reference.clone(),
                loaded.profile.identity().clone(),
                loaded.config.source_sha256.clone(),
            ));
        }
        match sample_analysis::analyze(sample_reads.reads(), loaded, log, super::read_workers())? {
            Outcome::Rejected { rejections, .. } => rejected += rejections.len(),
            Outcome::Analyzed(analysis) => {
                let pairs = analysis.length_pairs();
                tracing::info!(
                    event = "run_context_sample",
                    sample_id = ?sample_id,
                    reads = analysis.read_identities().len(),
                    pairs = pairs.len(),
                );
                rejected += analysis.rejected_reads();
                reads.extend(analysis.read_identities());
                contributing += 1;
                statistics = Some(
                    statistics
                        .unwrap_or_else(|| analysis.empty_length_statistics())
                        .with(&pairs),
                );
            }
        }
    }
    let (Some(length), Some((reference, profile, configuration_sha256))) = (statistics, identity)
    else {
        return Err(RunContextError::Inconsistent("no sample of the batch was analysed").into());
    };
    tracing::info!(
        event = "run_context_completed",
        elapsed_ms = stage_started.elapsed().as_millis(),
        samples = contributing,
        reads = reads.len(),
        rejected_reads = rejected,
        run1_pairs = length.runs[0].pairs,
        run2_pairs = length.runs[1].pairs,
    );
    let result = report::build_run_context(CompletedRunContext {
        context: RunContext {
            run_id: run_id.to_owned(),
            reads,
            rejected_reads: rejected,
            reference_sha256: reference.sequence_sha256.clone(),
            windows,
            length,
            sha256: None,
        },
        samples: contributing,
        skipped,
        reference,
        profile,
        configuration_sha256,
        plugins: plugin::CONSENSUS,
    });
    let bytes = report::serialize(&result)?;
    tracing::info!(
        event = "run_context_ready_for_publication",
        total_elapsed_ms = started.elapsed().as_millis(),
        schema = result.schema_version,
        output_path = ?output.display().to_string(),
        bytes = bytes.len(),
    );
    log.sync()?;
    report::publish(&output, &bytes)
}
