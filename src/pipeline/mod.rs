//! End-to-end orchestration for DNA operations.

mod analyze;
mod basecall;
mod call;
mod consensus;
mod evidence;
mod notation;
mod path;
pub(crate) mod plugins;
mod represent;
mod review;
mod run_context;
mod sample;
mod sample_call;
mod sample_metrics;
mod signal;

use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;
use std::thread;
use std::time::Instant;

use crate::error::{Error, Result};
use crate::input::evidence::{self as input_evidence, Recorded};
use crate::input::sanger;
use crate::operation_log::{OperationLog, Records};
use crate::sample_analysis::Loaded;
use crate::sample_analysis::reads::{ReadScope, Reads};
use dna_sanger::model::sanger::Chromatogram;

/// The reads of a sample command: AB1 traces, or the read evidence documents
/// `evidence` wrote for them (PROP-0011).
#[derive(Debug, Clone, Copy)]
pub(crate) enum SampleInput<'a> {
    /// AB1 traces.
    Traces(&'a [PathBuf]),
    /// `dna.read_evidence/v1` documents.
    Evidence(&'a [PathBuf]),
}

impl SampleInput<'_> {
    /// Number of reads.
    const fn len(self) -> usize {
        match self {
            Self::Traces(paths) | Self::Evidence(paths) => paths.len(),
        }
    }
}

/// One sample's loaded configuration, profile, reference, and reads.
struct LoadedSample {
    loaded: Loaded,
    reads: LoadedReads,
}

/// A sample's reads as loaded: decoded traces, or recorded evidence.
enum LoadedReads {
    Traces(Vec<Chromatogram>),
    Recorded(Vec<Recorded>),
}

impl LoadedReads {
    fn reads(&self) -> Reads<'_> {
        match self {
            Self::Traces(traces) => Reads::Traces(traces),
            Self::Recorded(reads) => Reads::Recorded(reads),
        }
    }

    fn len(&self) -> usize {
        match self {
            Self::Traces(traces) => traces.len(),
            Self::Recorded(reads) => reads.len(),
        }
    }
}

/// Validates and loads one sample's inputs; recorded evidence must come from
/// the loaded configuration and this build's evidence plugins.
fn load_sample(
    input: SampleInput<'_>,
    reference: &Path,
    config_path: &Path,
) -> Result<LoadedSample> {
    Ok(match input {
        SampleInput::Traces(paths) => {
            let inputs = sanger::load_sample(paths, reference, config_path)?;
            LoadedSample {
                loaded: Loaded {
                    config: inputs.config,
                    profile: inputs.profile,
                    reference: inputs.reference,
                    run_context: None,
                },
                reads: LoadedReads::Traces(inputs.traces),
            }
        }
        SampleInput::Evidence(paths) => {
            let inputs =
                input_evidence::load_sample(paths, reference, config_path, plugins::EVIDENCE)?;
            LoadedSample {
                loaded: Loaded {
                    config: inputs.config,
                    profile: inputs.profile,
                    reference: inputs.reference,
                    run_context: None,
                },
                reads: LoadedReads::Recorded(inputs.reads),
            }
        }
    })
}

/// The run context at `path`, when one is given (PROP-0014), with its
/// identity logged.
fn load_run_context(path: Option<&Path>) -> Result<Option<crate::sample_analysis::RunContext>> {
    let Some(path) = path else {
        return Ok(None);
    };
    let context = crate::input::run_context::load(path)?;
    tracing::info!(
        event = "run_context_loaded",
        context_run_id = ?context.run_id,
        sha256 = ?context.sha256,
        reads = context.reads.len(),
        run1_pairs = context.length.runs[0].pairs,
        run2_pairs = context.length.runs[1].pairs,
    );
    Ok(Some(context))
}

/// The read evidence document of one trace, in memory (PROP-0011).
pub(crate) fn evidence_document(trace: &Path, config_path: &Path) -> Result<Vec<u8>> {
    evidence::document(trace, config_path)
}

/// Runs one AB1-to-reference analysis.
pub(crate) fn analyze(trace: &Path, reference: &Path, config_path: &Path) -> Result<()> {
    analyze::run(trace, reference, config_path)
}

/// Runs one analysis of a read from its recorded evidence (PROP-0011).
pub(crate) fn analyze_evidence(
    document: &Path,
    reference: &Path,
    config_path: &Path,
) -> Result<()> {
    analyze::run_evidence(document, reference, config_path)
}

/// Runs one reference-free AB1 basecall operation.
pub(crate) fn basecall(trace: &Path, config_path: &Path) -> Result<()> {
    basecall::run(trace, config_path)
}

/// Runs one reference-free read-evidence operation (PROP-0011).
pub(crate) fn evidence(trace: &Path, config_path: &Path) -> Result<()> {
    evidence::run(trace, config_path)
}

/// Runs one sample consensus operation (PROP-0003).
pub(crate) fn consensus(
    sample_id: &str,
    input: SampleInput<'_>,
    reference: &Path,
    run_context: Option<&Path>,
    config_path: &Path,
) -> Result<()> {
    consensus::run(sample_id, input, reference, run_context, config_path)
}

/// Runs one run context operation over a laboratory batch (PROP-0014), from
/// its list of traces or of read evidence documents (exactly one).
pub(crate) fn run_context(
    run_id: &str,
    (traces, evidence): (Option<&Path>, Option<&Path>),
    reference: &Path,
    config_path: &Path,
) -> Result<()> {
    let list = match (traces, evidence) {
        (Some(list), None) => run_context::List::Traces(list),
        (None, Some(list)) => run_context::List::Evidence(list),
        _ => {
            return Err(crate::error::RunContextError::ReadList {
                line: 0,
                reason: "give exactly one of --traces and --evidence",
            }
            .into());
        }
    };
    run_context::run(run_id, list, reference, config_path)
}

/// Runs one core-only call over consensus sequences.
pub(crate) fn call(
    sample_id: &str,
    sequences: &[PathBuf],
    reference: &Path,
    config_path: &Path,
) -> Result<()> {
    call::run(sample_id, sequences, reference, config_path)
}

/// Derives the notation of one sample's variants document.
pub(crate) fn notation(
    sample_id: &str,
    document: &Path,
    reference: &Path,
    config_path: &Path,
) -> Result<()> {
    notation::run(sample_id, document, reference, config_path)
}

/// Runs one run signal operation (PROP-0009).
pub(crate) fn signal(run_id: &str, traces: &[PathBuf], config_path: &Path) -> Result<()> {
    signal::run(run_id, traces, config_path)
}

/// Runs one sample review operation (PROP-0004).
pub(crate) fn review(
    sample_id: &str,
    input: SampleInput<'_>,
    reference: &Path,
    signal: Option<&Path>,
    run_context: Option<&Path>,
    config_path: &Path,
) -> Result<()> {
    review::run(
        sample_id,
        input,
        reference,
        signal,
        run_context,
        config_path,
    )
}

/// Runs one multi-read sample evidence operation.
pub(crate) fn sample(
    sample_id: &str,
    input: SampleInput<'_>,
    reference: &Path,
    config_path: &Path,
) -> Result<()> {
    sample::run(sample_id, input, reference, config_path)
}

/// Worker threads for a sample's reads: the available parallelism, read once
/// per process.
fn read_workers() -> NonZeroUsize {
    static WORKERS: OnceLock<NonZeroUsize> = OnceLock::new();
    *WORKERS.get_or_init(|| thread::available_parallelism().unwrap_or(NonZeroUsize::MIN))
}

/// Reads processed on worker threads write to the operation log in read
/// order, so its content does not depend on scheduling.
impl ReadScope for OperationLog {
    type Records = Records;

    fn record<T>(&self, work: impl FnOnce() -> T) -> (T, Records) {
        self.capture(work)
    }

    fn replay(&self, records: Records) {
        self.append(records);
    }
}

/// One CLI operation bound to its own append-only operation log.
struct Operation {
    log: OperationLog,
    started: Instant,
    failure_event: &'static str,
}

impl Operation {
    /// Opens `<stem>.log`, records the start event, and fails fast when even
    /// that first record cannot be written.
    fn begin(stem: &str, failure_event: &'static str, announce: impl FnOnce()) -> Result<Self> {
        let log = OperationLog::open(stem)?;
        let started = Instant::now();
        log.in_scope(announce);
        log.check()?;
        Ok(Self {
            log,
            started,
            failure_event,
        })
    }

    /// Runs `body` under the log and records a terminal failure record if it fails.
    fn run(self, body: impl FnOnce(&OperationLog, Instant) -> Result<()>) -> Result<()> {
        match self.log.in_scope(|| body(&self.log, self.started)) {
            Ok(()) => Ok(()),
            Err(error) => Err(self.record_failure(error)),
        }
    }

    /// Records the failure in the stage that was running, keeping both errors
    /// if the record itself cannot be persisted.
    fn record_failure(&self, operation: Error) -> Error {
        // Every operation body enters a stage span before it can fail.
        let stage = self.log.stage().unwrap_or("operation");
        self.log.in_scope(|| {
            tracing::error!(
                event = self.failure_event,
                stage,
                elapsed_ms = self.started.elapsed().as_millis(),
                error = ?operation.to_string(),
            );
        });
        match self.log.sync() {
            Ok(()) => operation,
            Err(logging) => Error::OperationAndLog {
                operation: Box::new(operation),
                logging: Box::new(logging),
            },
        }
    }
}
