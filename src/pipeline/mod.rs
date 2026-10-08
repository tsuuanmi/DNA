//! End-to-end orchestration for DNA operations.

mod analyze;
mod basecall;
mod path;
mod sample;
mod sample_metrics;
mod sample_reads;

use std::path::{Path, PathBuf};
use std::time::Instant;

use crate::error::{Error, Result};
use crate::operation_log::OperationLog;

/// Runs one AB1-to-reference analysis.
pub(crate) fn analyze(trace: &Path, reference: &Path, config_path: &Path) -> Result<()> {
    analyze::run(trace, reference, config_path)
}

/// Runs one reference-free AB1 basecall operation.
pub(crate) fn basecall(trace: &Path, config_path: &Path) -> Result<()> {
    basecall::run(trace, config_path)
}

/// Runs one multi-read sample evidence operation.
pub(crate) fn sample(
    sample_id: &str,
    traces: &[PathBuf],
    reference: &Path,
    config_path: &Path,
) -> Result<()> {
    sample::run(sample_id, traces, reference, config_path)
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
