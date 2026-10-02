//! End-to-end orchestration for DNA operations.

mod analyze;
mod basecall;
mod input;
mod observation;
mod read;
mod sample;
mod sample_metrics;
mod sample_reads;

use std::path::{Path, PathBuf};
use std::time::Instant;

use crate::error::{Error, Result};
use crate::logger::Logger;

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

/// Records a terminal operation failure without discarding either error.
fn record_failure(
    logger: &mut Logger,
    event: &'static str,
    stage: &'static str,
    started: Instant,
    operation: Error,
) -> Error {
    let logging = logger
        .error(
            module_path!(),
            line!(),
            format_args!(
                "event={event} stage={stage} elapsed_ms={} error={:?}",
                started.elapsed().as_millis(),
                operation.to_string()
            ),
        )
        .and_then(|()| logger.sync());
    match logging {
        Ok(()) => operation,
        Err(logging) => Error::OperationAndLog {
            operation: Box::new(operation),
            logging: Box::new(logging),
        },
    }
}
