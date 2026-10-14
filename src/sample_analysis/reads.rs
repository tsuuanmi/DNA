//! One placement pass over a sample's traces, with reads processed concurrently.

use std::collections::BTreeMap;
use std::num::NonZeroUsize;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;

use crate::config::Config;
use crate::model::read_observation::{ReadObservation, RejectedRead, RejectionCause};
use crate::variant_analysis;
use dna_core::model::called_read::CalledRead;
use dna_core::read_call::PlacementRejection;
use dna_kernel::error::Result;
use dna_kernel::model::reference::Reference;
use dna_kernel::profile::Profile;
use dna_sanger::model::attachment::{SangerAttachment, SangerRejection};
use dna_sanger::model::sanger::Chromatogram;
use dna_sanger::read_processing;

/// Where the instrumentation of a read processed on a worker thread goes.
///
/// Each read runs inside [`record`](Self::record) on some thread; the caller's
/// thread then hands back what it recorded through [`replay`](Self::replay),
/// in read order and only up to the first read that fails, so what a
/// recording scope keeps is what a sequential pass would have produced.
pub(crate) trait ReadScope: Sync {
    /// What one read recorded.
    type Records: Send;

    /// Runs one read's work on the current thread.
    fn record<T>(&self, work: impl FnOnce() -> T) -> (T, Self::Records);

    /// Hands one read's records back on the caller's thread.
    fn replay(&self, records: Self::Records);
}

/// Every trace of a sample after one placement pass.
pub(crate) struct PlacedReads {
    /// Admitted reads, in trace order.
    pub(crate) called: Vec<CalledRead>,
    /// Sanger evidence of the admitted reads, by read content identity.
    pub(crate) sanger: BTreeMap<String, SangerAttachment>,
    /// Rejected reads, in trace order.
    pub(crate) rejected: Vec<RejectedRead>,
    /// Read-level warnings over the admitted reads.
    pub(crate) warning_total: usize,
}

/// What one trace became.
enum Placement {
    Admitted {
        read: Box<ReadObservation>,
        warnings: usize,
    },
    Rejected(Box<RejectedRead>),
}

/// Processes every trace independently and concurrently; a read with too few
/// callable calls, or one the core cannot place, is recorded as rejected and
/// the remaining reads continue. The first other failure, in trace order,
/// fails the pass.
pub(crate) fn place(
    traces: &[Chromatogram],
    reference: &Reference,
    config: &Config,
    profile: &Profile,
    scope: &impl ReadScope,
) -> Result<PlacedReads> {
    let outcomes = concurrently(traces.len(), |index| {
        scope.record(|| place_one(index, &traces[index], reference, config, profile))
    });
    let mut placed = PlacedReads {
        called: Vec::with_capacity(traces.len()),
        sanger: BTreeMap::new(),
        rejected: Vec::new(),
        warning_total: 0,
    };
    for (outcome, records) in outcomes {
        scope.replay(records);
        match outcome? {
            Placement::Admitted { read, warnings } => {
                let read = *read;
                placed.warning_total += warnings;
                placed
                    .sanger
                    .insert(read.called.input_sha256.clone(), read.sanger);
                placed.called.push(read.called);
            }
            Placement::Rejected(read) => placed.rejected.push(*read),
        }
    }
    Ok(placed)
}

/// Runs `work` for every index in `0..count` on up to the available number
/// of scoped threads, and returns the results in index order.
///
/// A panic in `work` is resumed on the caller's thread once every worker has
/// stopped.
fn concurrently<T: Send>(count: usize, work: impl Fn(usize) -> T + Sync) -> Vec<T> {
    let workers = thread::available_parallelism()
        .map_or(1, NonZeroUsize::get)
        .min(count);
    if workers <= 1 {
        return (0..count).map(work).collect();
    }
    let next = AtomicUsize::new(0);
    let mut indexed = thread::scope(|threads| {
        let handles = (0..workers)
            .map(|_| {
                threads.spawn(|| {
                    let mut done = Vec::new();
                    loop {
                        let index = next.fetch_add(1, Ordering::Relaxed);
                        if index >= count {
                            return done;
                        }
                        done.push((index, work(index)));
                    }
                })
            })
            .collect::<Vec<_>>();
        handles
            .into_iter()
            .flat_map(|handle| {
                handle
                    .join()
                    .unwrap_or_else(|panic| std::panic::resume_unwind(panic))
            })
            .collect::<Vec<_>>()
    });
    indexed.sort_unstable_by_key(|(index, _)| *index);
    indexed.into_iter().map(|(_, value)| value).collect()
}

fn place_one(
    index: usize,
    trace: &Chromatogram,
    reference: &Reference,
    config: &Config,
    profile: &Profile,
) -> Result<Placement> {
    tracing::info!(
        event = "sample_read_started",
        read_index = index,
        trace_name = ?trace.source_name,
        trace_sha256 = %trace.source_sha256,
    );
    let prepared = read_processing::prepare(trace, &config.sanger)?;
    if let Some(rejection) = prepared.rejection(&config.sanger) {
        tracing::warn!(
            event = "sample_read_rejected",
            read_index = index,
            trace_sha256 = %trace.source_sha256,
            reason = "callable_calls_below_minimum",
            callable_calls = rejection.callable_calls,
            minimum_callable_calls = rejection.minimum_callable_calls,
        );
        return Ok(Placement::Rejected(Box::new(RejectedRead {
            sanger: SangerRejection {
                input_name: trace.source_name.clone(),
                input_sha256: trace.source_sha256.clone(),
                integrity: prepared.signal.integrity,
                callability: prepared.callability,
            },
            cause: RejectionCause::Callability(rejection),
        })));
    }
    let processed = read_processing::finish(
        trace,
        prepared,
        &config.sanger,
        config.core.variant_calling.read_end_margin,
    )?;
    let sanger = SangerRejection {
        input_name: trace.source_name.clone(),
        input_sha256: trace.source_sha256.clone(),
        integrity: processed.signal.integrity.clone(),
        callability: processed.callability.clone(),
    };
    let completed = match variant_analysis::observation::observe(
        trace, processed, reference, config, profile,
    ) {
        Ok(completed) => completed,
        Err(error) => {
            let Some(placement) = PlacementRejection::of(&error) else {
                return Err(error);
            };
            tracing::warn!(
                event = "sample_read_rejected",
                read_index = index,
                trace_sha256 = %trace.source_sha256,
                reason = placement.label(),
                detail = %error,
            );
            return Ok(Placement::Rejected(Box::new(RejectedRead {
                sanger,
                cause: RejectionCause::Placement(placement),
            })));
        }
    };
    tracing::info!(
        event = "sample_read_completed",
        read_index = index,
        trace_sha256 = %completed.read.called.input_sha256,
        orientation = ?completed.read.called.alignment.orientation,
        segments = completed.read.called.alignment.reference_segments.len(),
        variants = completed.read.called.variants.reported.len(),
    );
    Ok(Placement::Admitted {
        read: Box::new(completed.read),
        warnings: completed.warning_total,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn concurrent_results_keep_index_order() {
        let squares = concurrently(97, |index| index * index);
        assert_eq!(
            squares,
            (0..97).map(|index| index * index).collect::<Vec<_>>()
        );
        assert!(concurrently(0, |index| index).is_empty());
    }
}
