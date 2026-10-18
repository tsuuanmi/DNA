//! One placement pass over a sample's traces, with reads processed concurrently.

use std::collections::BTreeMap;
use std::num::NonZeroUsize;
use std::panic::{self, AssertUnwindSafe};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;

use crate::config::Config;
use crate::error::{Result, SampleError};
use crate::model::read_observation::{ReadObservation, RejectedRead, RejectionCause};
use crate::variant_analysis;
use crate::variant_analysis::observation::Observation;
use dna_core::model::called_read::CalledRead;
use dna_kernel::model::reference::Reference;
use dna_kernel::profile::Profile;
use dna_sanger::model::attachment::{SangerAttachment, SangerRejection};
use dna_sanger::model::sanger::Chromatogram;
use dna_sanger::model::signal::ShortRead;
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

/// Processes every trace independently on up to `workers` threads; a read
/// shorter than one signal window, one with too few callable calls, or one
/// the core cannot place is recorded as rejected and the remaining reads
/// continue.
///
/// The outcome equals a sequential pass: the first other failure in trace
/// order, a typed error or a panic, fails the pass, and only the records of
/// the reads up to it are replayed. Workers stop taking reads after the
/// lowest failed one.
pub(crate) fn place(
    traces: &[Chromatogram],
    reference: &Reference,
    config: &Config,
    profile: &Profile,
    scope: &impl ReadScope,
    workers: NonZeroUsize,
) -> Result<PlacedReads> {
    let slots = concurrently(
        traces.len(),
        workers,
        |index| {
            scope.record(|| {
                panic::catch_unwind(AssertUnwindSafe(|| {
                    place_one(index, &traces[index], reference, config, profile)
                }))
            })
        },
        |(outcome, _)| !matches!(outcome, Ok(Ok(_))),
    );
    let placements = settle(slots, |records| scope.replay(records))?;
    let mut placed = PlacedReads {
        called: Vec::with_capacity(placements.len()),
        sanger: BTreeMap::new(),
        rejected: Vec::new(),
        warning_total: 0,
    };
    for placement in placements {
        match placement {
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

/// One read's outcome, a panic caught as a value, and what it recorded.
type Slot<T, R> = (thread::Result<Result<T>>, R);

/// Takes outcomes in index order, replaying each one's records, up to the
/// first failure, which it returns or resumes as a sequential pass would.
fn settle<T, R>(slots: Vec<Option<Slot<T, R>>>, mut replay: impl FnMut(R)) -> Result<Vec<T>> {
    let mut values = Vec::with_capacity(slots.len());
    for slot in slots {
        let (outcome, records) = slot.ok_or(SampleError::Inconsistent(
            "a read before the first failure was not processed",
        ))?;
        replay(records);
        match outcome {
            Ok(result) => values.push(result?),
            Err(payload) => panic::resume_unwind(payload),
        }
    }
    Ok(values)
}

/// Runs `work` for the indexes `0..count` on up to `workers` scoped threads
/// (the calling thread when one suffices) and returns the results by index.
///
/// Workers stop taking indexes above the lowest one whose result `failed`;
/// every index up to it has a result, the ones above it may not.
fn concurrently<T: Send>(
    count: usize,
    workers: NonZeroUsize,
    work: impl Fn(usize) -> T + Sync,
    failed: impl Fn(&T) -> bool + Sync,
) -> Vec<Option<T>> {
    let next = AtomicUsize::new(0);
    let first_failed = AtomicUsize::new(usize::MAX);
    let worker = || {
        let mut done = Vec::new();
        loop {
            let index = next.fetch_add(1, Ordering::Relaxed);
            if index >= count || index > first_failed.load(Ordering::Relaxed) {
                return done;
            }
            let value = work(index);
            if failed(&value) {
                first_failed.fetch_min(index, Ordering::Relaxed);
            }
            done.push((index, value));
        }
    };
    let threads = workers.get().min(count);
    let done = if threads <= 1 {
        worker()
    } else {
        thread::scope(|scope| {
            let handles = (0..threads)
                .map(|_| scope.spawn(worker))
                .collect::<Vec<_>>();
            handles
                .into_iter()
                .flat_map(|handle| {
                    handle
                        .join()
                        .unwrap_or_else(|payload| panic::resume_unwind(payload))
                })
                .collect::<Vec<_>>()
        })
    };
    let mut slots = (0..count).map(|_| None).collect::<Vec<_>>();
    for (index, value) in done {
        if let Some(slot) = slots.get_mut(index) {
            *slot = Some(value);
        }
    }
    slots
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
    let reject = |cause| {
        Ok(Placement::Rejected(Box::new(RejectedRead {
            input_name: trace.source_name.clone(),
            input_sha256: trace.source_sha256.clone(),
            cause,
        })))
    };
    let prepared = match read_processing::prepare(trace, &config.sanger) {
        Ok(prepared) => prepared,
        Err(error) => {
            let Some(short) = ShortRead::of(&error) else {
                return Err(error.into());
            };
            tracing::warn!(
                event = "sample_read_rejected",
                read_index = index,
                trace_sha256 = %trace.source_sha256,
                reason = "too_short",
                calls = short.calls,
                minimum_calls = short.minimum_calls,
            );
            return reject(RejectionCause::TooShort(short));
        }
    };
    if let Some(rejection) = prepared.rejection(&config.sanger) {
        tracing::warn!(
            event = "sample_read_rejected",
            read_index = index,
            trace_sha256 = %trace.source_sha256,
            reason = "callable_calls_below_minimum",
            callable_calls = rejection.callable_calls,
            minimum_callable_calls = rejection.minimum_callable_calls,
        );
        let sanger = SangerRejection {
            integrity: prepared.signal.integrity,
            callability: prepared.callability,
        };
        return reject(RejectionCause::Callability { sanger, rejection });
    }
    let processed = read_processing::finish(
        trace,
        prepared,
        &config.sanger,
        config.core.variant_calling.read_end_margin,
    )?;
    let completed =
        match variant_analysis::observation::observe(trace, processed, reference, config, profile)?
        {
            Observation::Placed(completed) => *completed,
            Observation::Unplaced {
                rejection,
                error,
                sanger,
            } => {
                tracing::warn!(
                    event = "sample_read_rejected",
                    read_index = index,
                    trace_sha256 = %trace.source_sha256,
                    reason = rejection.label(),
                    detail = %error,
                );
                return reject(RejectionCause::Placement {
                    sanger: *sanger,
                    rejection,
                });
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
    use std::sync::Mutex;

    use crate::error::{Error, KernelError};

    use super::*;

    fn workers(count: usize) -> NonZeroUsize {
        NonZeroUsize::new(count).unwrap_or(NonZeroUsize::MIN)
    }

    #[test]
    fn concurrent_results_keep_index_order() {
        for count in [1, 4, 16] {
            let squares = concurrently(97, workers(count), |index| index * index, |_| false);
            assert_eq!(
                squares,
                (0..97).map(|index| Some(index * index)).collect::<Vec<_>>()
            );
        }
        assert!(concurrently(0, workers(4), |index| index, |_| false).is_empty());
    }

    #[test]
    fn workers_stop_after_the_lowest_failure_and_keep_every_result_below_it() {
        for count in [1, 4, 16] {
            let ran = AtomicUsize::new(0);
            let slots = concurrently(
                200,
                workers(count),
                |index| {
                    ran.fetch_add(1, Ordering::Relaxed);
                    index
                },
                |index| *index == 5,
            );
            assert!(slots[..=5].iter().all(Option::is_some));
            if count == 1 {
                assert_eq!(ran.load(Ordering::Relaxed), 6);
                assert!(slots[6..].iter().all(Option::is_none));
            }
        }
    }

    /// Reads 0 and 2 succeed; `failures` fail the others: `Some(true)` with a
    /// panic, `Some(false)` with a typed error.
    fn run(
        threads: usize,
        failures: [Option<bool>; 4],
    ) -> (thread::Result<Result<Vec<usize>>>, Vec<usize>) {
        let slots = concurrently(
            4,
            workers(threads),
            |index| {
                let outcome = panic::catch_unwind(|| match failures[index] {
                    Some(true) => panic!("read {index} panicked"),
                    Some(false) => Err(Error::from(SampleError::NoReads)),
                    None => Ok(index),
                });
                (outcome, index)
            },
            |(outcome, _)| !matches!(outcome, Ok(Ok(_))),
        );
        let replayed = Mutex::new(Vec::new());
        let settled = panic::catch_unwind(AssertUnwindSafe(|| {
            settle(slots, |records| {
                replayed
                    .lock()
                    .unwrap_or_else(std::sync::PoisonError::into_inner)
                    .push(records);
            })
        }));
        let replayed = replayed
            .into_inner()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        (settled, replayed)
    }

    #[test]
    fn an_earlier_error_wins_over_a_later_panic() {
        for threads in [1, 4] {
            let (settled, replayed) = run(threads, [None, Some(false), None, Some(true)]);
            assert!(matches!(
                settled,
                Ok(Err(Error::Kernel(KernelError::Sample(
                    SampleError::NoReads
                ))))
            ));
            assert_eq!(replayed, [0, 1]);
        }
    }

    #[test]
    fn an_earlier_panic_wins_over_a_later_error() {
        for threads in [1, 4] {
            let (settled, replayed) = run(threads, [None, Some(true), None, Some(false)]);
            assert!(settled.is_err());
            assert_eq!(replayed, [0, 1]);
        }
    }

    #[test]
    fn every_outcome_is_kept_without_a_failure() {
        let (settled, replayed) = run(4, [None; 4]);
        assert!(matches!(settled, Ok(Ok(values)) if values == [0, 1, 2, 3]));
        assert_eq!(replayed, [0, 1, 2, 3]);
    }
}
