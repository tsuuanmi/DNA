//! One placement pass over a sample's reads, from their traces or their
//! recorded evidence (PROP-0011), with reads processed concurrently.

use std::collections::BTreeMap;
use std::num::NonZeroUsize;
use std::panic::{self, AssertUnwindSafe};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;

use crate::config::Config;
use crate::error::{Result, SampleError};
use crate::input::evidence::{Recorded, RecordedRead};
use crate::model::read_observation::{RejectedRead, RejectionCause};
use crate::model::sample_result::ReadRejectionResult;
use crate::model::sanger_read::{self, RejectedSanger, Rejection, SangerRead};
use crate::variant_analysis;
use crate::variant_analysis::observation::{Called, Observation};
use dna_call::model::called_read::CalledRead;
use dna_call::read_call::ReadIdentity;
use dna_kernel::model::reference::Reference;
use dna_kernel::profile::Profile;
use dna_sanger::model::attachment::SangerAttachment;
use dna_sanger::model::sanger::Chromatogram;
use dna_sanger::read_processing::{self, Admission};

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

/// Every read of a sample after one placement pass.
pub(crate) struct PlacedReads {
    /// Admitted reads, in input order.
    pub(crate) called: Vec<CalledRead>,
    /// Sanger evidence of the admitted reads as the documents report it, by
    /// read content identity.
    pub(crate) sanger: BTreeMap<String, SangerRead>,
    /// Rejected reads as the documents report them, in input order.
    pub(crate) rejected: Vec<Rejection>,
    /// The Sanger products of the reads processed from their traces, which
    /// the public sample analysis carries; empty for reads from their
    /// recorded evidence, which keeps no trace products.
    pub(crate) products: Products,
    /// Read-level warnings over the admitted reads.
    pub(crate) warning_total: usize,
}

/// The Sanger products of the reads of a pass over traces.
#[derive(Default)]
pub(crate) struct Products {
    /// Sanger products of the admitted reads, by read content identity.
    pub(crate) reads: BTreeMap<String, SangerAttachment>,
    /// Rejected reads with their causes, in trace order.
    pub(crate) rejected: Vec<RejectedRead>,
}

/// What one read became.
enum Placement {
    Admitted {
        called: Box<CalledRead>,
        sanger: Box<SangerRead>,
        /// The trace's Sanger products, for a read processed from its trace.
        products: Option<Box<SangerAttachment>>,
        warnings: usize,
    },
    Rejected {
        rejection: Box<Rejection>,
        /// The rejection with its Sanger products, for a read processed from
        /// its trace.
        products: Option<Box<RejectedRead>>,
    },
}

/// The reads of one sample pass: traces, or their recorded evidence.
#[derive(Clone, Copy)]
pub(crate) enum Reads<'a> {
    /// Traces, processed by the Sanger plugin.
    Traces(&'a [Chromatogram]),
    /// Read evidence documents the Sanger plugin wrote (PROP-0011).
    Recorded(&'a [Recorded]),
}

impl Reads<'_> {
    const fn len(self) -> usize {
        match self {
            Self::Traces(traces) => traces.len(),
            Self::Recorded(reads) => reads.len(),
        }
    }
}

/// Processes every read independently on up to `workers` threads; a read
/// shorter than one signal window, one with too few callable calls, or one
/// the core cannot place is recorded as rejected and the remaining reads
/// continue. A recorded read is placed exactly as the read of its trace.
///
/// The outcome equals a sequential pass: the first other failure in input
/// order, a typed error or a panic, fails the pass, and only the records of
/// the reads up to it are replayed. Workers stop taking reads after the
/// lowest failed one.
pub(crate) fn place(
    reads: Reads<'_>,
    reference: &Reference,
    config: &Config,
    profile: &Profile,
    scope: &impl ReadScope,
    workers: NonZeroUsize,
) -> Result<PlacedReads> {
    let slots = concurrently(
        reads.len(),
        workers,
        |index| {
            scope.record(|| {
                panic::catch_unwind(AssertUnwindSafe(|| match reads {
                    Reads::Traces(traces) => {
                        place_trace(index, &traces[index], reference, config, profile)
                    }
                    Reads::Recorded(recorded) => {
                        place_recorded(index, &recorded[index], reference, config, profile)
                    }
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
        products: Products::default(),
        warning_total: 0,
    };
    for placement in placements {
        match placement {
            Placement::Admitted {
                called,
                sanger,
                products,
                warnings,
            } => {
                placed.warning_total += warnings;
                placed.sanger.insert(called.input_sha256.clone(), *sanger);
                if let Some(products) = products {
                    placed
                        .products
                        .reads
                        .insert(called.input_sha256.clone(), *products);
                }
                placed.called.push(*called);
            }
            Placement::Rejected {
                rejection,
                products,
            } => {
                placed.rejected.push(*rejection);
                if let Some(products) = products {
                    placed.products.rejected.push(*products);
                }
            }
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

fn place_trace(
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
        let read = RejectedRead {
            input_name: trace.source_name.clone(),
            input_sha256: trace.source_sha256.clone(),
            cause,
        };
        Ok(Placement::Rejected {
            rejection: Box::new(Rejection::of(&read)),
            products: Some(Box::new(read)),
        })
    };
    let processed = match read_processing::admit(trace, &config.sanger)? {
        Admission::Processed(processed) => *processed,
        Admission::TooShort(short) => {
            rejected_too_short(
                index,
                &trace.source_sha256,
                short.calls,
                short.minimum_calls,
            );
            return reject(RejectionCause::TooShort(short));
        }
        Admission::TooFewCallable { sanger, rejection } => {
            rejected_callability(
                index,
                &trace.source_sha256,
                rejection.callable_calls,
                rejection.minimum_callable_calls,
            );
            return reject(RejectionCause::Callability {
                sanger: *sanger,
                rejection,
            });
        }
    };
    let completed =
        match variant_analysis::observation::observe(trace, processed, reference, config, profile)?
        {
            Observation::Placed(completed) => *completed,
            Observation::Unplaced {
                rejection,
                error,
                sanger,
            } => {
                rejected_placement(index, &trace.source_sha256, rejection.label(), &error);
                return reject(RejectionCause::Placement {
                    sanger: *sanger,
                    rejection,
                });
            }
        };
    let read = completed.read;
    completed_read(index, &read.called);
    Ok(Placement::Admitted {
        sanger: Box::new(SangerRead::of(&read.sanger)?),
        called: Box::new(read.called),
        products: Some(Box::new(read.sanger)),
        warnings: completed.warning_total,
    })
}

/// Places one read from its recorded evidence: a read the Sanger plugin
/// rejected is rejected again with its recorded cause, and every other read
/// goes to the core as the read of its trace does.
fn place_recorded(
    index: usize,
    recorded: &Recorded,
    reference: &Reference,
    config: &Config,
    profile: &Profile,
) -> Result<Placement> {
    let (input_name, input_sha256) = match recorded {
        Recorded::Read(read) => (&read.input_name, &read.input_sha256),
        Recorded::Rejected(rejection) => (&rejection.input_name, &rejection.input_sha256),
    };
    tracing::info!(
        event = "sample_read_started",
        read_index = index,
        trace_name = ?input_name,
        trace_sha256 = %input_sha256,
    );
    let read = match recorded {
        Recorded::Rejected(rejection) => {
            match rejection.cause {
                ReadRejectionResult::TooShort {
                    calls,
                    minimum_calls,
                } => rejected_too_short(index, input_sha256, calls, minimum_calls),
                ReadRejectionResult::CallableCallsBelowMinimum {
                    callable_calls,
                    minimum_callable_calls,
                } => rejected_callability(
                    index,
                    input_sha256,
                    callable_calls,
                    minimum_callable_calls,
                ),
                _ => {}
            }
            return Ok(Placement::Rejected {
                rejection: rejection.clone(),
                products: None,
            });
        }
        Recorded::Read(read) => read,
    };
    let RecordedRead {
        input_name,
        input_sha256,
        evidence,
        sanger,
    } = &**read;
    let stage = tracing::info_span!("read_evidence").entered();
    variant_analysis::observation::record_evidence(evidence, config);
    drop(stage);
    let identity = ReadIdentity {
        input_name: input_name.clone(),
        input_sha256: input_sha256.clone(),
        declared: None,
    };
    let called = match variant_analysis::observation::place(
        identity,
        evidence.clone(),
        reference,
        config,
        profile,
    )? {
        Called::Placed(called) => *called,
        Called::Unplaced { rejection, error } => {
            rejected_placement(index, input_sha256, rejection.label(), &error);
            return Ok(Placement::Rejected {
                rejection: Box::new(Rejection {
                    input_name: input_name.clone(),
                    input_sha256: input_sha256.clone(),
                    cause: sanger_read::placement(rejection),
                    sanger: Some(RejectedSanger {
                        integrity: sanger.signal_quality.integrity.clone(),
                        callability: sanger.callability.clone(),
                    }),
                }),
                products: None,
            });
        }
    };
    let warnings = variant_analysis::observation::warning_total(&sanger.read_warnings(), &called);
    completed_read(index, &called);
    Ok(Placement::Admitted {
        called: Box::new(called),
        sanger: Box::new(sanger.clone()),
        products: None,
        warnings,
    })
}

fn rejected_too_short(index: usize, sha256: &str, calls: usize, minimum_calls: usize) {
    tracing::warn!(
        event = "sample_read_rejected",
        read_index = index,
        trace_sha256 = %sha256,
        reason = "too_short",
        calls,
        minimum_calls,
    );
}

fn rejected_callability(
    index: usize,
    sha256: &str,
    callable_calls: usize,
    minimum_callable_calls: usize,
) {
    tracing::warn!(
        event = "sample_read_rejected",
        read_index = index,
        trace_sha256 = %sha256,
        reason = "callable_calls_below_minimum",
        callable_calls,
        minimum_callable_calls,
    );
}

fn rejected_placement(index: usize, sha256: &str, reason: &str, error: &crate::error::Error) {
    tracing::warn!(
        event = "sample_read_rejected",
        read_index = index,
        trace_sha256 = %sha256,
        reason,
        detail = %error,
    );
}

fn completed_read(index: usize, called: &CalledRead) {
    tracing::info!(
        event = "sample_read_completed",
        read_index = index,
        trace_sha256 = %called.input_sha256,
        orientation = ?called.alignment.orientation(),
        segments = called.alignment.reference_segments().len(),
        variants = called.variants.reported.len(),
    );
}

#[cfg(test)]
mod tests {
    use std::sync::Mutex;

    use crate::error::{Error, KernelError};

    use super::*;

    fn workers(count: usize) -> NonZeroUsize {
        NonZeroUsize::new(count).unwrap_or(NonZeroUsize::MIN)
    }

    // SRS-SAMPLE-001
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

    // SRS-SAMPLE-001
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

    // SRS-SAMPLE-001
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

    // SRS-SAMPLE-001
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
