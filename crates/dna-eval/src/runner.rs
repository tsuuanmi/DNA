//! Runs DNA over the panel's samples on a bounded pool of threads.
//!
//! Each sample runs the recommended `consensus` then `call` workflow in
//! process: [`dna::sample_analysis::analyze_sanger`] with one worker, since
//! the pool already runs samples in parallel, then
//! [`SampleAnalysis::call`](dna::sample_analysis::SampleAnalysis::call). Repeated trace content is
//! removed first, keeping the first path of each content in manifest order.
//! A sample that fails, or a DNA panic, yields a [`Failure`]; other samples
//! continue. Outcomes are returned in sample order, whatever the thread count.

use std::collections::BTreeSet;
use std::fs::File;
use std::io;
use std::panic::{self, AssertUnwindSafe};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;
use std::time::Instant;

use std::num::NonZeroUsize;

use dna::sample_analysis::{self, SampleAnalysis, SampleOutcome};
use sha2::{Digest, Sha256};

use crate::cause::loci::{self, Loci, Margins};
use crate::compare::{Allele, Failure};
use crate::error::{Error, Result};
use crate::evidence::{self, Evidence};
use crate::lod;
use crate::manifest::Sample;

/// Largest accepted thread count.
pub(crate) const MAX_THREADS: usize = 256;

/// What DNA returned for one sample.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Called {
    /// Traces set aside before the consensus.
    pub(crate) rejected: usize,
    /// The profile-notation calls as alleles.
    pub(crate) alleles: Vec<Allele>,
    /// Callable reference intervals of the consensus, closed and 1-based.
    pub(crate) intervals: Vec<(i64, i64)>,
    /// Where some read shows a second base at the limit of detection.
    pub(crate) signals: lod::Signals,
    /// What the reads show at each region position, for the cause of a
    /// discrepancy.
    pub(crate) loci: Loci,
    /// Identities of the configuration, reference, and profile.
    pub(crate) identity: Identity,
}

/// Scientific identities DNA reports for a sample.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, serde::Serialize)]
pub(crate) struct Identity {
    pub(crate) configuration_sha256: String,
    pub(crate) reference_name: String,
    pub(crate) reference_sha256: String,
    pub(crate) profile_id: String,
    pub(crate) profile_sha256: String,
}

/// One sample's run.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Outcome {
    /// Distinct traces given to DNA, once duplicates are removed.
    pub(crate) traces: Option<usize>,
    pub(crate) result: std::result::Result<Called, Failure>,
    /// The draft triage evidence, when the run asks for it and DNA analysed
    /// the sample or rejected every read; it is kept when the call fails. A
    /// failure to take it leaves the call as it is.
    pub(crate) triage: Option<Evidence>,
    /// Wall time of the sample, in milliseconds.
    pub(crate) millis: u128,
}

/// A worker count within `1..=MAX_THREADS`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Threads(usize);

impl Threads {
    pub(crate) const fn new(threads: usize) -> Result<Self> {
        if threads == 0 || threads > MAX_THREADS {
            return Err(Error::Threads {
                threads,
                maximum: MAX_THREADS,
            });
        }
        Ok(Self(threads))
    }

    pub(crate) const fn get(self) -> usize {
        self.0
    }
}

/// Runs DNA over every sample with `threads` workers on `reference`, the
/// FASTA path and the SHA-256 of the sequence the harness compares on (a
/// sample DNA analyses on another sequence fails), with DNA's configuration
/// and the margins the harness reads from it, taking each sample's
/// draft triage evidence when `triage` is set: from the analysis, even when
/// its call then fails, or from the rejected reads when none is admitted.
pub(crate) fn run(
    samples: &[Sample],
    (reference, reference_sha256): (&Path, &str),
    (config, margins): (&Path, Margins),
    threads: Threads,
    triage: bool,
) -> Result<Vec<Outcome>> {
    run_with(samples, threads, |sample, traces| {
        let failure = |error: dna::error::Error| Failure::Dna(error.to_string());
        let analysis = match sample_analysis::analyze_sanger_outcome(
            traces,
            reference,
            config,
            NonZeroUsize::MIN,
        ) {
            Ok(SampleOutcome::Analyzed(analysis)) => analysis,
            Ok(SampleOutcome::Rejected(rejected)) => {
                return (
                    Err(Failure::Rejected {
                        reads: rejected.len(),
                    }),
                    triage.then(|| evidence::rejected_only(rejected.values())),
                );
            }
            Err(error) => return (Err(failure(error)), None),
        };
        if analysis.reference.sha256 != reference_sha256 {
            return (
                Err(Failure::Reference {
                    analysed: analysis.reference.sha256.clone(),
                    compared: reference_sha256.to_owned(),
                }),
                None,
            );
        }
        let evidence = triage.then(|| evidence::extract(&analysis));
        (call(&analysis, &sample.id, margins), evidence)
    })
}

/// Calls an analysed sample's consensus and converts the calls into alleles
/// and the callable segments into intervals.
fn call(
    analysis: &SampleAnalysis,
    sample_id: &str,
    margins: Margins,
) -> std::result::Result<Called, Failure> {
    let result = analysis
        .call(sample_id)
        .map_err(|error| Failure::Dna(error.to_string()))?;
    let alleles = result
        .notation
        .ok_or(Failure::Notation)?
        .iter()
        .map(|call| {
            let position = if call.ordinal == 0 {
                call.position.to_string()
            } else {
                format!("{}.{}", call.position, call.ordinal)
            };
            let base = |base: Option<char>| base.map_or_else(|| "-".to_owned(), String::from);
            Allele::new(&position, &base(call.reference), &base(call.alternate))
        })
        .collect();
    // Half-open 0-based segments become closed 1-based intervals.
    let intervals = result
        .callable_segments
        .iter()
        .filter(|segment| segment.end_0based_exclusive > segment.start_0based)
        .filter_map(|segment| {
            Some((
                i64::try_from(segment.start_0based).ok()? + 1,
                i64::try_from(segment.end_0based_exclusive).ok()?,
            ))
        })
        .collect();
    let signals = lod::signals(analysis).map_err(|error| Failure::Dna(error.to_string()))?;
    let loci = loci::extract(analysis, margins).map_err(|error| Failure::Dna(error.to_string()))?;
    Ok(Called {
        rejected: analysis.rejected.len(),
        alleles,
        intervals,
        signals,
        loci,
        identity: Identity {
            configuration_sha256: analysis.evidence.configuration_sha256.clone(),
            reference_name: analysis.reference.name.clone(),
            reference_sha256: analysis.reference.sha256.clone(),
            profile_id: analysis.profile.id.clone(),
            profile_sha256: analysis.profile.sha256.clone(),
        },
    })
}

/// What running DNA on one sample gives: its calls or failure, and its draft
/// triage evidence.
pub(crate) type Ran = (std::result::Result<Called, Failure>, Option<Evidence>);

/// Runs `call` over every sample's distinct traces with `threads` workers.
pub(crate) fn run_with<F>(samples: &[Sample], threads: Threads, call: F) -> Result<Vec<Outcome>>
where
    F: Fn(&Sample, &[PathBuf]) -> Ran + Sync,
{
    let next = AtomicUsize::new(0);
    let done = AtomicUsize::new(0);
    let step = samples.len().div_ceil(20).max(1);
    let work = || {
        let mut outcomes = Vec::new();
        loop {
            let index = next.fetch_add(1, Ordering::Relaxed);
            let Some(sample) = samples.get(index) else {
                break;
            };
            outcomes.push((index, one(sample, &call)));
            let finished = done.fetch_add(1, Ordering::Relaxed) + 1;
            if finished.is_multiple_of(step) || finished == samples.len() {
                eprintln!("dna-eval: {finished}/{} samples", samples.len());
            }
        }
        outcomes
    };
    let mut indexed = thread::scope(|scope| {
        let workers = (0..threads.get().min(samples.len()))
            .map(|_| scope.spawn(work))
            .collect::<Vec<_>>();
        workers
            .into_iter()
            .map(|worker| worker.join().map_err(|_| Error::Worker))
            .collect::<Result<Vec<_>>>()
    })?
    .into_iter()
    .flatten()
    .collect::<Vec<_>>();
    indexed.sort_by_key(|(index, _)| *index);
    Ok(indexed.into_iter().map(|(_, outcome)| outcome).collect())
}

fn one<F>(sample: &Sample, call: &F) -> Outcome
where
    F: Fn(&Sample, &[PathBuf]) -> Ran,
{
    let started = Instant::now();
    let (traces, (result, triage)) = match distinct(&sample.traces) {
        Ok(traces) => {
            let ran = panic::catch_unwind(AssertUnwindSafe(|| call(sample, &traces)))
                .unwrap_or((Err(Failure::Panic), None));
            (Some(traces.len()), ran)
        }
        Err(failure) => (None, (Err(failure), None)),
    };
    Outcome {
        traces,
        result,
        triage,
        millis: started.elapsed().as_millis(),
    }
}

/// The traces with distinct content, in their first order.
fn distinct(traces: &[PathBuf]) -> std::result::Result<Vec<PathBuf>, Failure> {
    let mut seen = BTreeSet::new();
    let mut kept = Vec::with_capacity(traces.len());
    for path in traces {
        let digest = sha256_file(path).map_err(|error| Failure::Trace {
            path: path.clone(),
            message: error.to_string(),
        })?;
        if seen.insert(digest) {
            kept.push(path.clone());
        }
    }
    Ok(kept)
}

/// The SHA-256 of a file, read in a stream with bounded memory.
pub(crate) fn sha256_file(path: &Path) -> io::Result<String> {
    let mut hasher = Sha256::new();
    io::copy(&mut File::open(path)?, &mut hasher)?;
    Ok(format!("{:x}", hasher.finalize()))
}

#[cfg(test)]
mod tests {
    use std::fs;

    use super::*;
    use crate::manifest::Weight;

    fn sample(id: &str, traces: Vec<PathBuf>) -> Sample {
        let weight = Weight {
            value: 1.0,
            text: "1".into(),
        };
        Sample {
            id: id.into(),
            batch: "b1".into(),
            stratum: "clean".into(),
            weight: weight.clone(),
            weight_corpus: weight,
            anchor: false,
            traces,
            retired: false,
            truth_digest: None,
            fields: Vec::new(),
        }
    }

    fn identity() -> Identity {
        Identity {
            configuration_sha256: "c".into(),
            reference_name: "r".into(),
            reference_sha256: "s".into(),
            profile_id: "p".into(),
            profile_sha256: "q".into(),
        }
    }

    #[test]
    fn removes_repeated_trace_content_and_keeps_order()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        let directory = tempfile::tempdir()?;
        let paths = ["b", "a", "copy-of-b"].map(|name| directory.path().join(name));
        fs::write(&paths[0], b"trace one")?;
        fs::write(&paths[1], b"trace two")?;
        fs::write(&paths[2], b"trace one")?;
        assert_eq!(distinct(&paths)?, [paths[0].clone(), paths[1].clone()]);
        let missing = directory.path().join("missing");
        assert!(matches!(distinct(&[missing]), Err(Failure::Trace { .. })));
        Ok(())
    }

    #[test]
    fn outcomes_follow_sample_order_for_any_thread_count()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        let directory = tempfile::tempdir()?;
        let trace = directory.path().join("trace");
        fs::write(&trace, b"trace")?;
        let samples = (0..40)
            .map(|index| sample(&format!("S{index:02}"), vec![trace.clone(), trace.clone()]))
            .collect::<Vec<_>>();
        let called = |sample: &Sample, traces: &[PathBuf]| {
            let index = sample.id[1..]
                .parse::<usize>()
                .map_err(|_| Failure::Panic)?;
            match index % 3 {
                0 => Err(Failure::Dna(format!("fails with {} traces", traces.len()))),
                1 => panic!("a DNA defect"),
                _ => Ok(Called {
                    rejected: 0,
                    alleles: vec![Allele::new(&sample.id[1..], "A", "G")],
                    intervals: Vec::new(),
                    signals: lod::Signals::default(),
                    loci: Loci::default(),
                    identity: identity(),
                }),
            }
        };
        let call = |sample: &Sample, traces: &[PathBuf]| -> Ran { (called(sample, traces), None) };
        let strip = |outcomes: Vec<Outcome>| {
            outcomes
                .into_iter()
                .map(|outcome| (outcome.traces, outcome.result))
                .collect::<Vec<_>>()
        };
        let one = strip(run_with(&samples, Threads::new(1)?, call)?);
        let many = strip(run_with(&samples, Threads::new(7)?, call)?);
        assert_eq!(one, many);
        assert_eq!(one[0].1, Err(Failure::Dna("fails with 1 traces".into())));
        assert_eq!(one[1].1, Err(Failure::Panic));
        assert_eq!(
            one[2].1.as_ref().map(|called| called.alleles.clone()),
            Ok(vec![Allele::new("02", "A", "G")])
        );
        assert_eq!(one[0].0, Some(1));
        Ok(())
    }

    #[test]
    fn bounds_the_thread_count() {
        assert!(matches!(Threads::new(0), Err(Error::Threads { .. })));
        assert!(matches!(
            Threads::new(MAX_THREADS + 1),
            Err(Error::Threads { .. })
        ));
        assert!(Threads::new(MAX_THREADS).is_ok());
    }
}
