//! Runs DNA over the panel's samples on a bounded pool of threads.
//!
//! Each sample runs the recommended `consensus` then `call` workflow in
//! process: [`dna::sample_analysis::analyze_sanger`] with one worker, since
//! the pool already runs samples in parallel, then
//! [`SampleAnalysis::call`](dna::sample_analysis::SampleAnalysis::call). Repeated trace content is
//! removed first, keeping the first path of each content in manifest order.
//! A sample that fails, or a DNA panic, yields a [`Failure`]; other samples
//! continue. Outcomes are returned in sample order, whatever the thread count.

use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io;
use std::panic::{self, AssertUnwindSafe};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;
use std::time::Instant;

use std::num::NonZeroUsize;

use dna::recorded::{self, EvidenceOutcome};
use dna::sample_analysis::RunContext;
use dna::sample_analysis::{
    self, Consensus, ConventionLabel, EmpopExport, NotationCall, RunSignal, SampleAnalysis,
    SampleOutcome, TractNote,
};
use dna_consensus::model::consensus::SiteState;
use sha2::{Digest, Sha256};

use crate::cache::Caches;
use crate::cause::loci::{self, Loci};
use crate::compare::{Allele, Failure};
use crate::components;
use crate::error::{Error, Result};
use crate::lod;
use crate::manifest::Sample;
use crate::tracts;

/// Largest accepted thread count.
pub(crate) const MAX_THREADS: usize = 256;

/// What DNA returned for one sample.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Called {
    /// Traces set aside before the consensus.
    pub(crate) rejected: usize,
    /// The profile-notation calls as alleles.
    pub(crate) alleles: Vec<Allele>,
    /// The calls a profile convention labels `not_recorded` (SRS-NOM-017).
    pub(crate) not_recorded: BTreeSet<Allele>,
    /// Callable reference intervals of the consensus, closed and 1-based.
    pub(crate) intervals: Vec<(i64, i64)>,
    /// Where some read shows a second base at the limit of detection.
    pub(crate) signals: lod::Signals,
    /// What the reads show at each region position, for the cause of a
    /// discrepancy.
    pub(crate) loci: Loci,
    /// Identities of the configuration, reference, and profile.
    pub(crate) identity: Identity,
    /// The sample's EMPOP query, when the profile declares one.
    pub(crate) empop: Option<Empop>,
    /// The review list's check items, each `first..=last` (1-based), as
    /// `dna review` builds them.
    pub(crate) checks: Vec<(i64, i64)>,
    /// The same with each item's label names, space-separated.
    pub(crate) check_labels: Vec<(i64, i64, String)>,
    /// The review list's noted length tracts (`tract_length`).
    pub(crate) notes: Vec<TractNote>,
    /// Every read's phase-mixture components (ADR-0082).
    pub(crate) components: Vec<components::Row>,
    /// Every read's tract decodes (PROP-0012).
    pub(crate) tracts: Vec<tracts::Row>,
    /// Whether the run context of the sample's batch (PROP-0014), when it
    /// had one, calibrated τ for run 1 and run 2: it does not where the
    /// batch's pairs, the sample's own out, are too few.
    pub(crate) length_context: Option<[bool; 2]>,
    /// The stretches at a region's ends DNA did not analyse, 1-based and
    /// inclusive (SRS-CONS-012).
    pub(crate) not_analysed: Vec<(i64, i64)>,
}

/// A sample's EMPOP query (PROP-0007), and the consensus decisions its
/// ranges rest on.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(crate) struct Empop {
    /// `16024-16365 73-340`.
    pub(crate) ranges: String,
    /// `73G 263G 309.1C`.
    pub(crate) mitotype: String,
    /// Positions inside the ranges the consensus left undecided (`N`, or no
    /// call): none, if the export is sound. Where a call changes the length,
    /// its `N`s less the positions the ranges leave out.
    pub(crate) undecided: usize,
    /// Positions inside the ranges at a tie the consensus resolved: DNA's
    /// decision, exported as such.
    pub(crate) ties: usize,
}

impl Empop {
    fn of(export: &EmpopExport, consensus: &Consensus) -> Self {
        let inside = |position: usize| {
            export
                .ranges
                .iter()
                .any(|&[start, end]| (start..=end).contains(&position))
        };
        let mut empop = Self {
            ranges: export.ranges_text(),
            mitotype: export.mitotype_text(),
            ..Self::default()
        };
        for site in &consensus.sites {
            let positions = (site.first_0based + 1..=site.last_0based + 1).collect::<Vec<_>>();
            let within = positions
                .iter()
                .filter(|position| inside(**position))
                .count();
            match site.call.as_deref() {
                None => empop.undecided += within,
                // A call of the reference's length reads position by position.
                Some(call) if call.len() == positions.len() => {
                    empop.undecided += call
                        .bytes()
                        .zip(&positions)
                        .filter(|(base, position)| *base == b'N' && inside(**position))
                        .count();
                    if site.state == SiteState::Contested {
                        empop.ties += within;
                    }
                }
                // A length change places its `N`s nowhere in particular: the
                // ranges must leave out at least as many of its positions.
                Some(call) => {
                    let unresolved = call.bytes().filter(|base| *base == b'N').count();
                    empop.undecided += unresolved.saturating_sub(positions.len() - within);
                    if site.state == SiteState::Contested {
                        empop.ties += within;
                    }
                }
            }
        }
        empop
    }
}

/// Scientific identities DNA reports for a sample.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, serde::Serialize, serde::Deserialize)]
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

/// How samples are grouped into runs for their run signal (PROP-0009): the
/// panel's batch when it has at least `batch_reads` reads, else the batch
/// month; or the month alone.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Grouping {
    Batch { batch_reads: usize },
    Month,
}

/// Each batch's run context (`--run-context`, PROP-0014), by batch.
pub(crate) type Contexts = BTreeMap<String, RunContext>;

/// Reads `<directory>/<batch>.run_context.json` for every batch of
/// `samples` that has one, as `dna run-context` writes them over the
/// batches' full folders: a panel holds too few of a batch's samples to
/// build its context in process. With them, their digest: the SHA-256 of
/// one `<batch>\t<file SHA-256>` line per context, in batch order, so that
/// two runs on other contexts are told apart.
pub(crate) fn contexts(directory: &Path, samples: &[Sample]) -> Result<(Contexts, String)> {
    let mut contexts = Contexts::new();
    let mut files = BTreeMap::<String, String>::new();
    for batch in samples.iter().map(|sample| &sample.batch) {
        if contexts.contains_key(batch) {
            continue;
        }
        let path = directory.join(format!("{batch}.run_context.json"));
        if !path.is_file() {
            continue;
        }
        let context = RunContext::read(&path).map_err(|error| Error::RunContext {
            path: path.clone(),
            message: error.to_string(),
        })?;
        let sha256 = sha256_file(&path).map_err(|error| Error::RunContext {
            path: path.clone(),
            message: error.to_string(),
        })?;
        files.insert(batch.clone(), sha256);
        if context.run_id != *batch {
            return Err(Error::RunContext {
                path,
                message: format!("its run_id is {:?}, not the batch", context.run_id),
            });
        }
        contexts.insert(batch.clone(), context);
    }
    let mut digest = Sha256::new();
    for (batch, sha256) in &files {
        digest.update(format!("{batch}\t{sha256}\n"));
    }
    Ok((contexts, format!("{:x}", digest.finalize())))
}

/// A sample's run signal groups: its batch's and its month's.
fn groups(sample: &Sample) -> [String; 2] {
    let month = sample
        .batch_month
        .map_or_else(|| "unknown".to_owned(), |month| month.to_string());
    [format!("batch:{}", sample.batch), format!("month:{month}")]
}

/// What a run of DNA over a panel needs: the reference FASTA and the SHA-256
/// of the sequence the harness compares on (a sample DNA analyses on another
/// sequence fails), DNA's configuration, the target profile's reportable
/// regions (1-based, for the stretches DNA does not analyse), the workers,
/// the run signal grouping, the read evidence cache, and the batches' run
/// contexts.
pub(crate) struct Setup<'s> {
    pub(crate) reference: &'s Path,
    pub(crate) reference_sha256: &'s str,
    pub(crate) config: &'s Path,
    pub(crate) regions: &'s [[usize; 2]],
    pub(crate) threads: Threads,
    pub(crate) grouping: Option<Grouping>,
    pub(crate) cache: Option<&'s Caches>,
    pub(crate) contexts: &'s Contexts,
}

/// What a run signal's first pass built: one `signals.tsv` line per group,
/// and the samples whose analysis failed or rejected every read, which add
/// nothing to any signal.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct SignalPass {
    pub(crate) groups: Vec<String>,
    pub(crate) failed: usize,
}

/// Runs DNA over every sample as [`run_once`] does; with a grouping, first
/// builds each group's run signal from every sample's analysis, then runs
/// again with each sample's review judged against its group's signal.
pub(crate) fn run(
    samples: &[Sample],
    setup: &Setup<'_>,
) -> Result<(Vec<Outcome>, Option<SignalPass>)> {
    let Some(grouping) = setup.grouping else {
        return Ok((run_once(samples, setup, &|_| None)?, None));
    };
    let signals = std::sync::Mutex::new(BTreeMap::<String, RunSignal>::new());
    let failed = AtomicUsize::new(0);
    run_with(samples, setup.threads, |sample, traces| {
        match analyse(traces, setup.reference, setup.config, setup.cache) {
            Ok(Analysed::Analyzed(analysis)) => {
                let (reads, contributions) = analysis.signal_contributions();
                let count = reads.len();
                let mut signals = signals
                    .lock()
                    .unwrap_or_else(std::sync::PoisonError::into_inner);
                for group in groups(sample) {
                    signals
                        .entry(group.clone())
                        .or_insert_with(|| analysis.run_signal(group))
                        .include(reads.clone(), count, &contributions);
                }
            }
            Ok(Analysed::Rejected { .. }) | Err(_) => {
                failed.fetch_add(1, Ordering::Relaxed);
            }
        }
        Err(Failure::Panic)
    })?;
    let signals = signals
        .into_inner()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    let pass = SignalPass {
        groups: signals
            .iter()
            .map(|(group, signal)| {
                let counts = signal
                    .keys
                    .values()
                    .map(|shares| shares.values().sum::<u64>())
                    .collect::<Vec<_>>();
                let at_least =
                    |minimum: u64| counts.iter().filter(|&&calls| calls >= minimum).count();
                format!(
                    "{group}\t{}\t{}\t{}\t{}\t{}\t{}",
                    signal.reads,
                    signal.calls,
                    counts.len(),
                    at_least(20),
                    at_least(40),
                    at_least(100)
                )
            })
            .collect(),
        failed: failed.into_inner(),
    };
    let choose = |sample: &Sample| {
        let [batch, month] = groups(sample);
        match grouping {
            Grouping::Batch { batch_reads } => signals
                .get(&batch)
                .filter(|signal| signal.reads >= batch_reads)
                .or_else(|| signals.get(&month)),
            Grouping::Month => signals.get(&month),
        }
    };
    Ok((run_once(samples, setup, &choose)?, Some(pass)))
}

/// Runs DNA over every sample once, with its consensus calibrated by its
/// batch's run context, if any, and each review judged against the run
/// signal `signal` chooses for it, if any.
fn run_once<'s>(
    samples: &[Sample],
    setup: &Setup<'_>,
    signal: &(dyn Fn(&Sample) -> Option<&'s RunSignal> + Sync),
) -> Result<Vec<Outcome>> {
    run_with(samples, setup.threads, |sample, traces| {
        let analysis = match analyse(traces, setup.reference, setup.config, setup.cache) {
            Ok(Analysed::Analyzed(analysis)) => analysis,
            Ok(Analysed::Rejected { reads }) => return Err(Failure::Rejected { reads }),
            Err(failure) => return Err(failure),
        };
        let analysis = match setup.contexts.get(&sample.batch) {
            Some(context) => match analysis.with_run_context(context.clone()) {
                Ok(analysis) => Box::new(analysis),
                Err(error) => return Err(Failure::Dna(error.to_string())),
            },
            None => analysis,
        };
        if analysis.reference.sha256 != setup.reference_sha256 {
            return Err(Failure::Reference {
                analysed: analysis.reference.sha256.clone(),
                compared: setup.reference_sha256.to_owned(),
            });
        }
        call(&analysis, &sample.id, signal(sample), setup.regions)
    })
}

/// Calls an analysed sample's consensus and converts the calls into alleles
/// and the callable segments into intervals.
fn call(
    analysis: &SampleAnalysis,
    sample_id: &str,
    signal: Option<&RunSignal>,
    regions: &[[usize; 2]],
) -> std::result::Result<Called, Failure> {
    let result = analysis
        .call(sample_id)
        .map_err(|error| Failure::Dna(error.to_string()))?;
    let notation = result.notation.as_ref().ok_or(Failure::Notation)?;
    let allele = |call: &NotationCall| {
        let position = if call.ordinal == 0 {
            call.position.to_string()
        } else {
            format!("{}.{}", call.position, call.ordinal)
        };
        let base = |base: Option<char>| base.map_or_else(|| "-".to_owned(), String::from);
        Allele::new(&position, &base(call.reference), &base(call.alternate))
    };
    let alleles = notation.iter().map(allele).collect();
    let not_recorded = notation
        .iter()
        .filter(|call| {
            call.convention
                .as_ref()
                .is_some_and(|convention| convention.label == ConventionLabel::NotRecorded)
        })
        .map(allele)
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
    let loci = loci::extract(analysis).map_err(|error| Failure::Dna(error.to_string()))?;
    let empop = result
        .empop
        .as_ref()
        .map(|export| Empop::of(export, &analysis.consensus));
    let review = analysis
        .review_list(&result, signal)
        .map_err(|error| Failure::Dna(error.to_string()))?;
    let check_labels = review
        .checks
        .iter()
        .filter_map(|check| {
            Some((
                i64::try_from(check.first).ok()?,
                i64::try_from(check.last).ok()?,
                check.labels.join(" "),
            ))
        })
        .collect::<Vec<_>>();
    let checks = check_labels
        .iter()
        .map(|&(first, last, _)| (first, last))
        .collect();
    Ok(Called {
        rejected: analysis
            .evidence
            .rejected_reads
            .iter()
            .map(|read| read.input_sha256.as_str())
            .collect::<BTreeSet<_>>()
            .len(),
        empop,
        alleles,
        not_recorded,
        intervals,
        signals,
        loci,
        checks,
        check_labels,
        notes: review.notes,
        components: components::rows(analysis),
        tracts: tracts::rows(analysis),
        not_analysed: dna_consensus::decided::not_analysed(&analysis.consensus, regions)
            .iter()
            .filter_map(|stretch| {
                Some((
                    i64::try_from(stretch.first_0based).ok()? + 1,
                    i64::try_from(stretch.last_0based).ok()? + 1,
                ))
            })
            .collect(),
        length_context: analysis
            .run_context()
            .map(|used| used.scale.map(|scale| scale.is_some())),
        identity: Identity {
            configuration_sha256: analysis.evidence.configuration_sha256.clone(),
            reference_name: analysis.reference.name.clone(),
            reference_sha256: analysis.reference.sha256.clone(),
            profile_id: analysis.profile.id.clone(),
            profile_sha256: analysis.profile.sha256.clone(),
        },
    })
}

/// What running DNA on one sample gives: its calls or failure.
pub(crate) type Ran = std::result::Result<Called, Failure>;

/// What a sample's analysis yields: the analysis, or how many distinct
/// reads were set aside.
pub(crate) enum Analysed {
    Analyzed(Box<SampleAnalysis>),
    Rejected { reads: usize },
}

/// Analyses a sample from its traces, or, with an evidence cache, from
/// their read evidence documents, read from the cache or made and stored:
/// every read then comes from its document, whether the cache held it or
/// not, so a hit and a miss analyse alike.
pub(crate) fn analyse(
    traces: &[PathBuf],
    reference: &Path,
    config: &Path,
    cache: Option<&Caches>,
) -> std::result::Result<Analysed, Failure> {
    let failure = |error: dna::error::Error| Failure::Dna(error.to_string());
    let Some(cache) = cache else {
        return match sample_analysis::analyze_sanger_outcome(
            traces,
            reference,
            config,
            NonZeroUsize::MIN,
        )
        .map_err(failure)?
        {
            SampleOutcome::Analyzed(analysis) => Ok(Analysed::Analyzed(analysis)),
            SampleOutcome::Rejected(products) => Ok(Analysed::Rejected {
                reads: products.len(),
            }),
        };
    };
    let documents = traces
        .iter()
        .map(|trace| {
            let sha256 = sha256_file(trace).map_err(|error| Failure::Trace {
                path: trace.clone(),
                message: error.to_string(),
            })?;
            let name = trace.display().to_string();
            let evidence = cache
                .evidence
                .document(&sha256, || recorded::read_evidence_document(trace, config))
                .map_err(failure)?;
            let placed = cache
                .placed
                .document(&format!("{:x}", Sha256::digest(&evidence)), || {
                    recorded::placed_read_document(&name, &evidence, reference, config)
                })
                .map_err(failure)?;
            Ok(recorded::ReadDocuments {
                name,
                evidence,
                placed: Some(placed),
            })
        })
        .collect::<std::result::Result<Vec<_>, Failure>>()?;
    match recorded::analyze_placed_outcome(&documents, reference, config, NonZeroUsize::MIN)
        .map_err(failure)?
    {
        EvidenceOutcome::Analyzed(analysis) => Ok(Analysed::Analyzed(analysis)),
        EvidenceOutcome::Rejected { reads } => Ok(Analysed::Rejected { reads }),
    }
}

/// Runs `call` over every sample's traces with `threads` workers.
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
    let (traces, result) = match distinct(&sample.traces) {
        Ok(()) => {
            let ran = panic::catch_unwind(AssertUnwindSafe(|| call(sample, &sample.traces)))
                .unwrap_or(Err(Failure::Panic));
            (Some(sample.traces.len()), ran)
        }
        Err(failure) => (None, Err(failure)),
    };
    Outcome {
        traces,
        result,
        millis: started.elapsed().as_millis(),
    }
}

/// Checks that no two traces have one content: two inputs of one content
/// are one trace given twice, and the sample fails typed (ADR-0090), rather
/// than counting it twice or dropping a copy silently.
fn distinct(traces: &[PathBuf]) -> std::result::Result<(), Failure> {
    let mut seen = BTreeSet::new();
    for path in traces {
        let digest = sha256_file(path).map_err(|error| Failure::Trace {
            path: path.clone(),
            message: error.to_string(),
        })?;
        if !seen.insert(digest) {
            return Err(Failure::DuplicateTrace { path: path.clone() });
        }
    }
    Ok(())
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
    use dna_consensus::model::consensus::{ConsensusSite, ConsensusSummary};

    #[test]
    fn counts_the_consensus_decisions_an_empop_range_rests_on() {
        let site = |first: usize, reference: &str, call: Option<&str>, state| ConsensusSite {
            first_0based: first - 1,
            last_0based: first + reference.len() - 2,
            reference: reference.into(),
            call: call.map(Into::into),
            state,
            supporting: Vec::new(),
            opposing: Vec::new(),
            uninformative: Vec::new(),
            runs: Vec::new(),
        };
        let consensus = Consensus {
            segments: Vec::new(),
            sites: vec![
                // Undecided, outside the ranges: sound.
                site(10, "A", None, SiteState::Unresolved),
                // Undecided inside a range.
                site(20, "C", None, SiteState::Unresolved),
                // A stretch read position by position: its `N` at 31 lies
                // outside the ranges, its decided bases inside.
                site(30, "CCCT", Some("CNCC"), SiteState::Called),
                // A length change: two `N`s, one position left out.
                site(40, "CCT", Some("CCCNN"), SiteState::Called),
                // A tie the consensus resolved as the reference.
                site(50, "G", Some("G"), SiteState::Contested),
            ],
            summary: ConsensusSummary::default(),
            reference_only: BTreeSet::new(),
            unreliable_only: BTreeMap::new(),
        };
        let export = EmpopExport {
            ranges: vec![[11, 30], [32, 41], [50, 50]],
            mitotype: Vec::new(),
            cut: Vec::new(),
            ignored: Vec::new(),
        };
        let empop = Empop::of(&export, &consensus);
        assert_eq!((empop.undecided, empop.ties), (2, 1));
        assert_eq!(empop.ranges, "11-30 32-41 50-50");
    }

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
            batch_date: None,
            batch_month: None,
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
    fn a_trace_given_twice_fails_the_sample() -> std::result::Result<(), Box<dyn std::error::Error>>
    {
        let directory = tempfile::tempdir()?;
        let paths = ["b", "a", "copy-of-b"].map(|name| directory.path().join(name));
        fs::write(&paths[0], b"trace one")?;
        fs::write(&paths[1], b"trace two")?;
        fs::write(&paths[2], b"trace one")?;
        assert_eq!(distinct(&paths[..2]), Ok(()));
        assert_eq!(
            distinct(&paths),
            Err(Failure::DuplicateTrace {
                path: paths[2].clone()
            })
        );
        let missing = directory.path().join("missing");
        assert!(matches!(distinct(&[missing]), Err(Failure::Trace { .. })));
        Ok(())
    }

    #[test]
    fn outcomes_follow_sample_order_for_any_thread_count()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        let directory = tempfile::tempdir()?;
        let trace = directory.path().join("trace");
        let other = directory.path().join("other");
        fs::write(&trace, b"trace")?;
        fs::write(&other, b"other")?;
        let samples = (0..40)
            .map(|index| sample(&format!("S{index:02}"), vec![trace.clone(), other.clone()]))
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
                    not_recorded: BTreeSet::new(),
                    empop: None,
                    checks: Vec::new(),
                    check_labels: Vec::new(),
                    notes: Vec::new(),
                    components: Vec::new(),
                    tracts: Vec::new(),
                    length_context: None,
                    not_analysed: Vec::new(),
                    intervals: Vec::new(),
                    signals: lod::Signals::default(),
                    loci: Loci::default(),
                    identity: identity(),
                }),
            }
        };
        let call = |sample: &Sample, traces: &[PathBuf]| -> Ran { called(sample, traces) };
        let strip = |outcomes: Vec<Outcome>| {
            outcomes
                .into_iter()
                .map(|outcome| (outcome.traces, outcome.result))
                .collect::<Vec<_>>()
        };
        let one = strip(run_with(&samples, Threads::new(1)?, call)?);
        let many = strip(run_with(&samples, Threads::new(7)?, call)?);
        assert_eq!(one, many);
        assert_eq!(one[0].1, Err(Failure::Dna("fails with 2 traces".into())));
        assert_eq!(one[1].1, Err(Failure::Panic));
        assert_eq!(
            one[2].1.as_ref().map(|called| called.alleles.clone()),
            Ok(vec![Allele::new("02", "A", "G")])
        );
        assert_eq!(one[0].0, Some(2));
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
