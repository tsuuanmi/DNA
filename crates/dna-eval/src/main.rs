//! `dna-eval`: runs DNA over a test panel and compares it with reviewer truth
//! (PROP-0005, ADR-0075).
//!
//! `dna-eval run` reads a panel manifest, runs the recommended `consensus` then
//! `call` workflow for every sample in process, compares the calls with the
//! reviewer truth of the lab's statistics database, and writes the
//! comparisons, weighted estimates, anchors, and a summary into a new output
//! directory. `dna-eval diff` compares two runs and gates the candidate;
//! `dna-eval panel` derives a new panel version against the database;
//! `dna-eval regress` checks DNA's output on the frozen regression set
//! against its version ledger.
//! Nothing in DNA depends on this crate.

mod bench;
mod cache;
mod cause;
mod classes;
mod compare;
mod components;
#[cfg(test)]
mod contract;
mod conventions;
mod error;
mod estimate;
mod geometry;
mod haplotype;
mod lod;
mod manifest;
mod panel;
mod recheck;
mod regress;
mod report;
mod rundiff;
mod rundir;
mod runner;
mod same;
mod simbench;
mod subset;
mod tracts;
mod truth;
mod unrecorded;
mod workload;
mod zones;

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::fs::{self, File};
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::thread;
use std::time::Instant;

use clap::{Args, Parser, Subcommand};
use dna_kernel::model::reference::ReferenceTopology;
use sha2::{Digest, Sha256};

use crate::cause::forms::Forms;
use crate::compare::Measure;
use crate::error::{Error, Result};
use crate::geometry::Geometry;
use crate::manifest::{Manifest, Month};
use crate::recheck::Recheck;
use crate::report::{Report, Row, Status};
use crate::rundir::{
    Constants, Counts, ManifestHeader, RecheckHeader, RunContextHeader, RunHeader, Seconds,
    SignalGrouping, Source, TruthHeader, WorkloadHeader, ZoneHeader, prepare_output, write,
};
use crate::runner::{Outcome, Threads};
use crate::truth::TruthCheck;

/// The panels the harness reads: the development panel the rules are tuned
/// on, the readiness panels that test them, and the regression set, which is
/// development data (`dna-eval regress`). A hold-out panel stays reserved
/// for an approval path not yet designed.
const RUNNABLE_PANELS: [&str; 3] = ["dev", "readiness", "regress"];

/// Largest accepted panel manifest.
const MAX_MANIFEST_BYTES: u64 = 256 * 1024 * 1024;
/// Largest accepted forced-window list.
const MAX_FORCED_BYTES: u64 = 1024 * 1024;
/// Largest accepted re-check list.
const MAX_RECHECK_BYTES: u64 = 16 * 1024 * 1024;
/// Largest accepted bench case list.
const MAX_CASES_BYTES: u64 = 4 * 1024 * 1024;

/// Runs DNA over a test panel and compares it with reviewer truth.
#[derive(Debug, Parser)]
#[command(name = "dna-eval", version)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

// Parsed once per process; boxing the run's arguments would buy nothing.
#[allow(clippy::large_enum_variant)]
#[derive(Debug, Subcommand)]
enum Command {
    /// Runs DNA over a panel and compares it with the truth.
    Run(RunArgs),
    /// Derives a new panel version against the truth database.
    Panel(PanelArgs),
    /// Writes a quick stratified subset of a panel, for a first-pass `diff`
    /// before the full panel: a pre-filter, never the gate.
    Subset(SubsetArgs),
    /// Manages a read evidence cache (`run --evidence-cache`).
    Cache(CacheArgs),
    /// Compares a run with a base run of the same panel and gates it: exits
    /// 2 with `GATE FAIL` when, in the current period (all samples without
    /// `--current-since`), more samples, passed samples, or passed zones
    /// carry a non-het unreviewed block; het blocks and the earlier period
    /// are reported, not gated.
    Diff(DiffArgs),
    /// Says whether two runs' outputs are byte-identical (every file but
    /// `run.json`, by SHA-256) and which files differ: exits 2 with
    /// `DIFFERENT` when one does.
    Same(SameArgs),
    /// Runs DNA on the adjudicated artefact cases of a case list and checks
    /// each expected outcome: exits 2 with `BENCH FAIL` when a guarded case
    /// regresses or its sample cannot be analysed.
    Bench(BenchArgs),
    /// Runs DNA's read evidence stage on the frozen simulated tract-decode
    /// cases (PROP-0015) and reports the decoder's measures per tract and
    /// noise level, against a baseline when one is given.
    SimBench(SimBenchArgs),
    /// Checks DNA's output on the frozen regression set against the head of
    /// its version ledger (`check`, exits 2 with `REGRESS FAIL` when a change
    /// is undeclared or forbidden), or records a checked run of `main` in
    /// the ledger (`record`).
    Regress(regress::Args),
}

/// How DNA runs over a panel's samples: the options `run`, `bench`, and
/// `regress check` share.
#[derive(Debug, Args)]
pub(crate) struct DnaOptions {
    /// DNA configuration, as given to `dna` through `DNA_CONFIG`.
    #[arg(long)]
    config: PathBuf,
    /// Reference FASTA, as given to `dna --reference`.
    #[arg(long)]
    reference: PathBuf,
    /// Worker threads that run DNA (1 to 256).
    #[arg(long)]
    threads: usize,
    /// A read evidence cache: every read is analysed from its
    /// `dna.read_evidence/v1` document, read from the directory or made and
    /// stored there, so a run whose read stage is unchanged skips it. It
    /// needs a harness built from committed read-stage sources; `run.json`
    /// records its key and hits. The directory also caches the index of the
    /// truth's written forms, per truth database, reference, and repeat
    /// windows.
    #[arg(long)]
    evidence_cache: Option<PathBuf>,
    /// A directory of run contexts (PROP-0014), one
    /// `<batch>.run_context.json` per batch, as `dna run-context` writes them
    /// over the batches' full folders: a sample of a batch with one has its
    /// consensus calibrated by it. `run_context.tsv` reports per batch month
    /// how many samples a context calibrated.
    #[arg(long)]
    run_context: Option<PathBuf>,
}

/// How the review workload is counted: the options `run` and `regress
/// check` share.
#[derive(Debug, Args)]
pub(crate) struct ReviewOptions {
    /// Window of the review workload, in positions, aligned to each HV
    /// region's first position.
    #[arg(long, default_value_t = conventions::REVIEW_WINDOW)]
    review_window: i64,
    /// Review zones: positions on each side of a repeat window that follow
    /// it.
    #[arg(long, default_value_t = conventions::ZONE_AFTER)]
    zone_after: i64,
    /// Review zones: positions at each end of a region that are its
    /// amplicon's start or end.
    #[arg(long, default_value_t = conventions::ZONE_EDGE)]
    zone_edge: i64,
    /// Review zones: positions past a repeat window that an item overlapping
    /// the window still counts for the window.
    #[arg(long, default_value_t = conventions::ZONE_OVERRUN)]
    zone_overrun: i64,
}

impl ReviewOptions {
    /// The zone widths, once the window and every width are checked.
    fn widths(&self) -> Result<zones::Widths> {
        if self.review_window < 1 {
            return Err(Error::ReviewWindow {
                window: self.review_window,
            });
        }
        let widths = zones::Widths {
            after: self.zone_after,
            edge: self.zone_edge,
            overrun: self.zone_overrun,
        };
        for (name, value) in [
            ("--zone-after", widths.after),
            ("--zone-edge", widths.edge),
            ("--zone-overrun", widths.overrun),
        ] {
            if value < 0 {
                return Err(Error::ZoneWidth { name, value });
            }
        }
        Ok(widths)
    }

    /// The widths as `run.json` records them.
    const fn header(&self) -> ZoneHeader {
        ZoneHeader {
            after: self.zone_after,
            edge: self.zone_edge,
            overrun: self.zone_overrun,
        }
    }
}

/// Arguments of `dna-eval sim-bench`.
#[derive(Debug, Args)]
struct SimBenchArgs {
    /// DNA configuration, as given to `dna` through `DNA_CONFIG`.
    #[arg(long)]
    config: PathBuf,
    /// Reference FASTA the templates are read from.
    #[arg(long)]
    reference: PathBuf,
    /// Worker threads that run DNA (1 to 256).
    #[arg(long)]
    threads: usize,
    /// A new directory for `bench.md`, `bench.tsv`, and `cases.tsv`.
    #[arg(long)]
    output: PathBuf,
    /// An earlier `bench.tsv` to report the change from, in place of the
    /// committed baseline of the set's size.
    #[arg(long)]
    baseline: Option<PathBuf>,
    /// Repetitions per cell: 60, the gating set (the default), or 6, the
    /// smoke set, each with its committed baseline; any other size is a
    /// validation run without one.
    #[arg(long, default_value_t = simbench::GATE_REPS)]
    reps: usize,
    /// Keep the rendered traces in `<output>/traces`, for research checks of
    /// the renderer's realism.
    #[arg(long)]
    keep_traces: bool,
}

/// Runs `dna-eval sim-bench` and prints its report.
fn sim_bench(args: &SimBenchArgs) -> Result<String> {
    let threads = Threads::new(args.threads)?;
    let reference = dna_kernel::reference::load(&args.reference, ReferenceTopology::Linear)
        .map_err(|error| Error::Reference {
            path: args.reference.clone(),
            message: error.to_string(),
        })?;
    let table =
        dna_trace::simulate::Table::pooled().map_err(|error| Error::SimBench(error.to_string()))?;
    let baseline = match args.baseline.as_deref() {
        Some(path) => {
            let text = String::from_utf8(read(path, "baseline", MAX_CASES_BYTES)?)
                .map_err(|error| Error::SimBench(error.to_string()))?;
            Some(
                simbench::parse(&text)
                    .map_err(|message| Error::SimBench(format!("{}: {message}", path.display())))?,
            )
        }
        None => simbench::baseline(args.reps)
            .map(simbench::parse)
            .transpose()
            .map_err(|message| Error::SimBench(format!("the committed baseline: {message}")))?,
    };
    prepare_output(&args.output)?;
    let scratch = args.output.join("traces");
    fs::create_dir_all(&scratch).map_err(|source| Error::Write {
        path: scratch.clone(),
        source,
    })?;
    let cases = simbench::cases(args.reps);
    let outcomes = simbench::run(
        &cases,
        &table,
        reference.sequence.as_bytes(),
        &args.config,
        (&scratch, args.keep_traces),
        threads,
    );
    if !args.keep_traces {
        let _ = fs::remove_dir_all(&scratch);
    }
    let rows = simbench::measures(&outcomes);
    let source = source();
    let revision = format!(
        "{}{}",
        source.revision.as_deref().unwrap_or("unknown revision"),
        match source.dirty {
            Some(true) => " (dirty)",
            Some(false) => "",
            None => " (dirty state unknown)",
        }
    );
    let report = simbench::report(&revision, args.reps, &rows, &outcomes, baseline.as_ref());
    write(&args.output, "bench.md", &report)?;
    write(&args.output, "bench.tsv", &simbench::table(&rows))?;
    write(&args.output, "cases.tsv", &simbench::cases_table(&outcomes))?;
    Ok(report)
}

/// Arguments of `dna-eval bench`.
#[derive(Debug, Args)]
struct BenchArgs {
    /// The case list (TSV: `case`, `sample_id`, `panel`, `reads`, `class`,
    /// `first`, `last`, `expect`, `status`, `source`).
    #[arg(long)]
    cases: PathBuf,
    /// A panel manifest the cases name as `<panel>-<version>`; repeatable.
    #[arg(long, required = true)]
    manifest: Vec<PathBuf>,
    #[command(flatten)]
    dna: DnaOptions,
    /// A new directory for `bench.md` and `bench.tsv`.
    #[arg(long)]
    output: Option<PathBuf>,
}

/// Runs `dna-eval bench`: prints the report, writes its files when asked,
/// and exits 2 when a guarded case fails.
fn bench(args: &BenchArgs) -> ExitCode {
    match run_bench(args) {
        Ok(report) => {
            print!("{}", report.text);
            if report.failed {
                ExitCode::from(2)
            } else {
                ExitCode::SUCCESS
            }
        }
        Err(error) => {
            eprintln!("dna-eval: {error}");
            ExitCode::FAILURE
        }
    }
}

/// Runs DNA once on each sample the cases name and judges every case.
fn run_bench(args: &BenchArgs) -> Result<bench::Report> {
    let threads = Threads::new(args.dna.threads)?;
    let invalid = |message: String| Error::Cases {
        path: args.cases.clone(),
        message,
    };
    let text = String::from_utf8(read(&args.cases, "case list", MAX_CASES_BYTES)?)
        .map_err(|error| invalid(error.to_string()))?;
    let cases = bench::parse(&text).map_err(invalid)?;
    let mut panels = BTreeMap::new();
    for path in &args.manifest {
        let manifest = manifest::parse(&read(path, "manifest", MAX_MANIFEST_BYTES)?)?;
        runnable(&manifest)?;
        let (Some(panel), Some(version)) =
            (manifest.header.get("panel"), manifest.header.get("version"))
        else {
            return Err(invalid(format!(
                "manifest {} declares no panel and version",
                path.display()
            )));
        };
        let label = format!("{panel}-{version}");
        if panels.insert(label.clone(), manifest).is_some() {
            return Err(invalid(format!("two manifests are panel {label}")));
        }
    }
    // Each sample runs once, however many cases read it.
    let mut samples = Vec::new();
    let mut index = BTreeMap::new();
    for case in &cases {
        let key = (case.panel.as_str(), case.sample.as_str());
        if index.contains_key(&key) {
            continue;
        }
        let manifest = panels.get(&case.panel).ok_or_else(|| {
            invalid(format!(
                "case {} names panel {}, which no --manifest is",
                case.name, case.panel
            ))
        })?;
        let sample = manifest
            .samples
            .iter()
            .find(|sample| sample.id == case.sample)
            .ok_or_else(|| {
                invalid(format!(
                    "case {}: panel {} has no sample {}",
                    case.name, case.panel, case.sample
                ))
            })?;
        index.insert(key, samples.len());
        samples.push(sample.clone());
    }
    let reference = load_reference(&args.dna.reference)?;
    let cache = args.dna.caches(&reference.sequence_sha256)?;
    let (contexts, contexts_sha256) = args.dna.contexts(&samples)?;
    let outcomes = runner::run(
        &samples,
        (&args.dna.reference, &reference.sequence_sha256),
        &args.dna.config,
        threads,
        None,
        (cache.as_ref(), &contexts),
        &mut Vec::new(),
    )?;
    let judged = cases
        .iter()
        .map(|case| {
            let result = &outcomes[index[&(case.panel.as_str(), case.sample.as_str())]].result;
            let (verdict, observed) = bench::judge(case, result);
            (
                verdict,
                observed,
                result.as_ref().err().map(ToString::to_string),
            )
        })
        .collect::<Vec<_>>();
    let source = source();
    let mut revision = format!(
        "{}{}",
        source.revision.as_deref().unwrap_or("unknown revision"),
        match source.dirty {
            Some(true) => " (dirty)",
            Some(false) => "",
            None => " (dirty state unknown)",
        }
    );
    if let Some(directory) = &args.dna.run_context {
        let _ = write!(
            revision,
            "; run contexts {} ({} batches, sha256 {contexts_sha256})",
            directory.display(),
            contexts.len()
        );
    }
    let report = bench::report(&revision, &cases, &judged);
    if let Some(output) = &args.output {
        prepare_output(output)?;
        write(output, "bench.md", &report.text)?;
        write(output, "bench.tsv", &report.table)?;
    }
    Ok(report)
}

/// Arguments of `dna-eval diff`.
#[derive(Debug, Args)]
struct DiffArgs {
    /// The base run directory.
    #[arg(long)]
    base: PathBuf,
    /// The run directory compared with it.
    #[arg(long)]
    run: PathBuf,
    /// A new directory for `diff.md`, `samples.tsv`, and `statuses.tsv`.
    #[arg(long)]
    output: Option<PathBuf>,
    /// Compare runs that measure different things (another manifest, truth,
    /// re-check list, current month, window, or zones, or a dirty build),
    /// with a warning in the report.
    #[arg(long)]
    force: bool,
}

/// Arguments of `dna-eval same`.
#[derive(Debug, Args)]
struct SameArgs {
    /// A run directory.
    a: PathBuf,
    /// The run directory compared with it.
    b: PathBuf,
}

/// Runs `dna-eval same`: prints the comparison, and exits 2 when an output
/// differs.
fn same(args: &SameArgs) -> ExitCode {
    match same::compare(&args.a, &args.b) {
        Ok(sameness) => {
            print!("{}", sameness.report);
            if sameness.same {
                ExitCode::SUCCESS
            } else {
                ExitCode::from(2)
            }
        }
        Err(error) => {
            eprintln!("dna-eval: {error}");
            ExitCode::FAILURE
        }
    }
}

/// Arguments of `dna-eval cache`.
#[derive(Debug, Args)]
struct CacheArgs {
    #[command(subcommand)]
    action: CacheAction,
}

#[derive(Debug, Subcommand)]
enum CacheAction {
    /// Removes every generation of the cache but the most recently used.
    Prune {
        /// The cache directory.
        #[arg(long)]
        dir: PathBuf,
        /// Generations to keep.
        #[arg(long)]
        keep: usize,
    },
}

/// Runs `dna-eval diff`: prints the comparison, writes its files when asked,
/// and exits 2 when the gate fails.
fn diff(args: &DiffArgs) -> ExitCode {
    let compared = rundiff::compare(&args.base, &args.run, args.force)
        .map_err(Error::from)
        .and_then(|comparison| {
            if let Some(output) = &args.output {
                prepare_output(output)?;
                for (name, contents) in [
                    ("diff.md", &comparison.report),
                    ("samples.tsv", &comparison.samples),
                    ("statuses.tsv", &comparison.statuses),
                ] {
                    write(output, name, contents)?;
                }
            }
            Ok(comparison)
        });
    match compared {
        Ok(comparison) => {
            print!("{}", comparison.report);
            if comparison.gate_failed {
                ExitCode::from(2)
            } else {
                ExitCode::SUCCESS
            }
        }
        Err(error) => {
            eprintln!("dna-eval: {error}");
            ExitCode::FAILURE
        }
    }
}

/// Arguments of `dna-eval panel`.
#[derive(Debug, Args)]
struct PanelArgs {
    /// The parent panel manifest.
    #[arg(long)]
    manifest: PathBuf,
    /// The lab's statistics database (`statistics.sqlite`), opened read-only.
    #[arg(long)]
    truth: PathBuf,
    /// The new version label, for example `v2`.
    #[arg(long)]
    version: String,
    /// The creation date to record, for example `2026-10-15`.
    #[arg(long)]
    created: String,
    /// The new manifest file; it must not exist.
    #[arg(long)]
    output: PathBuf,
}

/// Arguments of `dna-eval subset`.
#[derive(Debug, Args)]
struct SubsetArgs {
    /// The parent panel manifest.
    #[arg(long)]
    manifest: PathBuf,
    /// Rows to share among the strata, in proportion to their staying rows.
    #[arg(long)]
    rows: usize,
    /// Fewest rows a stratum keeps (all of them when it has fewer).
    #[arg(long)]
    minimum: usize,
    /// Salt of the deterministic order within a stratum.
    #[arg(long)]
    salt: String,
    /// The new manifest file; it must not exist.
    #[arg(long)]
    output: PathBuf,
}

/// Arguments of `dna-eval run`.
#[derive(Debug, Args)]
struct RunArgs {
    /// Panel manifest (TSV with a `#` comment header).
    #[arg(long)]
    manifest: PathBuf,
    /// The lab's statistics database (`statistics.sqlite`), opened read-only.
    #[arg(long)]
    truth: PathBuf,
    /// New or empty output directory.
    #[arg(long)]
    output: PathBuf,
    #[command(flatten)]
    dna: DnaOptions,
    #[command(flatten)]
    review: ReviewOptions,
    /// Windows always sent to review, one `start-end` (1-based, inclusive)
    /// per line; the workload is reported again with them reviewed. A list
    /// learned on this run's panel overstates its effect.
    #[arg(long)]
    forced_windows: Option<PathBuf>,
    /// Difference blocks whose truth is under re-check (TSV: `sample_id`,
    /// `first`, `last`, `reason`, `source`): the pass measures count a block
    /// that overlaps a listed span of its sample as re-check, not as a DNA
    /// error.
    #[arg(long)]
    recheck: Option<PathBuf>,
    /// First batch month (`YYYY-MM`) of the current data: the summary reports
    /// the samples whose manifest `batch_date` is in or after it first, and
    /// `diff` gates on them alone.
    #[arg(long)]
    current_since: Option<Month>,
    /// Judge each review against a run signal (PROP-0009) built in process
    /// from the panel: `batch` takes a sample's batch when it has at least
    /// `--run-signal-batch-reads` reads, else its batch month; `month` takes
    /// the month alone. The run takes two passes.
    #[arg(long, value_enum)]
    run_signal: Option<SignalGrouping>,
    /// Fewest reads of a batch for `--run-signal batch` to use it.
    #[arg(long, default_value_t = conventions::RUN_SIGNAL_BATCH_READS)]
    run_signal_batch_reads: usize,
}

fn main() -> ExitCode {
    let result = match Cli::parse().command {
        Command::Diff(args) => return diff(&args),
        Command::Same(args) => return same(&args),
        Command::Bench(args) => return bench(&args),
        Command::SimBench(args) => sim_bench(&args),
        Command::Regress(args) => return regress::main(&args),
        Command::Run(args) => run(&args),
        Command::Panel(args) => derive_panel(&args),
        Command::Subset(args) => write_subset(&args),
        Command::Cache(CacheArgs {
            action: CacheAction::Prune { dir, keep },
        }) => cache::prune(&dir, keep)
            .map_err(Error::from)
            .map(|removed| {
                let bytes = removed.iter().map(|(_, bytes)| bytes).sum::<u64>();
                format!(
                    "dna-eval: removed {} cache generation(s), {bytes} bytes",
                    removed.len()
                )
            }),
    };
    match result {
        Ok(summary) => {
            println!("{summary}");
            ExitCode::SUCCESS
        }
        Err(error) => {
            eprintln!("dna-eval: {error}");
            ExitCode::FAILURE
        }
    }
}

/// Whether the harness may read the panel: a development or readiness panel,
/// as its header names it. A manifest that names no panel could be a
/// hold-out's, so it is refused too.
pub(crate) fn runnable(manifest: &Manifest) -> Result<()> {
    match manifest.header.get("panel") {
        Some(panel) if RUNNABLE_PANELS.contains(&panel.as_str()) => Ok(()),
        Some(panel) => Err(Error::HoldOut {
            panel: panel.clone(),
        }),
        None => Err(Error::UnnamedPanel),
    }
}

/// A panel's inputs, loaded and checked before DNA runs: what `run` and
/// `regress check` share.
pub(crate) struct Loaded {
    pub(crate) manifest: Manifest,
    pub(crate) truth: Truth,
    pub(crate) reference: dna_kernel::model::reference::Reference,
    pub(crate) profile: dna_kernel::profile::Profile,
    pub(crate) geometry: Geometry,
    pub(crate) forms: Forms,
    /// The forms index's cache use, with an evidence cache directory.
    pub(crate) forms_cache: Option<cache::FormsUsage>,
    /// The re-check list and its identity, when given.
    pub(crate) recheck: Option<(Recheck, RecheckHeader)>,
    pub(crate) widths: zones::Widths,
}

impl Loaded {
    /// Reads the manifest at `manifest`, which `accept` must allow, its
    /// whole truth from `truth`, the reference and the target profile, the
    /// written forms (from the evidence cache directory when stored there,
    /// else indexed over the whole database), and the re-check list when
    /// given; checks the review options.
    pub(crate) fn read(
        (manifest, accept): (&Path, fn(&Manifest) -> Result<()>),
        (truth, recheck): (&Path, Option<&Path>),
        dna: &DnaOptions,
        review: &ReviewOptions,
    ) -> Result<Self> {
        let manifest = manifest::parse(&read(manifest, "manifest", MAX_MANIFEST_BYTES)?)?;
        accept(&manifest)?;
        // The whole truth of the panel loads before DNA runs, so a database
        // that cannot be opened, has another schema, or fails a query stops
        // the run before any sample is called.
        let wanted = manifest
            .samples
            .iter()
            .map(|sample| sample.id.as_str())
            .collect::<BTreeSet<_>>();
        let truth_records = load_truth(truth, &wanted)?;
        // The haplotype measures apply alleles to this sequence; a sample DNA
        // analyses on another sequence fails.
        let reference = load_reference(&dna.reference)?;
        // The geometry the run compares and counts in, from the target
        // profile.
        let profile = target_profile(&dna.config)?;
        let geometry = Geometry::of(&profile);
        // How the reviewers write each window sequence, over the whole truth,
        // for the cause of a representation discrepancy.
        let (forms, forms_cache) = forms(
            (truth, &truth_records.sha256),
            dna.evidence_cache.as_deref(),
            (&reference.sequence, &reference.sequence_sha256),
            &geometry,
        )?;
        let widths = review.widths()?;
        let recheck = recheck
            .map(|path| -> Result<_> {
                let failure = |message: String| Error::Recheck {
                    path: path.to_path_buf(),
                    message,
                };
                let bytes = read(path, "re-check list", MAX_RECHECK_BYTES)?;
                let sha256 = format!("{:x}", Sha256::digest(&bytes));
                let text = String::from_utf8(bytes).map_err(|error| failure(error.to_string()))?;
                let list = Recheck::parse(&text).map_err(failure)?;
                let header = RecheckHeader {
                    path: absolute(path),
                    sha256,
                    rows: list.rows(),
                };
                Ok((list, header))
            })
            .transpose()?;
        Ok(Self {
            manifest,
            truth: truth_records,
            reference,
            profile,
            geometry,
            forms,
            forms_cache,
            recheck,
            widths,
        })
    }

    /// The profile's length tracts, as the review workload notes them.
    pub(crate) fn tracts(&self) -> Vec<workload::Tract> {
        self.profile
            .length_tracts
            .iter()
            .map(|tract| workload::Tract {
                name: tract.name.clone(),
                first: i64::try_from(tract.window[0]).unwrap_or(i64::MAX),
                last: i64::try_from(tract.window[1]).unwrap_or(i64::MAX),
                base: char::from(tract.base),
            })
            .collect()
    }

    /// The review zones of the profile's regions, placed with the widths.
    pub(crate) fn zones(&self) -> zones::Zones {
        zones::Zones::new(
            self.geometry.regions(),
            self.geometry.repeats(),
            self.widths,
        )
    }
}

/// DNA's outcomes over a loaded panel, in manifest order, with the caches
/// and run contexts it read.
pub(crate) struct Ran {
    pub(crate) outcomes: Vec<Outcome>,
    pub(crate) cache: Option<cache::Caches>,
    pub(crate) contexts: runner::Contexts,
    pub(crate) contexts_sha256: String,
    /// The run signal groups, with a grouping.
    pub(crate) signal_groups: Vec<String>,
    /// Wall time of DNA, in seconds.
    pub(crate) seconds: f64,
}

impl Ran {
    /// Runs DNA over every sample of `loaded` as `dna` says, with each
    /// review judged against its group's run signal under `grouping`, and
    /// checks that DNA analysed with the configuration's target profile.
    pub(crate) fn run(
        loaded: &Loaded,
        dna: &DnaOptions,
        grouping: Option<runner::Grouping>,
    ) -> Result<Self> {
        let threads = Threads::new(dna.threads)?;
        let started = Instant::now();
        let cache = dna.caches(&loaded.reference.sequence_sha256)?;
        let (contexts, contexts_sha256) = dna.contexts(&loaded.manifest.samples)?;
        let mut signal_groups = Vec::new();
        let outcomes = runner::run(
            &loaded.manifest.samples,
            (&dna.reference, &loaded.reference.sequence_sha256),
            &dna.config,
            threads,
            grouping,
            (cache.as_ref(), &contexts),
            &mut signal_groups,
        )?;
        let seconds = started.elapsed().as_secs_f64();
        let analysed = &loaded.profile.identity().sha256;
        if let Some(other) = outcomes
            .iter()
            .filter_map(|outcome| outcome.result.as_ref().ok())
            .map(|called| &called.identity.profile_sha256)
            .find(|sha256| *sha256 != analysed)
        {
            return Err(Error::Profile {
                path: dna.config.clone(),
                message: format!("DNA analysed with profile {other}, not {analysed}"),
            });
        }
        Ok(Self {
            outcomes,
            cache,
            contexts,
            contexts_sha256,
            signal_groups,
            seconds,
        })
    }
}

impl DnaOptions {
    /// The evidence and placement caches, when a cache directory is given.
    fn caches(&self, reference_sha256: &str) -> Result<Option<cache::Caches>> {
        Ok(self
            .evidence_cache
            .as_deref()
            .map(|directory| {
                Ok::<_, cache::CacheError>(cache::Caches {
                    evidence: cache::EvidenceCache::open(directory, &self.config)?,
                    placed: cache::EvidenceCache::open_placed(
                        directory,
                        &self.config,
                        reference_sha256,
                    )?,
                })
            })
            .transpose()?)
    }

    /// The run contexts of the samples' batches and their digest, when a
    /// run-context directory is given; none otherwise.
    fn contexts(&self, samples: &[manifest::Sample]) -> Result<(runner::Contexts, String)> {
        Ok(self
            .run_context
            .as_deref()
            .map(|directory| runner::contexts(directory, samples))
            .transpose()?
            .unwrap_or_default())
    }

    /// The run contexts as `run.json` records them.
    fn context_header(&self, ran: &Ran) -> Option<RunContextHeader> {
        self.run_context
            .as_deref()
            .map(|directory| RunContextHeader {
                directory: absolute(directory),
                batches: ran.contexts.len(),
                sha256: ran.contexts_sha256.clone(),
            })
    }
}

/// Loads the reference FASTA at `path`.
fn load_reference(path: &Path) -> Result<dna_kernel::model::reference::Reference> {
    dna_kernel::reference::load(path, ReferenceTopology::Linear).map_err(|error| Error::Reference {
        path: path.to_path_buf(),
        message: error.to_string(),
    })
}

fn run(cli: &RunArgs) -> Result<String> {
    let started = Instant::now();
    let mut loaded = Loaded::read(
        (&cli.manifest, runnable),
        (&cli.truth, cli.recheck.as_deref()),
        &cli.dna,
        &cli.review,
    )?;
    let forced = cli
        .forced_windows
        .as_ref()
        .map(|path| -> Result<_> {
            let text = read(path, "forced windows", MAX_FORCED_BYTES)?;
            let text = String::from_utf8(text).map_err(|error| Error::ForcedWindows {
                path: path.clone(),
                message: error.to_string(),
            })?;
            workload::parse_forced(&text).map_err(|message| Error::ForcedWindows {
                path: path.clone(),
                message,
            })
        })
        .transpose()?;
    prepare_output(&cli.output)?;
    let ran = Ran::run(
        &loaded,
        &cli.dna,
        cli.run_signal.map(|grouping| match grouping {
            SignalGrouping::Month => runner::Grouping::Month,
            SignalGrouping::Batch => runner::Grouping::Batch {
                batch_reads: cli.run_signal_batch_reads,
            },
        }),
    )?;
    let (recheck, recheck_header) = loaded.recheck.take().unzip();
    let zones = loaded.zones();
    let tracts = loaded.tracts();
    let Loaded {
        manifest,
        truth,
        reference,
        profile,
        geometry,
        forms,
        ..
    } = &loaded;

    let rows = report::rows(
        &manifest.samples,
        &ran.outcomes,
        &truth.database.samples,
        &reference.sequence,
        (geometry, &profile.conventions),
    );
    let check = TruthCheck::new(manifest, &truth.database);
    let report = Report::new(
        &rows,
        (geometry, &check),
        (forms, reference.sequence.as_str()),
    )
    .with_review(
        cli.review.review_window,
        forced,
        Some(zones),
        recheck,
        tracts,
    )
    .with_current_since(cli.current_since);

    let mut header = header(
        cli,
        manifest,
        truth,
        &report,
        &ran.outcomes,
        Times {
            started,
            dna: ran.seconds,
        },
    );
    header.recheck = recheck_header;
    header.evidence_cache = ran.cache.as_ref().map(|cache| cache.evidence.usage());
    header.forms_cache.clone_from(&loaded.forms_cache);
    header.placement_cache = ran.cache.as_ref().map(|cache| cache.placed.usage());
    header.run_context = cli.dna.context_header(&ran);
    for (name, contents) in report.files(&context(&header)) {
        write(&cli.output, name, &contents)?;
    }
    let (samples, outcomes) = (&manifest.samples, &ran.outcomes);
    if cli.dna.run_context.is_some() {
        write(
            &cli.output,
            rundir::RUN_CONTEXT,
            &rundir::context_coverage(samples, outcomes, &ran.contexts),
        )?;
    }
    if cli.run_signal.is_some() {
        write(
            &cli.output,
            rundir::SIGNALS,
            &rundir::signals(&ran.signal_groups),
        )?;
    }
    write(
        &cli.output,
        rundir::CHECKS,
        &rundir::checks(samples, outcomes),
    )?;
    write(
        &cli.output,
        rundir::NOTES,
        &rundir::notes(samples, outcomes),
    )?;
    write(
        &cli.output,
        rundir::COMPONENTS,
        &rundir::components(
            samples,
            outcomes,
            &truth.database.samples,
            reference.sequence.as_bytes(),
        ),
    )?;
    write(
        &cli.output,
        rundir::TRACTS,
        &rundir::tracts(samples, outcomes),
    )?;
    write(&cli.output, rundir::FORMS, &forms.table(&truth.sha256))?;
    // `run.json` comes last, with the digest of every file before it.
    header.outputs = rundir::digests(&cli.output)?;
    write(&cli.output, rundir::RUN_JSON, &header.text()?)?;

    let count = |measure: Measure, status: Status| count(&rows, measure, status);
    Ok(format!(
        "dna-eval: truth {} ({} samples changed); {} samples; primary {} concordant, {} discordant, {} not applicable, {} \
         failed; strict {} concordant; IUPAC-compatible {} concordant; haplotype {} concordant; \
         haplotype IUPAC-compatible {} concordant; {:.1} s; outputs in {}",
        check.state().label(),
        check.changed.len(),
        rows.len(),
        count(Measure::Primary, Status::Concordant),
        count(Measure::Primary, Status::Discordant),
        count(Measure::Primary, Status::NotApplicable),
        count(Measure::Primary, Status::Failed),
        count(Measure::Strict, Status::Concordant),
        count(Measure::IupacCompatible, Status::Concordant),
        count(Measure::Haplotype, Status::Concordant),
        count(Measure::HaplotypeIupacCompatible, Status::Concordant),
        started.elapsed().as_secs_f64(),
        cli.output.display()
    ))
}

/// The index of how the reviewers write each window sequence, over the
/// whole truth database at `truth` (whose file SHA-256 is `truth_sha256`),
/// on the reference `sequence` (whose SHA-256 is given), over the repeat
/// windows of `geometry`. With an evidence cache directory it is read from
/// the directory when stored under the same key, and stored when built;
/// without a code identity for the index, it is built and not stored.
fn forms(
    (truth, truth_sha256): (&Path, &str),
    cache_root: Option<&Path>,
    (sequence, sequence_sha256): (&str, &str),
    geometry: &Geometry,
) -> Result<(Forms, Option<cache::FormsUsage>)> {
    let build = || -> Result<Forms> {
        let corpus = truth::sample_ids(truth)?;
        let corpus = truth::load(truth, &corpus.iter().map(String::as_str).collect())?;
        Ok(Forms::new(corpus.samples.values(), sequence, geometry))
    };
    let key = cache_root.and_then(|root| {
        let windows = geometry
            .window_bounds()
            .map(|(window, first, last)| (window.label().to_owned(), first, last))
            .collect::<Vec<_>>();
        cache::FormsKey::new(
            truth_sha256,
            sequence_sha256,
            windows
                .iter()
                .map(|(name, first, last)| (name.as_str(), *first, *last)),
        )
        .map(|key| (root, key))
    });
    let Some((root, key)) = key else {
        return Ok((build()?, None));
    };
    let cache = cache::FormsCache::open(root, key)?;
    if let Some(forms) = cache
        .read()
        .and_then(|text| Forms::parse(&text, truth_sha256, geometry).ok())
    {
        return Ok((forms, Some(cache.usage(true))));
    }
    let forms = build()?;
    cache.store(&forms.table(truth_sha256));
    Ok((forms, Some(cache.usage(false))))
}

/// The target profile `config` names, resolved against the configuration
/// file's directory when relative, as DNA resolves it (SRS-CFG-002).
pub(crate) fn target_profile(config: &Path) -> Result<dna_kernel::profile::Profile> {
    let failure = |message: String| Error::Profile {
        path: config.to_path_buf(),
        message,
    };
    let text = fs::read_to_string(config).map_err(|error| failure(error.to_string()))?;
    let value = text
        .parse::<toml::Table>()
        .map_err(|error| failure(error.to_string()))?;
    let name = value
        .get("profile")
        .and_then(toml::Value::as_str)
        .ok_or_else(|| failure("the configuration names no profile".to_owned()))?;
    let path = config.parent().unwrap_or(Path::new(".")).join(name);
    dna_kernel::profile::Profile::load(&path).map_err(|error| Error::Profile {
        path,
        message: error.to_string(),
    })
}

/// `dna-eval panel`: writes the new manifest and returns its checksum line.
fn write_subset(args: &SubsetArgs) -> Result<String> {
    if args.salt.is_empty() || args.salt.contains(char::is_whitespace) {
        return Err(Error::Argument {
            name: "salt",
            value: args.salt.clone(),
        });
    }
    for (name, value) in [("rows", args.rows), ("minimum", args.minimum)] {
        if value == 0 {
            return Err(Error::Argument {
                name,
                value: value.to_string(),
            });
        }
    }
    let parent = manifest::parse(&read(&args.manifest, "manifest", MAX_MANIFEST_BYTES)?)?;
    let quick = subset::subset(&parent, args.rows, args.minimum, &args.salt);
    File::options()
        .write(true)
        .create_new(true)
        .open(&args.output)
        .and_then(|mut file| file.write_all(quick.text.as_bytes()))
        .map_err(|source| Error::Write {
            path: args.output.clone(),
            source,
        })?;
    let strata = quick
        .strata
        .iter()
        .map(|(stratum, (kept, of))| format!("{stratum} {kept}/{of}"))
        .collect::<Vec<_>>()
        .join(", ");
    Ok(format!(
        "dna-eval: {} rows in {}; per stratum: {strata}",
        quick.rows,
        args.output.display()
    ))
}

fn derive_panel(args: &PanelArgs) -> Result<String> {
    for (name, value) in [("version", &args.version), ("created", &args.created)] {
        if value.is_empty()
            || value.contains(|symbol: char| symbol.is_whitespace() || symbol == '=')
        {
            return Err(Error::Argument {
                name,
                value: value.clone(),
            });
        }
    }
    let parent = manifest::parse(&read(&args.manifest, "manifest", MAX_MANIFEST_BYTES)?)?;
    let wanted = parent
        .samples
        .iter()
        .map(|sample| sample.id.as_str())
        .collect::<BTreeSet<_>>();
    let truth = load_truth(&args.truth, &wanted)?;
    let truth_path = absolute(&args.truth);
    let derived = panel::derive(
        &parent,
        &panel::Source {
            path: &truth_path,
            sha256: &truth.sha256,
            database: &truth.database,
        },
        &args.version,
        &args.created,
    );
    let write_error = |source| Error::Write {
        path: args.output.clone(),
        source,
    };
    File::options()
        .write(true)
        .create_new(true)
        .open(&args.output)
        .and_then(|mut file| file.write_all(derived.text.as_bytes()))
        .map_err(write_error)?;
    eprintln!(
        "dna-eval: {} rows, {} newly retired, {} with another batch in the database",
        derived.rows, derived.retired, derived.batch_differs
    );
    if !derived.stranded.is_empty() {
        eprintln!(
            "dna-eval: no staying sample carries weight in {}; their retired weight is lost",
            derived.stranded.join(", ")
        );
    }
    let name = args
        .output
        .file_name()
        .map_or_else(String::new, |name| name.to_string_lossy().into_owned());
    Ok(format!("{}  {name} (body)", derived.body_sha256))
}

/// The panel's truth records, the database's hash and `info`, and the load
/// time.
pub(crate) struct Truth {
    pub(crate) sha256: String,
    pub(crate) database: truth::TruthDatabase,
    pub(crate) seconds: f64,
}

fn load_truth(path: &Path, wanted: &BTreeSet<&str>) -> Result<Truth> {
    let started = Instant::now();
    check_size(path, "truth database", truth::MAX_TRUTH_BYTES)?;
    let (sha256, database) = thread::scope(|scope| {
        let sha256 = scope.spawn(|| {
            runner::sha256_file(path).map_err(|source| Error::Read {
                kind: "truth database",
                path: path.to_path_buf(),
                source,
            })
        });
        let database = truth::load(path, wanted);
        (sha256.join().map_err(|_| Error::Worker), database)
    });
    Ok(Truth {
        sha256: sha256??,
        database: database?,
        seconds: started.elapsed().as_secs_f64(),
    })
}

/// Reads a whole input file of at most `maximum` bytes.
pub(crate) fn read(path: &Path, kind: &'static str, maximum: u64) -> Result<Vec<u8>> {
    check_size(path, kind, maximum)?;
    fs::read(path).map_err(|source| Error::Read {
        kind,
        path: path.to_path_buf(),
        source,
    })
}

/// Fails unless `path` is readable metadata of at most `maximum` bytes.
fn check_size(path: &Path, kind: &'static str, maximum: u64) -> Result<()> {
    let bytes = fs::metadata(path)
        .map_err(|source| Error::Read {
            kind,
            path: path.to_path_buf(),
            source,
        })?
        .len();
    if bytes > maximum {
        return Err(Error::TooLarge {
            kind,
            path: path.to_path_buf(),
            bytes,
            maximum,
        });
    }
    Ok(())
}

fn count(rows: &[Row<'_>], measure: Measure, status: Status) -> usize {
    rows.iter()
        .filter(|row| row.status(measure) == status)
        .count()
}

/// When the run started and how long DNA took, in seconds.
#[derive(Debug, Clone, Copy)]
struct Times {
    started: Instant,
    dna: f64,
}

fn header(
    cli: &RunArgs,
    manifest: &Manifest,
    truth: &Truth,
    report: &Report<'_>,
    outcomes: &[Outcome],
    times: Times,
) -> RunHeader {
    let (rows, check) = (report.rows(), report.check());
    let (all, current) = report.passes();
    let counts = |measure: Measure| Counts {
        total: rows.len(),
        concordant: count(rows, measure, Status::Concordant),
        discordant: count(rows, measure, Status::Discordant),
        not_applicable: count(rows, measure, Status::NotApplicable),
        failed: count(rows, measure, Status::Failed),
    };
    let mut millis = outcomes
        .iter()
        .map(|outcome| outcome.millis)
        .collect::<Vec<_>>();
    millis.sort_unstable();
    let seconds = |index: usize| {
        millis
            .get(index)
            .map_or(0.0, |value| *value as f64 / 1000.0)
    };
    let identities = rows
        .iter()
        .filter_map(|row| row.outcome.result.as_ref().ok())
        .map(|called| called.identity.clone())
        .collect::<BTreeSet<_>>();
    RunHeader {
        harness: env!("CARGO_PKG_NAME").to_owned(),
        harness_version: env!("CARGO_PKG_VERSION").to_owned(),
        source: source(),
        threads: cli.dna.threads,
        manifest: ManifestHeader {
            path: absolute(&cli.manifest),
            fields: manifest.header.clone(),
            body_sha256: manifest.body_sha256.clone(),
            rows: manifest.samples.len(),
        },
        truth: TruthHeader {
            path: absolute(&cli.truth),
            sha256: truth.sha256.clone(),
            info: truth.database.info.clone(),
            digest: check.current.clone(),
            recorded_digest: check.recorded.clone(),
            state: check.state().label().to_owned(),
            changed_samples: check.changed.len(),
            panel_samples_found: truth.database.samples.len(),
        },
        configuration_path: absolute(&cli.dna.config),
        reference_path: absolute(&cli.dna.reference),
        identities: identities.into_iter().collect(),
        review_window: cli.review.review_window,
        run_signal: cli.run_signal,
        run_signal_batch_reads: (cli.run_signal == Some(SignalGrouping::Batch))
            .then_some(cli.run_signal_batch_reads),
        current_since: cli.current_since,
        forced_windows: cli.forced_windows.as_deref().map(absolute),
        recheck: None,
        evidence_cache: None,
        forms_cache: None,
        placement_cache: None,
        run_context: None,
        constants: Constants::current(),
        zones: cli.review.header(),
        samples: Measure::ALL
            .into_iter()
            .map(|measure| (measure.label().to_owned(), counts(measure)))
            .collect(),
        workload: WorkloadHeader { all, current },
        seconds: Seconds {
            total: times.started.elapsed().as_secs_f64(),
            truth: truth.seconds,
            dna: times.dna,
            sample_median: seconds(millis.len() / 2),
            sample_p95: seconds(millis.len() * 95 / 100),
            sample_max: seconds(millis.len().saturating_sub(1)),
            sample_sum: millis.iter().sum::<u128>() as f64 / 1000.0,
        },
        outputs: BTreeMap::new(),
    }
}

/// The deterministic part of the header, for the summary.
fn context(header: &RunHeader) -> Vec<(&'static str, String)> {
    let field = |key: &str| {
        header
            .manifest
            .fields
            .get(key)
            .map_or("-", String::as_str)
            .to_owned()
    };
    let identity = header.identities.first();
    vec![
        (
            "panel",
            format!(
                "{} {} ({} rows, body sha256 {})",
                field("panel"),
                field("version"),
                header.manifest.rows,
                header.manifest.body_sha256
            ),
        ),
        (
            "configuration sha256",
            identity
                .map_or("-", |identity| identity.configuration_sha256.as_str())
                .to_owned(),
        ),
        (
            "reference",
            identity.map_or_else(
                || "-".to_owned(),
                |identity| {
                    format!(
                        "{} ({})",
                        identity.reference_name, identity.reference_sha256
                    )
                },
            ),
        ),
        (
            "profile",
            identity.map_or_else(
                || "-".to_owned(),
                |identity| format!("{} ({})", identity.profile_id, identity.profile_sha256),
            ),
        ),
    ]
    .into_iter()
    .chain(header.recheck.as_ref().map(|recheck| {
        (
            "re-check list",
            format!("{} rows, sha256 {}", recheck.rows, recheck.sha256),
        )
    }))
    .collect()
}

/// The git revision of the checkout the harness was built from, recorded by
/// the build script.
pub(crate) fn source() -> Source {
    Source {
        revision: option_env!("DNA_EVAL_REVISION")
            .filter(|revision| !revision.is_empty())
            .map(str::to_owned),
        dirty: match option_env!("DNA_EVAL_DIRTY") {
            Some("1") => Some(true),
            Some("0") => Some(false),
            _ => None,
        },
    }
}

/// The path made absolute against the working directory, or as given when
/// that fails.
pub(crate) fn absolute(path: &Path) -> PathBuf {
    std::path::absolute(path).unwrap_or_else(|_| path.to_path_buf())
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use super::*;

    fn panel(name: Option<&str>) -> Manifest {
        Manifest {
            header: name
                .map(|name| BTreeMap::from([("panel".to_owned(), name.to_owned())]))
                .unwrap_or_default(),
            body_sha256: String::new(),
            columns: Vec::new(),
            samples: Vec::new(),
        }
    }

    #[test]
    fn run_reads_development_readiness_and_regression_panels_only() {
        for name in ["dev", "readiness", "regress"] {
            assert!(runnable(&panel(Some(name))).is_ok(), "{name}");
        }
        assert!(matches!(runnable(&panel(None)), Err(Error::UnnamedPanel)));
        assert!(matches!(
            runnable(&panel(Some("validation"))),
            Err(Error::HoldOut { panel }) if panel == "validation"
        ));
    }
}
