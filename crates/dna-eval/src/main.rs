//! `dna-eval`: runs DNA over a test panel and compares it with reviewer truth
//! (PROP-0005, ADR-0075).
//!
//! `dna-eval run` reads a panel manifest, runs the recommended `consensus` then
//! `call` workflow for every sample in process, compares the calls with the
//! reviewer truth of the lab's statistics database, and writes the
//! comparisons, weighted estimates, anchors, and a summary into a new output
//! directory. `dna-eval diff` compares two runs and gates the candidate;
//! `dna-eval panel` derives a new panel version against the database.
//! Nothing in DNA depends on this crate.

mod bench;
mod cache;
mod cause;
mod classes;
mod compare;
mod components;
#[cfg(test)]
mod contract;
mod error;
mod estimate;
mod haplotype;
mod lod;
mod manifest;
mod panel;
mod recheck;
mod report;
mod rundiff;
mod rundir;
mod runner;
mod simbench;
mod subset;
mod tracts;
mod truth;
mod unrecorded;
mod workload;
mod zones;

use std::collections::{BTreeMap, BTreeSet};
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
use crate::manifest::Manifest;
use crate::recheck::Recheck;
use crate::report::{Report, Row, Status};
use crate::rundir::{
    Constants, Counts, ManifestHeader, RecheckHeader, RunContextHeader, RunHeader, Seconds, Source,
    TruthHeader, WorkloadHeader, ZoneHeader, prepare_output, write,
};
use crate::runner::{Outcome, Threads};
use crate::truth::TruthCheck;

/// The panels the harness reads: the development panel the rules are tuned
/// on, and the readiness panels that test them. A hold-out panel stays
/// reserved for an approval path not yet designed.
const RUNNABLE_PANELS: [&str; 2] = ["dev", "readiness"];

/// Largest accepted panel manifest.
const MAX_MANIFEST_BYTES: u64 = 256 * 1024 * 1024;
/// Largest accepted `forms.tsv`.
const MAX_FORMS_BYTES: u64 = 256 * 1024 * 1024;
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
    /// Runs DNA on the adjudicated artefact cases of a case list and checks
    /// each expected outcome: exits 2 with `BENCH FAIL` when a guarded case
    /// regresses or its sample cannot be analysed.
    Bench(BenchArgs),
    /// Runs DNA's read evidence stage on the frozen simulated tract-decode
    /// cases (PROP-0015) and reports the decoder's measures per tract and
    /// noise level, against a baseline when one is given.
    SimBench(SimBenchArgs),
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
    /// DNA configuration, as given to `dna` through `DNA_CONFIG`.
    #[arg(long)]
    config: PathBuf,
    /// Reference FASTA, as given to `dna --reference`.
    #[arg(long)]
    reference: PathBuf,
    /// Worker threads that run DNA (1 to 256).
    #[arg(long)]
    threads: usize,
    /// A read evidence cache, as for `run --evidence-cache`.
    #[arg(long)]
    evidence_cache: Option<PathBuf>,
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
    let threads = Threads::new(args.threads)?;
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
    let reference = dna_kernel::reference::load(&args.reference, ReferenceTopology::Linear)
        .map_err(|error| Error::Reference {
            path: args.reference.clone(),
            message: error.to_string(),
        })?;
    let cache = args
        .evidence_cache
        .as_deref()
        .map(|directory| {
            Ok::<_, cache::CacheError>(cache::Caches {
                evidence: cache::EvidenceCache::open(directory, &args.config)?,
                placed: cache::EvidenceCache::open_placed(
                    directory,
                    &args.config,
                    &reference.sequence_sha256,
                )?,
            })
        })
        .transpose()?;
    let outcomes = runner::run(
        &samples,
        (&args.reference, &reference.sequence_sha256),
        &args.config,
        threads,
        None,
        (cache.as_ref(), &runner::Contexts::new()),
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
    let revision = format!(
        "{}{}",
        source.revision.as_deref().unwrap_or("unknown revision"),
        match source.dirty {
            Some(true) => " (dirty)",
            Some(false) => "",
            None => " (dirty state unknown)",
        }
    );
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
    /// DNA configuration, as given to `dna` through `DNA_CONFIG`.
    #[arg(long)]
    config: PathBuf,
    /// Reference FASTA, as given to `dna --reference`.
    #[arg(long)]
    reference: PathBuf,
    /// New or empty output directory.
    #[arg(long)]
    output: PathBuf,
    /// Worker threads that run DNA (1 to 256).
    #[arg(long)]
    threads: usize,
    /// An earlier run's `forms.tsv` for the same truth database, read
    /// instead of indexing the whole database.
    #[arg(long)]
    forms: Option<PathBuf>,
    /// Window of the review workload, in positions, aligned to each HV
    /// region's first position.
    #[arg(long, default_value_t = report::REVIEW_WINDOW)]
    review_window: i64,
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
    current_since: Option<String>,
    /// Judge each review against a run signal (PROP-0009) built in process
    /// from the panel: `batch` takes a sample's batch when it has at least
    /// `--run-signal-batch-reads` reads, else its batch month; `month` takes
    /// the month alone. The run takes two passes.
    #[arg(long, value_parser = ["batch", "month"])]
    run_signal: Option<String>,
    /// Fewest reads of a batch for `--run-signal batch` to use it.
    #[arg(long, default_value_t = 40)]
    run_signal_batch_reads: usize,
    /// A read evidence cache: every read is analysed from its
    /// `dna.read_evidence/v1` document, read from the directory or made and
    /// stored there, so a run whose read stage is unchanged skips it. It
    /// needs a harness built from committed read-stage sources; `run.json`
    /// records its key and hits.
    #[arg(long)]
    evidence_cache: Option<PathBuf>,
    /// A directory of run contexts (PROP-0014), one
    /// `<batch>.run_context.json` per batch, as `dna run-context` writes them
    /// over the batches' full folders: a sample of a batch with one has its
    /// consensus calibrated by it. `run_context.tsv` reports per batch month
    /// how many samples a context calibrated.
    #[arg(long)]
    run_context: Option<PathBuf>,
    /// Review zones: positions on each side of a repeat window that follow
    /// it.
    #[arg(long, default_value_t = 20)]
    zone_after: i64,
    /// Review zones: positions at each end of a region that are its
    /// amplicon's start or end.
    #[arg(long, default_value_t = 40)]
    zone_edge: i64,
    /// Review zones: positions past a repeat window that an item overlapping
    /// the window still counts for the window.
    #[arg(long, default_value_t = 2)]
    zone_overrun: i64,
}

fn main() -> ExitCode {
    let result = match Cli::parse().command {
        Command::Diff(args) => return diff(&args),
        Command::Bench(args) => return bench(&args),
        Command::SimBench(args) => sim_bench(&args),
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

/// Whether the harness may read the panel: one its header does not name, or
/// a development or readiness panel.
fn runnable(manifest: &Manifest) -> Result<()> {
    match manifest.header.get("panel") {
        Some(panel) if !RUNNABLE_PANELS.contains(&panel.as_str()) => Err(Error::HoldOut {
            panel: panel.clone(),
        }),
        _ => Ok(()),
    }
}

/// Whether `value` is a month written `YYYY-MM`.
fn is_month(value: &str) -> bool {
    value.len() == 7
        && value.bytes().enumerate().all(|(index, byte)| {
            if index == 4 {
                byte == b'-'
            } else {
                byte.is_ascii_digit()
            }
        })
        && matches!(
            &value[5..],
            "01" | "02" | "03" | "04" | "05" | "06" | "07" | "08" | "09" | "10" | "11" | "12"
        )
}

fn run(cli: &RunArgs) -> Result<String> {
    let started = Instant::now();
    let threads = Threads::new(cli.threads)?;
    let manifest = manifest::parse(&read(&cli.manifest, "manifest", MAX_MANIFEST_BYTES)?)?;
    runnable(&manifest)?;

    // The whole truth of the panel loads before DNA runs, so a database that
    // cannot be opened, has another schema, or fails a query stops the run
    // before any sample is called.
    let wanted = manifest
        .samples
        .iter()
        .map(|sample| sample.id.as_str())
        .collect::<BTreeSet<_>>();
    let truth = load_truth(&cli.truth, &wanted)?;
    // The haplotype measures apply alleles to this sequence; a sample DNA
    // analyses on another sequence fails.
    let reference = dna_kernel::reference::load(&cli.reference, ReferenceTopology::Linear)
        .map_err(|error| Error::Reference {
            path: cli.reference.clone(),
            message: error.to_string(),
        })?;
    // How the reviewers write each window sequence, over the whole truth,
    // for the cause of a representation discrepancy.
    let forms = if let Some(path) = &cli.forms {
        let text = read(path, "forms", MAX_FORMS_BYTES)?;
        let text = String::from_utf8(text).map_err(|error| Error::Forms {
            path: path.clone(),
            message: error.to_string(),
        })?;
        Forms::parse(&text, &truth.sha256).map_err(|message| Error::Forms {
            path: path.clone(),
            message,
        })?
    } else {
        let corpus = truth::sample_ids(&cli.truth)?;
        let corpus = truth::load(&cli.truth, &corpus.iter().map(String::as_str).collect())?;
        Forms::new(corpus.samples.values(), &reference.sequence)
    };
    if let Some(since) = &cli.current_since
        && !is_month(since)
    {
        return Err(Error::CurrentSince {
            value: since.clone(),
        });
    }
    if cli.review_window < 1 {
        return Err(Error::ReviewWindow {
            window: cli.review_window,
        });
    }
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
    let (recheck, recheck_header) = match &cli.recheck {
        Some(path) => {
            let failure = |message: String| Error::Recheck {
                path: path.clone(),
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
            (Some(list), Some(header))
        }
        None => (None, None),
    };
    let widths = zones::Widths {
        after: cli.zone_after,
        edge: cli.zone_edge,
        overrun: cli.zone_overrun,
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
    let profile = target_profile(&cli.config)?;
    prepare_output(&cli.output)?;
    let dna_started = Instant::now();
    let cache = cli
        .evidence_cache
        .as_deref()
        .map(|directory| {
            Ok::<_, cache::CacheError>(cache::Caches {
                evidence: cache::EvidenceCache::open(directory, &cli.config)?,
                placed: cache::EvidenceCache::open_placed(
                    directory,
                    &cli.config,
                    &reference.sequence_sha256,
                )?,
            })
        })
        .transpose()?;
    let contexts = cli
        .run_context
        .as_deref()
        .map(|directory| runner::contexts(directory, &manifest.samples))
        .transpose()?
        .unwrap_or_default();
    let mut signal_groups = Vec::new();
    let outcomes = runner::run(
        &manifest.samples,
        (&cli.reference, &reference.sequence_sha256),
        &cli.config,
        threads,
        cli.run_signal.as_deref().map(|grouping| match grouping {
            "month" => runner::Grouping::Month,
            _ => runner::Grouping::Batch {
                batch_reads: cli.run_signal_batch_reads,
            },
        }),
        (cache.as_ref(), &contexts),
        &mut signal_groups,
    )?;
    let dna_seconds = dna_started.elapsed().as_secs_f64();
    let analysed = &profile.identity().sha256;
    if let Some(other) = outcomes
        .iter()
        .filter_map(|outcome| outcome.result.as_ref().ok())
        .map(|called| &called.identity.profile_sha256)
        .find(|sha256| *sha256 != analysed)
    {
        return Err(Error::Profile {
            path: cli.config.clone(),
            message: format!("DNA analysed with profile {other}, not {analysed}"),
        });
    }
    let span = |first: usize, last: usize| {
        (
            i64::try_from(first).unwrap_or(i64::MAX),
            i64::try_from(last).unwrap_or(i64::MAX),
        )
    };
    let zones = zones::Zones::new(
        &profile
            .regions
            .iter()
            .map(|&[first, last]| span(first, last))
            .collect::<Vec<_>>(),
        &profile
            .windows
            .iter()
            .map(|window| span(window.start + 1, window.end()))
            .collect::<Vec<_>>(),
        widths,
    );

    let rows = report::rows(
        &manifest.samples,
        &outcomes,
        &truth.database.samples,
        &reference.sequence,
        &profile.conventions,
    );
    let check = TruthCheck::new(&manifest, &truth.database);
    let report = Report::new(&rows, &check, (&forms, reference.sequence.as_str()))
        .with_review(
            cli.review_window,
            forced,
            Some(zones),
            recheck,
            profile
                .length_tracts
                .iter()
                .map(|tract| workload::Tract {
                    name: tract.name.clone(),
                    first: i64::try_from(tract.window[0]).unwrap_or(i64::MAX),
                    last: i64::try_from(tract.window[1]).unwrap_or(i64::MAX),
                    base: char::from(tract.base),
                })
                .collect(),
        )
        .with_current_since(cli.current_since.clone());

    let mut header = header(
        cli,
        &manifest,
        &truth,
        &report,
        &outcomes,
        Times {
            started,
            dna: dna_seconds,
        },
    );
    header.recheck = recheck_header;
    header.evidence_cache = cache.as_ref().map(|cache| cache.evidence.usage());
    header.placement_cache = cache.as_ref().map(|cache| cache.placed.usage());
    header.run_context = cli
        .run_context
        .as_deref()
        .map(|directory| RunContextHeader {
            directory: absolute(directory),
            batches: contexts.len(),
        });
    for (name, contents) in report.files(&context(&header)) {
        write(&cli.output, name, &contents)?;
    }
    if cli.run_context.is_some() {
        write(
            &cli.output,
            rundir::RUN_CONTEXT,
            &rundir::context_coverage(&manifest.samples, &outcomes, &contexts),
        )?;
    }
    if cli.run_signal.is_some() {
        write(
            &cli.output,
            rundir::SIGNALS,
            &rundir::signals(&signal_groups),
        )?;
    }
    write(
        &cli.output,
        rundir::CHECKS,
        &rundir::checks(&manifest.samples, &outcomes),
    )?;
    write(
        &cli.output,
        rundir::NOTES,
        &rundir::notes(&manifest.samples, &outcomes),
    )?;
    write(
        &cli.output,
        rundir::COMPONENTS,
        &rundir::components(
            &manifest.samples,
            &outcomes,
            &truth.database.samples,
            reference.sequence.as_bytes(),
        ),
    )?;
    write(
        &cli.output,
        rundir::TRACTS,
        &rundir::tracts(&manifest.samples, &outcomes),
    )?;
    write(&cli.output, rundir::FORMS, &forms.table(&truth.sha256))?;
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
struct Truth {
    sha256: String,
    database: truth::TruthDatabase,
    seconds: f64,
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
fn read(path: &Path, kind: &'static str, maximum: u64) -> Result<Vec<u8>> {
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
        threads: cli.threads,
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
        configuration_path: absolute(&cli.config),
        reference_path: absolute(&cli.reference),
        identities: identities.into_iter().collect(),
        review_window: cli.review_window,
        run_signal: cli.run_signal.clone(),
        current_since: cli.current_since.clone(),
        forced_windows: cli.forced_windows.as_deref().map(absolute),
        recheck: None,
        evidence_cache: None,
        placement_cache: None,
        run_context: None,
        constants: Constants::current(),
        zones: ZoneHeader {
            after: cli.zone_after,
            edge: cli.zone_edge,
            overrun: cli.zone_overrun,
        },
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
fn source() -> Source {
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
fn absolute(path: &Path) -> PathBuf {
    std::path::absolute(path).unwrap_or_else(|_| path.to_path_buf())
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use super::*;

    #[test]
    fn current_data_months_are_year_dash_month() {
        assert!(is_month("2026-01") && is_month("2025-12"));
        assert!(!is_month("2026-1") && !is_month("2026-13") && !is_month("2026/01"));
    }

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
    fn run_reads_development_and_readiness_panels_only() {
        for name in [None, Some("dev"), Some("readiness")] {
            assert!(runnable(&panel(name)).is_ok(), "{name:?}");
        }
        assert!(matches!(
            runnable(&panel(Some("validation"))),
            Err(Error::HoldOut { panel }) if panel == "validation"
        ));
    }
}
