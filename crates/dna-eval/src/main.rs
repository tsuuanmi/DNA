//! `dna-eval`: runs DNA over a test panel and compares it with reviewer truth
//! (PROP-0005, ADR-0075).
//!
//! `dna-eval run` reads a panel manifest, runs the recommended `consensus` then
//! `call` workflow for every sample in process, compares the calls with the
//! reviewer truth of the lab's statistics database, and writes the
//! comparisons, weighted estimates, anchors, the diff against an optional
//! baseline run, and a summary into a new output directory; with `--triage`
//! it also scores the draft stage-1 triage of PROP-0004. `dna-eval panel`
//! derives a new panel version against the database. Nothing in DNA depends on
//! this crate.

mod baseline;
mod compare;
mod error;
mod estimate;
mod evidence;
mod manifest;
mod panel;
mod report;
mod runner;
mod triage;
mod truth;

use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File};
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::thread;
use std::time::Instant;

use clap::{Args, Parser, Subcommand};
use serde::Serialize;

use crate::compare::Measure;
use crate::error::{Error, Result};
use crate::manifest::Manifest;
use crate::report::{Report, Row, Status};
use crate::runner::{Identity, Outcome, Threads};
use crate::truth::TruthCheck;

/// Largest accepted panel manifest.
const MAX_MANIFEST_BYTES: u64 = 256 * 1024 * 1024;
/// Largest accepted baseline scores file.
const MAX_BASELINE_BYTES: u64 = 256 * 1024 * 1024;

/// Runs DNA over a test panel and compares it with reviewer truth.
#[derive(Debug, Parser)]
#[command(name = "dna-eval", version)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Runs DNA over a panel and compares it with the truth.
    Run(RunArgs),
    /// Derives a new panel version against the truth database.
    Panel(PanelArgs),
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
    /// Earlier run directory whose `scores.tsv` is the baseline.
    #[arg(long)]
    baseline: Option<PathBuf>,
    /// Worker threads that run DNA (1 to 256).
    #[arg(long)]
    threads: usize,
    /// Also run the draft stage-1 triage of PROP-0004 phase 0 and write its
    /// `triage*` outputs; research only.
    #[arg(long)]
    triage: bool,
}

fn main() -> ExitCode {
    let result = match Cli::parse().command {
        Command::Run(args) => run(&args),
        Command::Panel(args) => derive_panel(&args),
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

fn run(cli: &RunArgs) -> Result<String> {
    let started = Instant::now();
    let threads = Threads::new(cli.threads)?;
    let manifest = manifest::parse(&read(&cli.manifest, "manifest", MAX_MANIFEST_BYTES)?)?;
    let baseline = cli
        .baseline
        .as_ref()
        .map(|directory| -> Result<_> {
            let path = directory.join(baseline::SCORES);
            Ok(baseline::parse(&read(
                &path,
                "baseline scores",
                MAX_BASELINE_BYTES,
            )?)?)
        })
        .transpose()?;

    // The whole truth of the panel loads before DNA runs, so a database that
    // cannot be opened, has another schema, or fails a query stops the run
    // before any sample is called.
    let wanted = manifest
        .samples
        .iter()
        .map(|sample| sample.id.as_str())
        .collect::<BTreeSet<_>>();
    let truth = load_truth(&cli.truth, &wanted)?;
    prepare_output(&cli.output)?;
    let dna_started = Instant::now();
    let outcomes = runner::run(
        &manifest.samples,
        &cli.reference,
        &cli.config,
        threads,
        cli.triage,
    )?;
    let dna_seconds = dna_started.elapsed().as_secs_f64();

    let rows = report::rows(&manifest.samples, &outcomes, &truth.database.samples);
    let check = TruthCheck::new(&manifest, &truth.database);
    let report = Report::new(&rows, baseline.as_ref(), &check);
    let triage_started = Instant::now();
    let triage = cli
        .triage
        .then(|| triage::score::Triage::new(&rows, triage::score::sweep()));
    let triage_files = triage.as_ref().map(triage::score::Triage::files);
    let triage_seconds = triage_started.elapsed().as_secs_f64();

    let mut header = header(
        cli,
        &manifest,
        &truth,
        &check,
        &rows,
        &outcomes,
        Times {
            started,
            dna: dna_seconds,
        },
    );
    header.triage = triage.as_ref().map(|triage| TriageHeader {
        policy: triage
            .policy()
            .map(|policy| format!("{}: {}", policy.name, policy.parameters()))
            .unwrap_or_default(),
        policies: triage.policies(),
        seconds: triage_seconds,
    });
    for (name, contents) in report
        .files(&context(&header))
        .into_iter()
        .chain(triage_files.into_iter().flatten())
    {
        write(&cli.output, name, &contents)?;
    }
    let mut json = serde_json::to_string_pretty(&header).map_err(Error::Serialize)?;
    json.push('\n');
    write(&cli.output, "run.json", &json)?;

    let count = |measure: Measure, status: Status| count(&rows, measure, status);
    Ok(format!(
        "dna-eval: truth {} ({} samples changed); {} samples; primary {} concordant, {} discordant, {} not applicable, {} \
         failed; strict {} concordant; IUPAC-compatible {} concordant; {:.1} s; outputs in {}",
        check.state().label(),
        check.changed.len(),
        rows.len(),
        count(Measure::Primary, Status::Concordant),
        count(Measure::Primary, Status::Discordant),
        count(Measure::Primary, Status::NotApplicable),
        count(Measure::Primary, Status::Failed),
        count(Measure::Strict, Status::Concordant),
        count(Measure::IupacCompatible, Status::Concordant),
        started.elapsed().as_secs_f64(),
        cli.output.display()
    ))
}

/// `dna-eval panel`: writes the new manifest and returns its checksum line.
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

/// Creates the output directory, which must not hold files yet.
fn prepare_output(directory: &Path) -> Result<()> {
    let failure = |source| Error::Write {
        path: directory.to_path_buf(),
        source,
    };
    fs::create_dir_all(directory).map_err(failure)?;
    if fs::read_dir(directory).map_err(failure)?.next().is_some() {
        return Err(Error::OutputNotEmpty {
            path: directory.to_path_buf(),
        });
    }
    Ok(())
}

fn write(directory: &Path, name: &str, contents: &str) -> Result<()> {
    let path = directory.join(name);
    fs::write(&path, contents).map_err(|source| Error::Write { path, source })
}

/// `run.json`: what was run, on what, and how long it took.
#[derive(Debug, Serialize)]
struct Header {
    harness: &'static str,
    harness_version: &'static str,
    source: Source,
    threads: usize,
    manifest: ManifestHeader,
    truth: TruthHeader,
    configuration_path: PathBuf,
    reference_path: PathBuf,
    /// Distinct identities DNA reported; one unless the inputs changed during
    /// the run.
    identities: Vec<Identity>,
    baseline: Option<PathBuf>,
    /// Outcome counts per measure.
    samples: BTreeMap<&'static str, Counts>,
    seconds: Seconds,
    /// The draft triage, when the run asked for it.
    #[serde(skip_serializing_if = "Option::is_none")]
    triage: Option<TriageHeader>,
}

/// The draft triage of a run: the reported policy, the policies swept, and
/// the time to decide and score every policy.
#[derive(Debug, Serialize)]
struct TriageHeader {
    policy: String,
    policies: usize,
    seconds: f64,
}

#[derive(Debug, Serialize)]
struct Source {
    /// Git revision of the checkout the harness was built from, if available.
    revision: Option<String>,
    /// Whether tracked files of that checkout differed from the revision.
    dirty: Option<bool>,
}

#[derive(Debug, Serialize)]
struct ManifestHeader {
    path: PathBuf,
    fields: BTreeMap<String, String>,
    body_sha256: String,
    rows: usize,
}

#[derive(Debug, Serialize)]
struct TruthHeader {
    path: PathBuf,
    /// The database file's SHA-256, for provenance only.
    sha256: String,
    /// The database's `info` table, for provenance only.
    info: BTreeMap<String, String>,
    /// The panel truth digest of the current database.
    digest: String,
    /// The panel truth digest the manifest records.
    recorded_digest: Option<String>,
    /// `unchanged`, `changed`, or `not_recorded`.
    state: &'static str,
    /// Samples whose truth differs from the manifest's record.
    changed_samples: usize,
    panel_samples_found: usize,
}

#[derive(Debug, Serialize)]
struct Counts {
    total: usize,
    concordant: usize,
    discordant: usize,
    not_applicable: usize,
    failed: usize,
}

fn count(rows: &[Row<'_>], measure: Measure, status: Status) -> usize {
    rows.iter()
        .filter(|row| row.status(measure) == status)
        .count()
}

#[derive(Debug, Serialize)]
struct Seconds {
    total: f64,
    truth: f64,
    dna: f64,
    sample_median: f64,
    sample_p95: f64,
    sample_max: f64,
    sample_sum: f64,
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
    check: &TruthCheck,
    rows: &[Row<'_>],
    outcomes: &[Outcome],
    times: Times,
) -> Header {
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
    Header {
        harness: env!("CARGO_PKG_NAME"),
        harness_version: env!("CARGO_PKG_VERSION"),
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
            state: check.state().label(),
            changed_samples: check.changed.len(),
            panel_samples_found: truth.database.samples.len(),
        },
        configuration_path: absolute(&cli.config),
        reference_path: absolute(&cli.reference),
        identities: identities.into_iter().collect(),
        baseline: cli.baseline.as_deref().map(absolute),
        samples: Measure::ALL
            .into_iter()
            .map(|measure| (measure.label(), counts(measure)))
            .collect(),
        seconds: Seconds {
            total: times.started.elapsed().as_secs_f64(),
            truth: truth.seconds,
            dna: times.dna,
            sample_median: seconds(millis.len() / 2),
            sample_p95: seconds(millis.len() * 95 / 100),
            sample_max: seconds(millis.len().saturating_sub(1)),
            sample_sum: millis.iter().sum::<u128>() as f64 / 1000.0,
        },
        triage: None,
    }
}

/// The deterministic part of the header, for the summary.
fn context(header: &Header) -> Vec<(&'static str, String)> {
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
        (
            "baseline",
            header
                .baseline
                .as_ref()
                .map_or_else(|| "none".to_owned(), |path| path.display().to_string()),
        ),
    ]
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
