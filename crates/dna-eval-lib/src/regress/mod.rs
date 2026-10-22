//! `dna-eval regress`: the versioned regression set.
//!
//! `check` runs DNA on every row of a `panel=regress` manifest exactly as
//! `run` does, writes its golden output (`golden.tsv`), and compares it with
//! the head of the version ledger: every changed line must be declared by
//! the pull request's committed intent file, a truth-free change needs a
//! new version key, and the forbidden rules fail whatever the intents.
//! `record` appends a checked run of a clean `main` build to the ledger: a
//! copy under `results/` and one row in each committed baseline file. The
//! outputs are `dna.eval_regress/v1`
//! (`docs/reference/schemas/eval-regress-v1.schema.json`).

mod golden;
mod intents;
mod judge;
mod ledger;

use std::collections::BTreeSet;
use std::fmt::Write as _;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

use clap::Subcommand;
use serde::{Deserialize, Serialize};

use crate::cache::Usage;
use crate::error::{Error, Result};
use crate::geometry::{Geometry, Span};
use crate::identity::Identity;
use crate::manifest::{self, Manifest};
use crate::report::{self, Report};
use crate::rundir::{
    ManifestHeader, RecheckHeader, RunContextHeader, Source, ZoneHeader, prepare_output, write,
};
use crate::truth::{TruthCheck, TruthState};
use crate::zones::{Widths, Zones};
use crate::{DnaOptions, Loaded, MAX_MANIFEST_BYTES, Ran, ReviewOptions, absolute, read, workload};
use golden::{Counts, Golden, Kind, Sources};
use intents::{Change, Intents};
use judge::{Context, Judgement, Label, Verdict};

/// Largest accepted intent file.
const MAX_INTENT_BYTES: u64 = 1024 * 1024;
/// Largest accepted golden output, run header, or baseline file.
const MAX_LEDGER_BYTES: u64 = 512 * 1024 * 1024;
/// The branch whose current commit a recorded run must have been built from.
const MAIN: &str = "origin/main";
/// The changed lines `regress.md` lists one by one.
const LISTED_LINES: usize = 40;

/// The files of a check's output directory, which `record` copies.
const GOLDEN: &str = "golden.tsv";
const CHANGES: &str = "changes.tsv";
const REPORT: &str = "regress.md";
const RUN_JSON: &str = "run.json";
/// The intent file's name in a ledger entry.
const INTENT: &str = "intent.tsv";

/// Arguments of `dna-eval regress`.
#[derive(Debug, clap::Args)]
pub(crate) struct Args {
    #[command(subcommand)]
    action: Action,
}

#[derive(Debug, Subcommand)]
enum Action {
    /// Runs DNA on every row of the regression set and compares its golden
    /// output with the ledger head; exits 2 with `REGRESS FAIL` when a rule
    /// fires. It writes nothing to the ledger.
    Check(CheckArgs),
    /// Appends a checked run of a clean build of `main` to the ledger and
    /// one row to each committed baseline file. It never runs DNA.
    Record(RecordArgs),
}

/// Arguments of `dna-eval regress check`.
#[derive(Debug, clap::Args)]
struct CheckArgs {
    /// The regression set's manifest (`panel=regress`).
    #[arg(long)]
    manifest: PathBuf,
    /// The pinned truth database (`statistics.sqlite`), opened read-only.
    #[arg(long)]
    truth: PathBuf,
    /// The frozen re-check list, as for `run --recheck`.
    #[arg(long)]
    recheck: PathBuf,
    #[command(flatten)]
    dna: DnaOptions,
    #[command(flatten)]
    review: ReviewOptions,
    /// The ledger directory (`results/regress/ledger`); its last
    /// `index.tsv` entry is the head. Absent or empty, the check is the
    /// series' first.
    #[arg(long)]
    ledger: PathBuf,
    /// The pull request's intent file, committed and unmodified at the
    /// harness's revision.
    #[arg(long)]
    intent: Option<PathBuf>,
    /// Accept another truth digest or re-check list than the head's: the
    /// truth-free lines must not change.
    #[arg(long)]
    rebase: bool,
    /// A new directory for `golden.tsv`, `changes.tsv`, `regress.md`, and
    /// `run.json`.
    #[arg(long)]
    output: PathBuf,
}

/// Arguments of `dna-eval regress record`.
#[derive(Debug, clap::Args)]
struct RecordArgs {
    /// The output directory of `regress check`.
    #[arg(long)]
    run: PathBuf,
    /// The ledger directory.
    #[arg(long)]
    ledger: PathBuf,
    /// The new entry, `r` and four or more digits, for example `r0001`.
    #[arg(long)]
    entry: String,
    /// The committed ledger baseline (`regress-v1-ledger.tsv`).
    #[arg(long)]
    baseline: PathBuf,
    /// The committed strata baseline (`regress-v1-strata.tsv`).
    #[arg(long)]
    strata: PathBuf,
}

/// Runs `dna-eval regress`: exits 0 for `SAME`, `EXPLAINED`, `INITIAL`,
/// `REBASE`, and a recorded entry, 2 for `FAIL`, and 1 for an input error.
pub(crate) fn main(args: &Args, identity: &Identity) -> ExitCode {
    let result = match &args.action {
        Action::Check(args) => check(args, &System, identity).map(|(report, verdict)| {
            print!("{report}");
            verdict
        }),
        Action::Record(args) => record(args, &System, &today()).map(|message| {
            println!("{message}");
            Verdict::Same
        }),
    };
    match result {
        Ok(Verdict::Fail) => ExitCode::from(2),
        Ok(_) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("dna-eval: {error}");
            ExitCode::FAILURE
        }
    }
}

/// The version key of a build: what must move when its truth-free output
/// does ([[governance/versioning#Plugins]]).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct VersionKey {
    /// The workspace version.
    pub(crate) dna_version: String,
    /// Every plugin of the static registry, `id=version` sorted by id and
    /// comma-joined.
    pub(crate) plugins: String,
    /// The result contracts DNA emits, sorted and comma-joined.
    pub(crate) contracts: String,
    /// The configuration file's SHA-256.
    pub(crate) configuration_sha256: String,
}

impl VersionKey {
    /// The key of the build `identity` records under the configuration of
    /// `configuration_sha256`.
    fn current(configuration_sha256: &str, identity: &Identity) -> Self {
        let mut plugins = dna::plugins();
        plugins.sort_unstable();
        Self {
            dna_version: identity.version.to_owned(),
            plugins: plugins
                .iter()
                .map(|(id, version)| format!("{id}={version}"))
                .collect::<Vec<_>>()
                .join(","),
            contracts: identity.contracts.to_owned(),
            configuration_sha256: configuration_sha256.to_owned(),
        }
    }

    /// The fields that differ from `other`'s.
    fn moved(&self, other: &Self) -> Vec<&'static str> {
        [
            ("dna_version", self.dna_version != other.dna_version),
            ("plugins", self.plugins != other.plugins),
            ("contracts", self.contracts != other.contracts),
            (
                "configuration_sha256",
                self.configuration_sha256 != other.configuration_sha256,
            ),
        ]
        .into_iter()
        .filter_map(|(name, moved)| moved.then_some(name))
        .collect()
    }
}

/// The truth a check measured against.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct TruthIdentity {
    pub(crate) path: PathBuf,
    /// The database file's SHA-256, for provenance.
    pub(crate) sha256: String,
    /// The panel truth digest of the database.
    pub(crate) digest: String,
    /// The panel truth digest the manifest records.
    pub(crate) recorded_digest: Option<String>,
    /// `unchanged`, `changed`, or `not_recorded`.
    pub(crate) state: String,
    pub(crate) changed_samples: usize,
}

/// The intent file a check read.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct IntentHeader {
    pub(crate) path: PathBuf,
    /// Its git blob id, committed at the harness's revision.
    pub(crate) blob: String,
    pub(crate) pr: String,
    pub(crate) rows: usize,
    pub(crate) posthoc_rows: usize,
}

/// One kind's digest.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct KindDigest {
    pub(crate) kind: String,
    pub(crate) lines: usize,
    pub(crate) sha256: String,
}

/// One stratum's digest and changed samples.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct StratumDigest {
    pub(crate) stratum: String,
    pub(crate) samples: usize,
    pub(crate) golden_sha256: String,
    pub(crate) changed_samples: usize,
}

/// The change counts against the head; zero for a first check.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct ChangeCounts {
    pub(crate) changed_samples: usize,
    pub(crate) changed_lines: usize,
    pub(crate) explained_lines: usize,
    pub(crate) posthoc_rows: usize,
}

/// Wall times, in seconds.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Seconds {
    pub(crate) total: f64,
    pub(crate) truth: f64,
    pub(crate) dna: f64,
}

/// `run.json` of a check: what was checked, against what, and the verdict.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Header {
    pub(crate) harness: String,
    pub(crate) harness_version: String,
    pub(crate) source: Source,
    pub(crate) threads: usize,
    /// `regress-<manifest version>`.
    pub(crate) series: String,
    pub(crate) manifest: ManifestHeader,
    pub(crate) truth: TruthIdentity,
    pub(crate) recheck: Option<RecheckHeader>,
    pub(crate) configuration_path: PathBuf,
    pub(crate) reference_path: PathBuf,
    pub(crate) profile_sha256: String,
    pub(crate) reference_sha256: String,
    pub(crate) review_window: i64,
    pub(crate) zones: ZoneHeader,
    pub(crate) run_context: Option<RunContextHeader>,
    pub(crate) version_key: VersionKey,
    /// The workload counting version ([`workload::COUNTING_VERSION`]) the
    /// evaluation lines follow. Absent in the checks recorded before it was
    /// kept, whose version [`counting_version`] reads from their revision.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) workload_version: Option<u32>,
    pub(crate) evidence_cache: Option<Usage>,
    pub(crate) placement_cache: Option<Usage>,
    /// The head entry compared with; none for the series' first check.
    pub(crate) head: Option<String>,
    pub(crate) intent: Option<IntentHeader>,
    pub(crate) rebase: bool,
    pub(crate) verdict: String,
    /// The rules that fired.
    pub(crate) rules: Vec<String>,
    pub(crate) warnings: Vec<String>,
    pub(crate) golden_sha256: String,
    pub(crate) kinds: Vec<KindDigest>,
    pub(crate) strata: Vec<StratumDigest>,
    pub(crate) counts: Counts,
    pub(crate) changes: ChangeCounts,
    pub(crate) seconds: Seconds,
}

impl Header {
    fn text(&self) -> Result<String> {
        let mut json = serde_json::to_string_pretty(self).map_err(Error::Serialize)?;
        json.push('\n');
        Ok(json)
    }
}

/// The git operations the ledger rests on.
trait Git {
    /// The blob id of the file at `path` as it is on disk.
    fn blob(&self, path: &Path) -> std::result::Result<String, String>;
    /// The blob id of the tracked file at `path` at `revision`.
    fn committed(&self, path: &Path, revision: &str) -> std::result::Result<String, String>;
    /// The commit [`MAIN`] names in the repository holding `directory`.
    fn main_revision(&self, directory: &Path) -> std::result::Result<String, String>;
}

/// The `git` command.
struct System;

impl System {
    fn git(directory: &Path, arguments: &[&str]) -> std::result::Result<String, String> {
        let output = Command::new("git")
            .arg("-C")
            .arg(directory)
            .args(arguments)
            .output()
            .map_err(|error| format!("cannot run git: {error}"))?;
        if output.status.success() {
            Ok(String::from_utf8_lossy(&output.stdout).trim().to_owned())
        } else {
            Err(format!(
                "git {}: {}",
                arguments.join(" "),
                String::from_utf8_lossy(&output.stderr).trim()
            ))
        }
    }

    fn split(path: &Path) -> std::result::Result<(&Path, &str), String> {
        let name = path
            .file_name()
            .and_then(|name| name.to_str())
            .ok_or_else(|| format!("{} names no file", path.display()))?;
        let directory = path
            .parent()
            .filter(|parent| !parent.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        Ok((directory, name))
    }
}

impl Git for System {
    fn blob(&self, path: &Path) -> std::result::Result<String, String> {
        let (directory, name) = Self::split(path)?;
        Self::git(directory, &["hash-object", "--", name])
    }

    fn committed(&self, path: &Path, revision: &str) -> std::result::Result<String, String> {
        let (directory, name) = Self::split(path)?;
        let tracked = Self::git(
            directory,
            &["ls-files", "--full-name", "--error-unmatch", "--", name],
        )
        .map_err(|_| format!("{} is not tracked by git", path.display()))?;
        Self::git(directory, &["rev-parse", &format!("{revision}:{tracked}")])
            .map_err(|_| format!("{} is not committed at {revision}", path.display()))
    }

    fn main_revision(&self, directory: &Path) -> std::result::Result<String, String> {
        let directory = if directory.as_os_str().is_empty() {
            Path::new(".")
        } else {
            directory
        };
        Self::git(
            directory,
            &[
                "rev-parse",
                "--verify",
                "--quiet",
                &format!("{MAIN}^{{commit}}"),
            ],
        )
    }
}

/// Today's UTC date.
fn today() -> String {
    ledger::date(
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map_or(0, |elapsed| elapsed.as_secs()),
    )
}

/// An input error of `regress`.
fn refused(message: impl Into<String>) -> Error {
    Error::Regress(message.into())
}

/// Accepts a regression-set manifest only.
fn regress_panel(manifest: &Manifest) -> Result<()> {
    match manifest.header.get("panel").map(String::as_str) {
        Some("regress") => Ok(()),
        Some(panel) => Err(refused(format!(
            "the manifest is panel {panel:?}; regress reads panel=regress manifests only"
        ))),
        None => Err(refused(
            "the manifest names no panel; regress reads panel=regress manifests only",
        )),
    }
}

/// The workload counting version of each check recorded before `run.json`
/// held it, by the revision whose build ran it: regress-v1's `r0001` to
/// `r0003` counted with version 1 (no counting version yet), `r0004`, the
/// merge of the hardened gate, with version 2.
const UNRECORDED_WORKLOAD: [(&str, u32); 4] = [
    ("74c9cbaf2a3c4e5e25cfe83b090394bc9bd75959", 1),
    ("c2fcd7715f6d2e8f75fd90315e64d69a3dcb5123", 1),
    ("39a2c582f293e32b8489b3d02ae0cbff2378c932", 1),
    ("1a23f5751e256aa0760ff39dedf692464f2d3b6d", 2),
];

/// The workload counting version `header`'s evaluation lines follow: the
/// recorded one, else the one [`UNRECORDED_WORKLOAD`] names for its
/// revision; a check that records neither is refused.
fn counting_version(header: &Header) -> Result<u32> {
    if let Some(version) = header.workload_version {
        return Ok(version);
    }
    let revision = header.source.revision.as_deref();
    UNRECORDED_WORKLOAD
        .iter()
        .find(|(known, _)| Some(*known) == revision)
        .map(|(_, version)| *version)
        .ok_or_else(|| {
            refused(format!(
                "the check of revision {} records no workload counting version",
                revision.unwrap_or("unknown")
            ))
        })
}

/// The ledger head: its index row and run header.
struct Head {
    row: ledger::IndexRow,
    header: Header,
}

/// The ledger's head entry, if any.
fn head(ledger_directory: &Path) -> Result<Option<Head>> {
    let index = ledger::index_text(ledger_directory).map_err(refused)?;
    let rows = ledger::read_index(index.as_deref()).map_err(refused)?;
    let Some(row) = rows.last().cloned() else {
        return Ok(None);
    };
    let header = read_header(&ledger_directory.join(&row.entry))?;
    if header.golden_sha256 != row.golden_sha256 {
        return Err(refused(format!(
            "head {}: run.json and index.tsv record different golden outputs",
            row.entry
        )));
    }
    Ok(Some(Head { row, header }))
}

impl Head {
    /// The head's golden output over `samples`, read once a check is known
    /// to be comparable with the head.
    fn golden(&self, ledger_directory: &Path, samples: &[String]) -> Result<Golden> {
        let path = ledger_directory.join(&self.row.entry).join(GOLDEN);
        let golden = Golden::parse(&text(&path, "head golden output")?, samples)
            .map_err(|message| refused(format!("head {}: {message}", self.row.entry)))?;
        let sha256 = golden.sha256();
        if sha256 != self.row.golden_sha256 {
            return Err(refused(format!(
                "head {}: golden.tsv hashes to {sha256}, not what index.tsv and run.json record",
                self.row.entry
            )));
        }
        Ok(golden)
    }
}

/// Reads a UTF-8 ledger file.
fn text(path: &Path, kind: &'static str) -> Result<String> {
    String::from_utf8(read(path, kind, MAX_LEDGER_BYTES)?)
        .map_err(|_| refused(format!("{} is not UTF-8", path.display())))
}

/// Reads the run header of the check in `directory`.
fn read_header(directory: &Path) -> Result<Header> {
    let path = directory.join(RUN_JSON);
    serde_json::from_str(&text(&path, "regress run header")?)
        .map_err(|error| refused(format!("{}: {error}", path.display())))
}

/// Refuses a check that measures something else than the head: another
/// series, manifest, profile, reference, review window, zones, or run
/// contexts; another truth or re-check list unless `rebase`; and, with
/// `rebase`, another configuration, or no change of truth or re-check list
/// at all.
fn comparable((entry, head): (&str, &Header), candidate: &Header, rebase: bool) -> Result<()> {
    if head.series != candidate.series {
        return Err(refused(format!(
            "the ledger's head {entry} is series {}, the manifest series {}: another series \
             takes its own ledger directory, which its first check seeds",
            head.series, candidate.series
        )));
    }
    let recheck = |header: &Header| {
        header
            .recheck
            .as_ref()
            .map(|recheck| recheck.sha256.clone())
    };
    let contexts = |header: &Header| {
        header
            .run_context
            .as_ref()
            .map(|context| context.sha256.clone())
    };
    let truth_moved =
        head.truth.digest != candidate.truth.digest || recheck(head) != recheck(candidate);
    let mut differ = [
        (
            "manifest body",
            head.manifest.body_sha256 != candidate.manifest.body_sha256,
        ),
        ("profile", head.profile_sha256 != candidate.profile_sha256),
        (
            "reference",
            head.reference_sha256 != candidate.reference_sha256,
        ),
        (
            "review window",
            head.review_window != candidate.review_window,
        ),
        ("zone widths", head.zones != candidate.zones),
        ("run contexts", contexts(head) != contexts(candidate)),
        (
            "truth digest (use --rebase)",
            !rebase && head.truth.digest != candidate.truth.digest,
        ),
        (
            "re-check list (use --rebase)",
            !rebase && recheck(head) != recheck(candidate),
        ),
        (
            "configuration (a rebase changes the truth only)",
            rebase
                && head.version_key.configuration_sha256
                    != candidate.version_key.configuration_sha256,
        ),
        (
            "nothing: --rebase needs another truth digest or re-check list",
            rebase && !truth_moved,
        ),
    ]
    .into_iter()
    .filter_map(|(name, differs)| differs.then_some(name))
    .peekable();
    if differ.peek().is_none() {
        return Ok(());
    }
    Err(refused(format!(
        "not comparable with head {entry}: {} differ",
        differ.collect::<Vec<_>>().join(", ")
    )))
}

/// Refuses an intent file that names a sample of the manifest.
fn names_no_sample(text: &str, samples: &BTreeSet<&str>) -> Result<()> {
    if text
        .split(|symbol: char| !symbol.is_ascii_alphanumeric() && symbol != '_')
        .any(|token| samples.contains(token))
    {
        return Err(refused("the intent file names a sample of the manifest"));
    }
    Ok(())
}

/// Reads the intent file at `path` against the zones and regions its
/// `where` may name, and refuses one that names a sample of `samples`.
fn read_intents(
    path: &Path,
    zones: &Zones,
    regions: &[Span],
    samples: &[String],
) -> Result<Intents> {
    let content = String::from_utf8(read(path, "intent file", MAX_INTENT_BYTES)?)
        .map_err(|_| refused(format!("{} is not UTF-8", path.display())))?;
    let labels = zones
        .instances()
        .into_iter()
        .map(|instance| zones.label(instance))
        .collect();
    let regions = regions.iter().map(|region| region.name.clone()).collect();
    let parsed = intents::parse(&content, &labels, &regions)
        .map_err(|message| refused(format!("intent file {}: {message}", path.display())))?;
    names_no_sample(&content, &samples.iter().map(String::as_str).collect())?;
    Ok(parsed)
}

/// What a check derives from its golden output, the head's, and the
/// intents: the judgement and the fields of `run.json` that record it.
struct Judged {
    judgement: Option<Judgement>,
    verdict: Verdict,
    rules: Vec<String>,
    golden_sha256: String,
    kinds: Vec<KindDigest>,
    strata: Vec<StratumDigest>,
    counts: Counts,
    changes: ChangeCounts,
}

impl Judged {
    /// Judges `golden` against `head` (none for a series' first check),
    /// with each sample in its manifest stratum.
    fn new(
        golden: &Golden,
        head: Option<&Golden>,
        intents: Option<&Intents>,
        context: &Context<'_>,
        strata: &[&str],
    ) -> Self {
        let judgement = head.map(|head| judge::judge(head, golden, intents, context));
        let changed = judgement
            .as_ref()
            .map(Judgement::changed_samples)
            .unwrap_or_default();
        Self {
            verdict: judgement
                .as_ref()
                .map_or(Verdict::Initial, |judgement| judgement.verdict),
            rules: judgement
                .iter()
                .flat_map(|judgement| &judgement.rules)
                .map(|rule| rule.label().to_owned())
                .collect(),
            golden_sha256: golden.sha256(),
            kinds: golden
                .kind_digests()
                .into_iter()
                .map(|(kind, sha256, lines)| KindDigest {
                    kind: kind.label().to_owned(),
                    lines,
                    sha256,
                })
                .collect(),
            strata: golden
                .stratum_digests(strata)
                .into_iter()
                .map(|(stratum, members, golden_sha256)| StratumDigest {
                    changed_samples: changed
                        .iter()
                        .filter(|sample| strata[**sample] == stratum)
                        .count(),
                    stratum,
                    samples: members,
                    golden_sha256,
                })
                .collect(),
            counts: golden.counts(),
            changes: judgement
                .as_ref()
                .map(|judgement| ChangeCounts {
                    changed_samples: changed.len(),
                    changed_lines: judgement.changes.len(),
                    explained_lines: judgement.explained_lines(),
                    posthoc_rows: intents.map_or(0, Intents::posthoc_rows),
                })
                .unwrap_or_default(),
            judgement,
        }
    }

    /// Writes the derived fields into `header`.
    fn record_in(&self, header: &mut Header) {
        self.verdict.label().clone_into(&mut header.verdict);
        header.rules.clone_from(&self.rules);
        header.golden_sha256.clone_from(&self.golden_sha256);
        header.kinds.clone_from(&self.kinds);
        header.strata.clone_from(&self.strata);
        header.counts = self.counts;
        header.changes = self.changes;
    }

    /// The fields of `header` that differ from the derived ones.
    fn disagreements(&self, header: &Header) -> Vec<&'static str> {
        [
            ("verdict", header.verdict != self.verdict.label()),
            ("rules", header.rules != self.rules),
            ("golden_sha256", header.golden_sha256 != self.golden_sha256),
            ("kinds", header.kinds != self.kinds),
            ("strata", header.strata != self.strata),
            ("counts", header.counts != self.counts),
            ("changes", header.changes != self.changes),
        ]
        .into_iter()
        .filter_map(|(name, differs)| differs.then_some(name))
        .collect()
    }
}

/// Runs `regress check`: returns `regress.md` and the verdict.
fn check(args: &CheckArgs, git: &dyn Git, identity: &Identity) -> Result<(String, Verdict)> {
    let started = Instant::now();
    let mut loaded = Loaded::read(
        (&args.manifest, regress_panel),
        (&args.truth, Some(&args.recheck)),
        (&args.dna, &args.review),
        identity,
    )?;
    let version = loaded
        .manifest
        .header
        .get("version")
        .ok_or_else(|| refused("the manifest declares no version"))?;
    let series = format!("regress-{version}");
    let samples = loaded
        .manifest
        .samples
        .iter()
        .map(|sample| sample.id.clone())
        .collect::<Vec<_>>();
    let strata = loaded
        .manifest
        .samples
        .iter()
        .map(|sample| sample.stratum.as_str())
        .collect::<Vec<_>>();
    let head = head(&args.ledger)?;
    if head.is_none() && (args.rebase || args.intent.is_some()) {
        return Err(refused(
            "the ledger has no head: the first check takes no --rebase and no --intent",
        ));
    }
    let zones = loaded.zones();
    let source = crate::source(identity);
    let intents = args
        .intent
        .as_deref()
        .map(|path| -> Result<(Intents, IntentHeader)> {
            let parsed = read_intents(path, &zones, loaded.geometry.regions(), &samples)?;
            let base = head.as_ref().map_or("-", |head| head.row.entry.as_str());
            if parsed.series != series || parsed.base != base {
                return Err(refused(format!(
                    "intent file {}: declared for series {} base {}, not {series} base {base}",
                    path.display(),
                    parsed.series,
                    parsed.base
                )));
            }
            let revision = source.revision.as_deref().ok_or_else(|| {
                refused("the harness records no revision to check the intent file against")
            })?;
            let committed = git.committed(path, revision).map_err(refused)?;
            let blob = git.blob(path).map_err(refused)?;
            if committed != blob {
                return Err(refused(format!(
                    "intent file {} differs from its commit at {revision}",
                    path.display()
                )));
            }
            let header = IntentHeader {
                path: absolute(path),
                blob,
                pr: parsed.pr.clone(),
                rows: parsed.rows.len(),
                posthoc_rows: parsed.posthoc_rows(),
            };
            Ok((parsed, header))
        })
        .transpose()?;

    let truth_check = TruthCheck::new(&loaded.manifest, &loaded.truth.database);
    // A series starts on the truth its manifest pins; later checks measure
    // against the head's truth (`comparable`), which only a rebase changes.
    if head.is_none() && truth_check.state() != TruthState::Unchanged {
        return Err(refused(format!(
            "the truth is not the one the manifest pins (truth {}, {} sample(s) changed): \
             the first check of a series needs the pinned database",
            truth_check.state().label(),
            truth_check.changed.len()
        )));
    }
    let configuration_sha256 =
        crate::runner::sha256_file(&args.dna.config).map_err(|source| Error::Read {
            kind: "configuration",
            path: args.dna.config.clone(),
            source,
        })?;
    let mut header = Header {
        harness: identity.name.to_owned(),
        harness_version: identity.version.to_owned(),
        source,
        threads: args.dna.threads,
        series,
        manifest: ManifestHeader {
            path: absolute(&args.manifest),
            fields: loaded.manifest.header.clone(),
            body_sha256: loaded.manifest.body_sha256.clone(),
            rows: loaded.manifest.samples.len(),
        },
        truth: TruthIdentity {
            path: absolute(&args.truth),
            sha256: loaded.truth.sha256.clone(),
            digest: truth_check.current.clone(),
            recorded_digest: truth_check.recorded.clone(),
            state: truth_check.state().label().to_owned(),
            changed_samples: truth_check.changed.len(),
        },
        recheck: loaded.recheck.as_ref().map(|(_, header)| header.clone()),
        configuration_path: absolute(&args.dna.config),
        reference_path: absolute(&args.dna.reference),
        profile_sha256: loaded.profile.identity().sha256.clone(),
        reference_sha256: loaded.reference.sequence_sha256.clone(),
        review_window: args.review.review_window,
        zones: args.review.header(),
        run_context: None,
        version_key: VersionKey::current(&configuration_sha256, identity),
        workload_version: Some(workload::COUNTING_VERSION),
        evidence_cache: None,
        placement_cache: None,
        head: head.as_ref().map(|head| head.row.entry.clone()),
        intent: intents.as_ref().map(|(_, header)| header.clone()),
        rebase: args.rebase,
        verdict: String::new(),
        rules: Vec::new(),
        warnings: Vec::new(),
        golden_sha256: String::new(),
        kinds: Vec::new(),
        strata: Vec::new(),
        counts: Counts::default(),
        changes: ChangeCounts::default(),
        seconds: Seconds {
            total: 0.0,
            truth: loaded.truth.seconds,
            dna: 0.0,
        },
    };
    let (_, contexts_sha256) = if args.dna.run_context.is_some() {
        args.dna.contexts(&loaded.manifest.samples)?
    } else {
        Default::default()
    };
    header.run_context = args
        .dna
        .run_context
        .as_deref()
        .map(|directory| RunContextHeader {
            directory: absolute(directory),
            batches: 0,
            sha256: contexts_sha256,
        });
    // Everything the comparison rests on is known before DNA runs.
    if let Some(head) = &head {
        comparable((&head.row.entry, &head.header), &header, args.rebase)?;
    }
    let head_golden = head
        .as_ref()
        .map(|head| head.golden(&args.ledger, &samples))
        .transpose()?;
    prepare_output(&args.output)?;

    let ran = Ran::run(&loaded, &args.dna, None, identity)?;
    if let Some(other) = ran
        .outcomes
        .iter()
        .filter_map(|outcome| outcome.result.as_ref().ok())
        .map(|called| &called.identity.configuration_sha256)
        .find(|sha256| **sha256 != configuration_sha256)
    {
        return Err(refused(format!(
            "DNA analysed with configuration {other}, not {configuration_sha256}"
        )));
    }
    header.run_context = args.dna.context_header(&ran);
    header.evidence_cache = ran.cache.as_ref().map(|cache| cache.evidence.usage());
    header.placement_cache = ran.cache.as_ref().map(|cache| cache.placed.usage());
    header.seconds.dna = ran.seconds;

    let (recheck, _) = loaded.recheck.take().unzip();
    let tracts = loaded.tracts();
    let windows = workload::windows(args.review.review_window, loaded.geometry.regions());
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
    let report = Report::new(
        &rows,
        (geometry, &truth_check),
        (forms, reference.sequence.as_str()),
    )
    .with_review(
        args.review.review_window,
        None,
        Some(zones.clone()),
        recheck,
        tracts,
    );
    let golden = golden::build(&Sources {
        samples: &manifest.samples,
        outcomes: &ran.outcomes,
        report: &report,
        windows: &windows,
        zones: &zones,
    });

    let head_counting = head
        .as_ref()
        .map_or(Ok(workload::COUNTING_VERSION), |head| {
            counting_version(&head.header)
        })?;
    let judged = Judged::new(
        &golden,
        head_golden.as_ref(),
        intents.as_ref().map(|(intents, _)| intents),
        &Context {
            zones: &zones,
            regions: geometry.regions(),
            rebase: args.rebase,
            version_moved: head
                .as_ref()
                .is_some_and(|head| head.header.version_key != header.version_key),
            workload: (head_counting, workload::COUNTING_VERSION),
        },
        &strata,
    );
    judged.record_in(&mut header);
    let (judgement, verdict) = (judged.judgement.as_ref(), judged.verdict);
    header.warnings = warnings(
        head.as_ref(),
        &header,
        judgement,
        intents.as_ref(),
        &golden.samples,
    );
    header.seconds.total = started.elapsed().as_secs_f64();

    let intents = intents.as_ref().map(|(intents, _)| intents);
    let report_text = summary(
        &header,
        head.as_ref().map(|head| &head.header),
        judgement,
        (&golden.samples, &strata),
        intents,
    );
    write(&args.output, GOLDEN, &golden.text())?;
    write(
        &args.output,
        CHANGES,
        &judgement.map_or_else(
            || judge::HEADER.to_owned(),
            |judgement| judge::table(judgement, &golden.samples, &strata, intents),
        ),
    )?;
    write(&args.output, REPORT, &report_text)?;
    write(&args.output, RUN_JSON, &header.text()?)?;
    Ok((report_text, verdict))
}

/// What a check reports beside the verdict: a moved version key without a
/// change, a moved workload counting version, intent rows that matched
/// nothing, rows added after a check, the
/// version-key fields that moved, the configuration that changed, escapes
/// a rebase measures, and escape lines that changed while their sample's
/// escape count in the zone did not rise.
fn warnings(
    head: Option<&Head>,
    header: &Header,
    judgement: Option<&Judgement>,
    intents: Option<&(Intents, IntentHeader)>,
    samples: &[String],
) -> Vec<String> {
    let mut warnings = Vec::new();
    let (Some(head), Some(judgement)) = (head, judgement) else {
        return warnings;
    };
    let moved = header.version_key.moved(&head.header.version_key);
    if !moved.is_empty() {
        warnings.push(format!("version key moved: {}", moved.join(", ")));
        if judgement.changes.is_empty() {
            warnings.push(
                "the version key moved but no line changed: the change may lie outside the set"
                    .to_owned(),
            );
        }
    }
    if let (Ok(before), Ok(after)) = (counting_version(&head.header), counting_version(header))
        && before != after
    {
        warnings.push(format!(
            "workload counting version moved: {before} -> {after}"
        ));
    }
    if let Some((intents, _)) = intents {
        for (index, (intent, (lines, _))) in intents.rows.iter().zip(&judgement.intents).enumerate()
        {
            if *lines == 0 && judgement.workload_intent != Some(index) {
                warnings.push(format!("unused intent {}", intent.id));
            }
            if intent.posthoc {
                warnings.push(format!(
                    "intent {} was declared after a check (posthoc)",
                    intent.id
                ));
            }
        }
    }
    if judgement.rebase_escapes > 0 {
        warnings.push(format!(
            "the rebase turns {} block(s) into non-het escapes",
            judgement.rebase_escapes
        ));
    }
    for change in judgement
        .reshaped_escapes
        .iter()
        .filter_map(|index| judgement.changes.get(*index))
    {
        warnings.push(format!(
            "escape changed, count not risen: {} block {} {} {} -> {}",
            samples[change.sample],
            change.key,
            change.change.label(),
            change.before.as_deref().unwrap_or("-"),
            change.after.as_deref().unwrap_or("-"),
        ));
    }
    warnings
}

/// A header's revision, marked when the build was dirty or its state unknown.
fn revision(header: &Header) -> String {
    format!(
        "{}{}",
        header.source.revision.as_deref().unwrap_or("unknown"),
        match header.source.dirty {
            Some(false) => "",
            Some(true) => " (dirty)",
            None => " (dirty state unknown)",
        }
    )
}

/// `regress.md`: the verdict, the rules that fired, the version keys, the
/// changes by kind, stratum, and intent, the first unexplained lines (sample
/// lines: the file stays local), and the counts.
fn summary(
    header: &Header,
    head: Option<&Header>,
    judgement: Option<&Judgement>,
    (samples, strata): (&[String], &[&str]),
    intents: Option<&Intents>,
) -> String {
    let mut text = format!("# Regress check: {}\n\n", header.series);
    let _ = writeln!(text, "**Verdict: {}**\n", header.verdict);
    if !header.rules.is_empty() {
        let _ = writeln!(text, "Rules fired: {}.\n", header.rules.join(", "));
    }
    for warning in &header.warnings {
        let _ = writeln!(text, "- Warning: {warning}");
    }
    if !header.warnings.is_empty() {
        text.push('\n');
    }
    let _ = writeln!(
        text,
        "| | Head {} | Candidate |\n|---|---|---|",
        header.head.as_deref().unwrap_or("(none)")
    );
    let key = |header: Option<&Header>, field: fn(&Header) -> String| {
        header.map_or_else(|| "-".to_owned(), field)
    };
    for (name, field) in [
        ("revision", revision as fn(&Header) -> String),
        ("dna_version", |header| {
            header.version_key.dna_version.clone()
        }),
        ("plugins", |header| header.version_key.plugins.clone()),
        ("contracts", |header| header.version_key.contracts.clone()),
        ("configuration sha256", |header| {
            header.version_key.configuration_sha256.clone()
        }),
        ("workload counting version", |header| {
            counting_version(header).map_or_else(|_| "-".to_owned(), |version| version.to_string())
        }),
        ("truth digest", |header| header.truth.digest.clone()),
        ("re-check list sha256", |header| {
            header
                .recheck
                .as_ref()
                .map_or_else(|| "-".to_owned(), |recheck| recheck.sha256.clone())
        }),
        ("golden sha256", |header| header.golden_sha256.clone()),
    ] {
        let _ = writeln!(
            text,
            "| {name} | {} | {} |",
            key(head, field),
            field(header)
        );
    }
    let empty = Judgement {
        changes: Vec::new(),
        rules: Vec::new(),
        verdict: Verdict::Initial,
        intents: Vec::new(),
        rebase_escapes: 0,
        reshaped_escapes: Vec::new(),
        workload_intent: None,
    };
    let judgement = judgement.unwrap_or(&empty);

    text.push_str(
        "\n## Changes by kind\n\n| Kind | Added | Removed | Changed |\n|---|---:|---:|---:|\n",
    );
    for kind in Kind::ALL {
        let count = |change: Change| {
            judgement
                .changes
                .iter()
                .filter(|line| line.kind == kind && line.change == change)
                .count()
        };
        let _ = writeln!(
            text,
            "| {} | {} | {} | {} |",
            kind.label(),
            count(Change::Added),
            count(Change::Removed),
            count(Change::Changed)
        );
    }

    text.push_str("\n## By stratum\n\n| Stratum | Changed samples | Members |\n|---|---:|---:|\n");
    for stratum in &header.strata {
        let _ = writeln!(
            text,
            "| {} | {} | {} |",
            stratum.stratum, stratum.changed_samples, stratum.samples
        );
    }

    if let Some(intents) = intents {
        text.push_str(
            "\n## By intent\n\n| Intent | Lines | Samples | Max samples | Posthoc |\n\
             |---|---:|---:|---:|---|\n",
        );
        for (intent, (lines, count)) in intents.rows.iter().zip(&judgement.intents) {
            let _ = writeln!(
                text,
                "| {} | {lines} | {count} | {} | {} |",
                intent.id,
                intent
                    .max_samples
                    .map_or_else(|| "-".to_owned(), |maximum| maximum.to_string()),
                if intent.posthoc { "yes" } else { "no" }
            );
        }
    }

    let unexplained = judgement
        .changes
        .iter()
        .filter(|change| !change.explained())
        .collect::<Vec<_>>();
    if !unexplained.is_empty() {
        let _ = writeln!(
            text,
            "\n## Unexplained or forbidden lines ({}, first {LISTED_LINES})\n\n\
             These lines name samples: this file stays local.\n\n\
             | Sample | Stratum | Kind | Key | Change | Before | After | Label |\n\
             |---|---|---|---|---|---|---|---|",
            unexplained.len()
        );
        for change in unexplained.into_iter().take(LISTED_LINES) {
            let label = match change.label {
                Label::Forbidden(rule) => rule.label(),
                _ => "unexplained",
            };
            let _ = writeln!(
                text,
                "| {} | {} | {} | {} | {} | {} | {} | {label} |",
                samples[change.sample],
                strata[change.sample],
                change.kind.label(),
                change.key,
                change.change.label(),
                change.before.as_deref().unwrap_or("-"),
                change.after.as_deref().unwrap_or("-"),
            );
        }
    }

    let counts = &header.counts;
    let changes = &header.changes;
    let _ = write!(
        text,
        "\n## Summary\n\n| Measure | Count |\n|---|---:|\n\
         | Samples | {} |\n| Failed | {} |\n| Clear | {} |\n| With a check | {} |\n\
         | Check items | {} |\n| Noted tract runs | {} |\n| Alleles | {} |\n\
         | Non-het blocks | {} |\n| Non-het blocks in an unchecked zone (escapes) | {} |\n\
         | Het blocks | {} |\n| Length-note blocks | {} |\n| Truth-N blocks | {} |\n\
         | Re-check blocks | {} |\n| Changed samples | {} |\n| Changed lines | {} |\n\
         | Explained lines | {} |\n| Posthoc intent rows | {} |\n",
        counts.samples,
        counts.failed,
        counts.clear,
        counts.with_check,
        counts.checks,
        counts.notes,
        counts.alleles,
        counts.blocks_non_het,
        counts.escapes,
        counts.blocks_het,
        counts.blocks_length_note,
        counts.blocks_truth_n,
        counts.blocks_recheck,
        changes.changed_samples,
        changes.changed_lines,
        changes.explained_lines,
        changes.posthoc_rows,
    );
    text
}

/// The ledger entry kind of a recorded verdict.
const fn entry_kind(verdict: Verdict) -> Option<&'static str> {
    match verdict {
        Verdict::Initial => Some("initial"),
        Verdict::Same | Verdict::Explained => Some("change"),
        Verdict::Rebase => Some("rebase"),
        Verdict::Fail => None,
    }
}

/// Runs `regress record` on `date`: returns what it recorded.
///
/// It refuses before writing anything: a dirty build, a build of another
/// commit than [`MAIN`]'s current one, a failed verdict, a moved head, an
/// existing entry, baselines that do not end at the head, and a `run.json`
/// that a re-judgement of its golden output against the head and the intent
/// file does not reproduce.
fn record(args: &RecordArgs, git: &dyn Git, date: &str) -> Result<String> {
    if !ledger::valid_entry(&args.entry) {
        return Err(refused(format!(
            "entry {:?} is not r and four or more digits",
            args.entry
        )));
    }
    let header = read_header(&args.run)?;
    let golden_text = text(&args.run.join(GOLDEN), "golden output")?;
    if golden::sha256(&golden_text) != header.golden_sha256 {
        return Err(refused(
            "the run's golden.tsv does not hash to the golden_sha256 its run.json records",
        ));
    }
    let revision = match (&header.source.revision, header.source.dirty) {
        (Some(revision), Some(false)) => revision.clone(),
        _ => {
            return Err(refused(
                "the run's build was dirty or records no revision: record a clean build of main",
            ));
        }
    };
    let verdict = Verdict::parse(&header.verdict)
        .ok_or_else(|| refused(format!("unknown verdict {:?}", header.verdict)))?;
    let kind = entry_kind(verdict).ok_or_else(|| refused("the run's verdict is REGRESS FAIL"))?;
    let index = ledger::index_text(&args.ledger).map_err(refused)?;
    let rows = ledger::read_index(index.as_deref()).map_err(refused)?;
    let current = rows.last().map(|row| row.entry.clone());
    if rows.iter().any(|row| row.entry == args.entry) || args.ledger.join(&args.entry).exists() {
        return Err(refused(format!("entry {} already exists", args.entry)));
    }
    if header.head != current {
        return Err(refused(format!(
            "the head moved: the run compared with {}, the ledger head is {}",
            header.head.as_deref().unwrap_or("none"),
            current.as_deref().unwrap_or("none")
        )));
    }
    let baseline = text(&args.baseline, "ledger baseline")?;
    let strata = text(&args.strata, "strata baseline")?;
    let entries = ledger::baseline_entries(
        &baseline,
        ("regress-ledger", &header.series),
        ledger::LEDGER_COLUMNS,
    )
    .map_err(|message| refused(format!("{}: {message}", args.baseline.display())))?;
    let stratum_entries = ledger::baseline_entries(
        &strata,
        ("regress-strata", &header.series),
        ledger::STRATA_COLUMNS,
    )
    .map_err(|message| refused(format!("{}: {message}", args.strata.display())))?;
    if entries.contains(&args.entry) || stratum_entries.contains(&args.entry) {
        return Err(refused(format!(
            "entry {} already exists in the baselines",
            args.entry
        )));
    }
    for (path, last) in [
        (&args.baseline, entries.last()),
        (&args.strata, stratum_entries.last()),
    ] {
        if last != current.as_ref() {
            return Err(refused(format!(
                "{} ends at {}, the local ledger at {}",
                path.display(),
                last.map_or("none", String::as_str),
                current.as_deref().unwrap_or("none")
            )));
        }
    }
    let repository = args.baseline.parent().unwrap_or(Path::new("."));
    let main = git.main_revision(repository).map_err(refused)?;
    if revision != main {
        return Err(refused(format!(
            "the run's revision {revision} is not {MAIN} ({main}): fetch, then check and record \
             the current main"
        )));
    }
    if let Some(intent) = &header.intent
        && git.blob(&intent.path).map_err(refused)? != intent.blob
    {
        return Err(refused(format!(
            "the intent file {} changed since the check",
            intent.path.display()
        )));
    }
    let head = head(&args.ledger)?;
    rejudge(&header, &golden_text, head.as_ref(), &args.ledger)?;

    let row = ledger::IndexRow {
        entry: args.entry.clone(),
        kind: kind.to_owned(),
        base: header.head.clone().unwrap_or_else(|| "-".to_owned()),
        revision,
        golden_sha256: header.golden_sha256.clone(),
        verdict: header.verdict.clone(),
        recorded: date.to_owned(),
    };
    let files = [
        (
            args.baseline.clone(),
            format!(
                "{baseline}{}",
                ledger::ledger_row(&args.entry, (date, kind), &header)
            ),
        ),
        (
            args.strata.clone(),
            format!("{strata}{}", ledger::strata_rows(&args.entry, &header)),
        ),
        // The index goes last: it names the head.
        (
            args.ledger.join(ledger::INDEX),
            format!(
                "{}{}",
                index.unwrap_or_else(|| ledger::INDEX_HEADER.to_owned()),
                row.line()
            ),
        ),
    ];
    install(args, &header, &files)?;
    Ok(format!(
        "dna-eval: recorded {} ({kind}, {}) on {date}: {} and {} have a new row; commit them{}",
        args.entry,
        header.verdict,
        args.baseline.display(),
        args.strata.display(),
        header
            .intent
            .as_ref()
            .map_or_else(String::new, |intent| format!(
                " beside the intent file {}",
                intent.path.display()
            ))
    ))
}

/// Re-derives a checked run's verdict and recorded fields from its golden
/// output, the head's, the manifest, the configuration's geometry, and the
/// intent file, and refuses a run whose `run.json` records anything else or
/// whose inputs changed since the check.
fn rejudge(
    header: &Header,
    golden_text: &str,
    head: Option<&Head>,
    ledger_directory: &Path,
) -> Result<()> {
    let manifest = manifest::parse(&read(
        &header.manifest.path,
        "manifest",
        MAX_MANIFEST_BYTES,
    )?)?;
    regress_panel(&manifest)?;
    if manifest.body_sha256 != header.manifest.body_sha256 {
        return Err(refused(format!(
            "the manifest {} changed since the check",
            header.manifest.path.display()
        )));
    }
    let configuration_sha256 =
        crate::runner::sha256_file(&header.configuration_path).map_err(|source| Error::Read {
            kind: "configuration",
            path: header.configuration_path.clone(),
            source,
        })?;
    let profile = crate::target_profile(&header.configuration_path)?;
    if configuration_sha256 != header.version_key.configuration_sha256
        || profile.identity().sha256 != header.profile_sha256
    {
        return Err(refused(format!(
            "the configuration {} changed since the check",
            header.configuration_path.display()
        )));
    }
    let geometry = Geometry::of(&profile);
    let zones = Zones::new(
        geometry.regions(),
        geometry.repeats(),
        Widths {
            after: header.zones.after,
            edge: header.zones.edge,
            overrun: header.zones.overrun,
        },
    );
    let samples = manifest
        .samples
        .iter()
        .map(|sample| sample.id.clone())
        .collect::<Vec<_>>();
    let strata = manifest
        .samples
        .iter()
        .map(|sample| sample.stratum.as_str())
        .collect::<Vec<_>>();
    if let Some(head) = head {
        comparable((&head.row.entry, &head.header), header, header.rebase)?;
    } else if header.truth.state != TruthState::Unchanged.label() {
        return Err(refused(
            "the series' first check ran on another truth than the manifest pins",
        ));
    }
    let head_golden = head
        .map(|head| head.golden(ledger_directory, &samples))
        .transpose()?;
    let intents = header
        .intent
        .as_ref()
        .map(|intent| read_intents(&intent.path, &zones, geometry.regions(), &samples))
        .transpose()?;
    let golden = Golden::parse(golden_text, &samples)
        .map_err(|message| refused(format!("the run's golden.tsv: {message}")))?;
    let counting = counting_version(header)?;
    let head_counting = head.map_or(Ok(counting), |head| counting_version(&head.header))?;
    let judged = Judged::new(
        &golden,
        head_golden.as_ref(),
        intents.as_ref(),
        &Context {
            zones: &zones,
            regions: geometry.regions(),
            rebase: header.rebase,
            version_moved: head.is_some_and(|head| head.header.version_key != header.version_key),
            workload: (head_counting, counting),
        },
        &strata,
    );
    let disagreements = judged.disagreements(header);
    if !disagreements.is_empty() {
        return Err(refused(format!(
            "the run's run.json disagrees with a re-judgement of its golden.tsv: {}",
            disagreements.join(", ")
        )));
    }
    Ok(())
}

/// The staging path of `path`: a hidden sibling.
fn staged(path: &Path) -> PathBuf {
    let name = path
        .file_name()
        .map_or_else(String::new, |name| name.to_string_lossy().into_owned());
    path.with_file_name(format!(".{name}.partial"))
}

/// Writes a recorded entry: stages the entry directory and every file
/// beside its target, then renames each into place, `files` in order, so
/// an interrupted record never leaves a half-written file and leaves the
/// head unchanged unless the index moved.
fn install(args: &RecordArgs, header: &Header, files: &[(PathBuf, String)]) -> Result<()> {
    let failure = |path: &Path| {
        let path = path.to_path_buf();
        move |source| Error::Write { path, source }
    };
    let directory = args.ledger.join(&args.entry);
    let staging = staged(&directory);
    for path in std::iter::once(&staging).chain(files.iter().map(|(path, _)| path)) {
        let partial = if *path == staging {
            staging.clone()
        } else {
            staged(path)
        };
        if partial.exists() {
            return Err(refused(format!(
                "{} is left from an interrupted record: remove it",
                partial.display()
            )));
        }
    }
    let stage = || -> Result<()> {
        fs::create_dir_all(&args.ledger).map_err(failure(&args.ledger))?;
        fs::create_dir(&staging).map_err(failure(&staging))?;
        for name in [GOLDEN, CHANGES, REPORT, RUN_JSON] {
            fs::copy(args.run.join(name), staging.join(name))
                .map_err(failure(&staging.join(name)))?;
        }
        if let Some(intent) = &header.intent {
            fs::copy(&intent.path, staging.join(INTENT)).map_err(failure(&staging.join(INTENT)))?;
        }
        for (path, content) in files {
            fs::write(staged(path), content).map_err(failure(&staged(path)))?;
        }
        Ok(())
    };
    if let Err(error) = stage() {
        // Best effort: the staged copies are not the ledger's.
        let _ = fs::remove_dir_all(&staging);
        for (path, _) in files {
            let _ = fs::remove_file(staged(path));
        }
        return Err(error);
    }
    fs::rename(&staging, &directory).map_err(failure(&directory))?;
    for (path, _) in files {
        fs::rename(staged(path), path).map_err(failure(path))?;
    }
    Ok(())
}

#[cfg(test)]
mod tests;
