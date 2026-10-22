//! The run directory of `dna-eval run`, `dna.eval_run/v1`: its file names,
//! the run header `run.json` as one typed record that `run` writes and
//! `diff` reads back, the tables `run` writes beside the report's, and the
//! strict reader `diff` reads a run's tables with.
//!
//! The reader fails on a table without a column it needs, with a column
//! named twice, or with a row of another width than its header, so a run of
//! another contract revision is refused rather than read as empty cells.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::fs;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::cache;
use crate::components;
use crate::conventions;
use crate::error::{DiffError, Error, Result};
use crate::estimate;
use crate::manifest::{Month, Sample};
use crate::runner::{Contexts, Identity, Outcome};
use crate::tracts;
use crate::truth::TruthSample;
use crate::workload::Passes;

/// The run header.
pub(crate) const RUN_JSON: &str = "run.json";
/// Per sample: the comparisons, the class and cause, and the review workload.
pub(crate) const SCORES: &str = "scores.tsv";
/// Per unreviewed difference block: its class and zone.
pub(crate) const UNREVIEWED: &str = "unreviewed.tsv";
/// Per check item of `dna review`.
pub(crate) const CHECKS: &str = "checks.tsv";
/// Per noted tract run of `dna review`.
pub(crate) const NOTES: &str = "notes.tsv";
/// Per run-signal group, with `--run-signal`.
pub(crate) const SIGNALS: &str = "signals.tsv";
/// The run-context coverage, with `--run-context`.
pub(crate) const RUN_CONTEXT: &str = "run_context.tsv";
/// Per phase-mixture component.
pub(crate) const COMPONENTS: &str = "components.tsv";
/// Per tract decode.
pub(crate) const TRACTS: &str = "tracts.tsv";
/// How the reviewers write each window sequence.
pub(crate) const FORMS: &str = "forms.tsv";

/// The header of `signals.tsv`.
pub(crate) const SIGNALS_HEADER: &str =
    "group\treads\tcalls\tkeys\tkeys_20_calls\tkeys_40_calls\tkeys_100_calls\n";
/// The header of `run_context.tsv`.
pub(crate) const RUN_CONTEXT_HEADER: &str = "month\tsamples\tanalysed\twith_context\trun1_calibrated\trun2_calibrated\tbatches\trun1_pairs\trun2_pairs\n";
/// The header of `checks.tsv`.
pub(crate) const CHECKS_HEADER: &str = "sample_id\tfirst\tlast\tlabels\n";
/// The header of `notes.tsv`.
pub(crate) const NOTES_HEADER: &str =
    "sample_id\tfirst\tlast\ttract\trun\tstatus\tcalled\treading\tconvention\n";

/// `run.json`: what was run, on what, and how long it took.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub(crate) struct RunHeader {
    pub(crate) harness: String,
    pub(crate) harness_version: String,
    pub(crate) source: Source,
    pub(crate) threads: usize,
    pub(crate) manifest: ManifestHeader,
    pub(crate) truth: TruthHeader,
    pub(crate) configuration_path: PathBuf,
    pub(crate) reference_path: PathBuf,
    /// Distinct identities DNA reported; one unless the inputs changed during
    /// the run.
    pub(crate) identities: Vec<Identity>,
    /// Window of the review workload, in positions.
    pub(crate) review_window: i64,
    /// First batch month of the current data, when given.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) current_since: Option<Month>,
    /// The run signal grouping, when the reviews were judged against one.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) run_signal: Option<SignalGrouping>,
    /// The fewest reads of a batch for its own run signal, with
    /// `run_signal` `batch`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) run_signal_batch_reads: Option<usize>,
    /// The forced-window list the workload also counts as reviewed.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) forced_windows: Option<PathBuf>,
    /// The re-check list the pass measures read.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) recheck: Option<RecheckHeader>,
    /// The read evidence cache the run read and filled.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) evidence_cache: Option<cache::Usage>,
    /// The placement cache the run read and filled, beside the evidence
    /// cache.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) placement_cache: Option<cache::Usage>,
    /// The cached index of the truth's written forms the run read or
    /// stored, with `--evidence-cache` and a build with the index's code
    /// identity.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) forms_cache: Option<cache::FormsUsage>,
    /// The run contexts the consensus was calibrated by.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) run_context: Option<RunContextHeader>,
    /// The widths that place the review zones, in positions.
    pub(crate) zones: ZoneHeader,
    /// The harness's own measurement constants.
    pub(crate) constants: Constants,
    /// Outcome counts per measure.
    pub(crate) samples: BTreeMap<String, Counts>,
    /// What passes without review.
    pub(crate) workload: WorkloadHeader,
    pub(crate) seconds: Seconds,
    /// The SHA-256 of every other file of the run, by file name, so two
    /// runs' outputs compare without reading them (`dna-eval same`). A run
    /// of an earlier revision records none.
    #[serde(default)]
    pub(crate) outputs: BTreeMap<String, String>,
}

/// How the samples are grouped for their run signal (`--run-signal`,
/// PROP-0009): by batch, falling back to the batch month for a small batch,
/// or by month alone.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, clap::ValueEnum)]
#[serde(rename_all = "lowercase")]
pub(crate) enum SignalGrouping {
    Batch,
    Month,
}

/// The harness's measurement constants, recorded so that a run's measures can
/// be reproduced from its outputs. Shares are in thousandths.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub(crate) struct Constants {
    /// Limit of detection of a second base (`haplotype_lod`).
    pub(crate) lod_share: u16,
    /// Smallest share of a faint second base, for discrepancy causes.
    pub(crate) faint_share: u16,
    /// Normal quantile of the two-sided 95 % intervals and Wilson bounds.
    pub(crate) interval_z: f64,
}

impl Constants {
    pub(crate) const fn current() -> Self {
        Self {
            lod_share: conventions::LOD,
            faint_share: conventions::FAINT,
            interval_z: estimate::Z_95,
        }
    }
}

/// The run contexts' directory, how many of the panel's batches had one,
/// and the digest of the contexts read: the SHA-256 of one
/// `<batch>\t<file SHA-256>` line per context, in batch order.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub(crate) struct RunContextHeader {
    pub(crate) directory: PathBuf,
    pub(crate) batches: usize,
    pub(crate) sha256: String,
}

/// The re-check list: its path, its SHA-256, and its rows.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub(crate) struct RecheckHeader {
    pub(crate) path: PathBuf,
    pub(crate) sha256: String,
    pub(crate) rows: usize,
}

/// What passes without review over the samples the review workload counts,
/// and over the current data when the run names it.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub(crate) struct WorkloadHeader {
    pub(crate) all: Passes,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub(crate) current: Option<Passes>,
}

/// The widths that place the review zones, from the run's options.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub(crate) struct ZoneHeader {
    pub(crate) after: i64,
    pub(crate) edge: i64,
    pub(crate) overrun: i64,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub(crate) struct Source {
    /// Git revision of the checkout the harness was built from, if available.
    pub(crate) revision: Option<String>,
    /// Whether tracked files of that checkout differed from the revision.
    pub(crate) dirty: Option<bool>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub(crate) struct ManifestHeader {
    pub(crate) path: PathBuf,
    pub(crate) fields: BTreeMap<String, String>,
    pub(crate) body_sha256: String,
    pub(crate) rows: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub(crate) struct TruthHeader {
    pub(crate) path: PathBuf,
    /// The database file's SHA-256, for provenance only.
    pub(crate) sha256: String,
    /// The database's `info` table, for provenance only.
    pub(crate) info: BTreeMap<String, String>,
    /// The panel truth digest of the current database.
    pub(crate) digest: String,
    /// The panel truth digest the manifest records.
    pub(crate) recorded_digest: Option<String>,
    /// `unchanged`, `changed`, or `not_recorded`.
    pub(crate) state: String,
    /// Samples whose truth differs from the manifest's record.
    pub(crate) changed_samples: usize,
    pub(crate) panel_samples_found: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub(crate) struct Counts {
    pub(crate) total: usize,
    pub(crate) concordant: usize,
    pub(crate) discordant: usize,
    pub(crate) not_applicable: usize,
    pub(crate) failed: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub(crate) struct Seconds {
    pub(crate) total: f64,
    pub(crate) truth: f64,
    pub(crate) dna: f64,
    pub(crate) sample_median: f64,
    pub(crate) sample_p95: f64,
    pub(crate) sample_max: f64,
    pub(crate) sample_sum: f64,
}

impl RunHeader {
    /// `run.json`'s text: the header as pretty JSON and a newline.
    pub(crate) fn text(&self) -> Result<String> {
        let mut json = serde_json::to_string_pretty(self).map_err(Error::Serialize)?;
        json.push('\n');
        Ok(json)
    }

    /// Reads the header of the run in `directory`.
    pub(crate) fn read(directory: &Path) -> std::result::Result<Self, DiffError> {
        let path = directory.join(RUN_JSON);
        let bytes = fs::read(&path).map_err(|error| failure(&path, error))?;
        serde_json::from_slice(&bytes).map_err(|error| failure(&path, error))
    }
}

/// `checks.tsv`: every analysed sample's check items, in sample order.
pub(crate) fn checks(samples: &[Sample], outcomes: &[Outcome]) -> String {
    let mut text = String::from(CHECKS_HEADER);
    for (sample, outcome) in samples.iter().zip(outcomes) {
        if let Ok(called) = &outcome.result {
            for (first, last, labels) in &called.check_labels {
                let _ = writeln!(text, "{}\t{first}\t{last}\t{labels}", sample.id);
            }
        }
    }
    text
}

/// `notes.tsv`: every run of every `tract_length` note, in sample order.
pub(crate) fn notes(samples: &[Sample], outcomes: &[Outcome]) -> String {
    let mut text = String::from(NOTES_HEADER);
    let cell =
        |value: Option<usize>| value.map_or_else(|| ".".to_owned(), |value| value.to_string());
    for (sample, outcome) in samples.iter().zip(outcomes) {
        if let Ok(called) = &outcome.result {
            for note in &called.notes {
                for run in &note.runs {
                    let _ = writeln!(
                        text,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        sample.id,
                        note.first,
                        note.last,
                        note.tract,
                        run.run,
                        run.status,
                        cell(run.called),
                        cell(run.reading),
                        cell(run.convention),
                    );
                }
            }
        }
    }
    text
}

/// `signals.tsv`: one line per run-signal group, as the runner gives them.
pub(crate) fn signals(groups: &[String]) -> String {
    let mut text = String::from(SIGNALS_HEADER);
    for line in groups {
        text.push_str(line);
        text.push('\n');
    }
    text
}

/// `components.tsv`: every analysed sample's phase-mixture components,
/// measured against its truth on `reference`.
pub(crate) fn components(
    samples: &[Sample],
    outcomes: &[Outcome],
    truth: &BTreeMap<String, TruthSample>,
    reference: &[u8],
) -> String {
    let mut text = String::from(components::HEADER);
    for (sample, outcome) in samples.iter().zip(outcomes) {
        if let Ok(called) = &outcome.result {
            components::append(
                &mut text,
                (&sample.id, sample.batch_date.as_deref().unwrap_or("")),
                &called.components,
                truth
                    .get(&sample.id)
                    .map(|record| record.alleles.as_slice()),
                reference,
            );
        }
    }
    text
}

/// `tracts.tsv`: every analysed sample's tract decodes.
pub(crate) fn tracts(samples: &[Sample], outcomes: &[Outcome]) -> String {
    let mut text = String::from(tracts::HEADER);
    for (sample, outcome) in samples.iter().zip(outcomes) {
        if let Ok(called) = &outcome.result {
            tracts::append(
                &mut text,
                (&sample.id, sample.batch_date.as_deref().unwrap_or("")),
                &called.tracts,
            );
        }
    }
    text
}

/// `run_context.tsv`: run context coverage per batch month, then over all
/// samples: of the analysed samples, how many had their batch's context and
/// in how many it calibrated τ of run 1 and of run 2; the rest fail closed
/// to the configured τ. With them, the batches that had a context and their
/// pairs per run, which show why: a batch below the minimum pairs calibrates
/// nothing.
pub(crate) fn context_coverage(
    samples: &[Sample],
    outcomes: &[Outcome],
    contexts: &Contexts,
) -> String {
    let mut months = BTreeMap::<String, [usize; 8]>::new();
    let mut counted = BTreeSet::<(String, &str)>::new();
    for (sample, outcome) in samples.iter().zip(outcomes) {
        let month = sample
            .batch_month
            .map_or_else(|| "unknown".to_owned(), |month| month.to_string());
        for key in [month.as_str(), "all"] {
            let counts = months.entry(key.to_owned()).or_default();
            if let Some(context) = contexts.get(&sample.batch)
                && counted.insert((key.to_owned(), &sample.batch))
            {
                counts[5] += 1;
                counts[6] += context.length.runs[0].pairs;
                counts[7] += context.length.runs[1].pairs;
            }
            counts[0] += 1;
            if let Ok(called) = &outcome.result {
                counts[1] += 1;
                if let Some([run1, run2]) = called.length_context {
                    counts[2] += 1;
                    counts[3] += usize::from(run1);
                    counts[4] += usize::from(run2);
                }
            }
        }
    }
    let all = months.remove("all").unwrap_or_default();
    let mut text = String::from(RUN_CONTEXT_HEADER);
    for (month, counts) in months.iter().chain([(&"all".to_owned(), &all)]) {
        let [
            samples,
            analysed,
            context,
            run1,
            run2,
            batches,
            pairs1,
            pairs2,
        ] = counts;
        let _ = writeln!(
            text,
            "{month}\t{samples}\t{analysed}\t{context}\t{run1}\t{run2}\t{batches}\t{pairs1}\t{pairs2}"
        );
    }
    text
}

/// Creates the output directory, which must not hold files yet.
pub(crate) fn prepare_output(directory: &Path) -> Result<()> {
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

/// The SHA-256 of every file of the run in `directory` but `run.json`, by
/// file name.
///
/// # Errors
///
/// Fails when the directory cannot be read, holds anything but files, a
/// name that is not UTF-8, or a file that cannot be read.
pub(crate) fn digests(
    directory: &Path,
) -> std::result::Result<BTreeMap<String, String>, DiffError> {
    let mut digests = BTreeMap::new();
    for entry in fs::read_dir(directory).map_err(|error| failure(directory, error))? {
        let entry = entry.map_err(|error| failure(directory, error))?;
        let path = entry.path();
        let kind = entry.file_type().map_err(|error| failure(&path, error))?;
        if !kind.is_file() {
            return Err(failure(&path, "not a file of a run"));
        }
        let name = entry
            .file_name()
            .into_string()
            .map_err(|_| failure(&path, "the file name is not UTF-8"))?;
        if name != RUN_JSON {
            let sha256 =
                crate::runner::sha256_file(&path).map_err(|error| failure(&path, error))?;
            digests.insert(name, sha256);
        }
    }
    Ok(digests)
}

/// Writes `contents` as the file `name` of `directory`.
pub(crate) fn write(directory: &Path, name: &str, contents: &str) -> Result<()> {
    let path = directory.join(name);
    fs::write(&path, contents).map_err(|source| Error::Write { path, source })
}

fn failure(path: &Path, message: impl std::fmt::Display) -> DiffError {
    DiffError::Run {
        path: path.to_path_buf(),
        message: message.to_string(),
    }
}

/// One row of a table read by [`read_table`]: its line number and the cells
/// of the columns asked for, in the order asked.
pub(crate) struct Line {
    pub(crate) number: usize,
    pub(crate) cells: Vec<String>,
}

/// Reads the table `name` of the run in `directory`, keeping the cells of
/// `columns` in that order. Comment lines (`# `) may come before the
/// header only. The table must name every column of `columns`, and no
/// column twice; every row must have as many cells as the header.
pub(crate) fn read_table(
    directory: &Path,
    name: &str,
    columns: &[&str],
) -> std::result::Result<Vec<Line>, DiffError> {
    let path = directory.join(name);
    let text = fs::read_to_string(&path).map_err(|error| failure(&path, error))?;
    parse_table(&path, &text, columns)
}

fn parse_table(
    path: &Path,
    text: &str,
    columns: &[&str],
) -> std::result::Result<Vec<Line>, DiffError> {
    let mut lines = text
        .lines()
        .enumerate()
        .skip_while(|(_, line)| line.starts_with("# "));
    let (_, header) = lines.next().ok_or_else(|| failure(path, "no header"))?;
    let header = header.split('\t').collect::<Vec<_>>();
    let mut seen = BTreeSet::new();
    if let Some(twice) = header.iter().find(|name| !seen.insert(**name)) {
        return Err(failure(path, format!("the column {twice} appears twice")));
    }
    let indices = columns
        .iter()
        .map(|column| {
            header
                .iter()
                .position(|name| name == column)
                .ok_or_else(|| failure(path, format!("no column {column}")))
        })
        .collect::<std::result::Result<Vec<_>, _>>()?;
    lines
        .map(|(index, line)| {
            let number = index + 1;
            let cells = line.split('\t').collect::<Vec<_>>();
            if cells.len() != header.len() {
                return Err(failure(
                    path,
                    format!(
                        "line {number} has {} cells, the header {}",
                        cells.len(),
                        header.len()
                    ),
                ));
            }
            Ok(Line {
                number,
                cells: indices.iter().map(|&at| cells[at].to_owned()).collect(),
            })
        })
        .collect()
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::compare::Measure;
    use crate::{truth, workload, zones};

    /// A run header with every optional field.
    pub(crate) fn fixture() -> RunHeader {
        let digest = "a".repeat(64);
        let counts = || Counts {
            total: 3,
            concordant: 1,
            discordant: 1,
            not_applicable: 0,
            failed: 1,
        };
        let zones = zones::Zones::new(
            &crate::geometry::tests::spans(&[(16024, 16365), (73, 340)]),
            &[(16184, 16193)],
            zones::Widths {
                after: 20,
                edge: 40,
                overrun: 2,
            },
        );
        let passes = || workload::passes(std::iter::empty(), Some(&zones));
        RunHeader {
            harness: env!("CARGO_PKG_NAME").to_owned(),
            harness_version: env!("CARGO_PKG_VERSION").to_owned(),
            source: Source {
                revision: Some("0".repeat(40)),
                dirty: Some(false),
            },
            threads: 4,
            manifest: ManifestHeader {
                path: PathBuf::from("/panel.tsv"),
                fields: BTreeMap::from([("panel".to_owned(), "dev".to_owned())]),
                body_sha256: digest.clone(),
                rows: 3,
            },
            truth: TruthHeader {
                path: PathBuf::from("/statistics.sqlite"),
                sha256: digest.clone(),
                info: BTreeMap::from([("schema_version".to_owned(), "1".to_owned())]),
                digest: digest.clone(),
                recorded_digest: Some(digest.clone()),
                state: truth::TruthState::Unchanged.label().to_owned(),
                changed_samples: 0,
                panel_samples_found: 3,
            },
            configuration_path: PathBuf::from("/dna.toml"),
            reference_path: PathBuf::from("/rCRS.fasta"),
            identities: vec![Identity {
                configuration_sha256: digest.clone(),
                reference_name: "rCRS".into(),
                reference_sha256: digest.clone(),
                profile_id: "human-mtdna-rcrs".into(),
                profile_sha256: digest,
            }],
            review_window: conventions::REVIEW_WINDOW,
            current_since: "2026-01".parse().ok(),
            run_signal: Some(SignalGrouping::Batch),
            run_signal_batch_reads: Some(40),
            forced_windows: Some(PathBuf::from("/forced.txt")),
            recheck: Some(RecheckHeader {
                path: PathBuf::from("/recheck.tsv"),
                sha256: "b".repeat(64),
                rows: 2,
            }),
            evidence_cache: Some(cache::Usage {
                directory: "/cache/evidence".into(),
                key: cache::Key {
                    sanger: 17,
                    configuration_sha256: "a".repeat(64),
                    code_sha256: "b".repeat(64),
                },
                hits: 5,
                misses: 1,
                replaced: 0,
                bytes: 270_000,
            }),
            placement_cache: Some(cache::Usage {
                directory: "/cache/evidence".into(),
                key: cache::Key {
                    sanger: 17,
                    configuration_sha256: "a".repeat(64),
                    code_sha256: "c".repeat(64),
                },
                hits: 4,
                misses: 2,
                replaced: 0,
                bytes: 90_000,
            }),
            forms_cache: Some(cache::FormsUsage {
                directory: "/cache/evidence".into(),
                key: cache::FormsKey {
                    truth: "a".repeat(64),
                    reference: "a".repeat(64),
                    windows: "d".repeat(64),
                    code: "e".repeat(64),
                },
                hit: true,
            }),
            run_context: None,
            zones: ZoneHeader {
                after: 20,
                edge: 40,
                overrun: 2,
            },
            constants: Constants::current(),
            samples: Measure::ALL
                .into_iter()
                .map(|measure| (measure.label().to_owned(), counts()))
                .collect(),
            workload: WorkloadHeader {
                all: passes(),
                current: Some(passes()),
            },
            seconds: Seconds {
                total: 2.5,
                truth: 0.5,
                dna: 1.5,
                sample_median: 0.5,
                sample_p95: 0.75,
                sample_max: 1.0,
                sample_sum: 1.5,
            },
            outputs: BTreeMap::from([(SCORES.to_owned(), "c".repeat(64))]),
        }
    }

    /// `run.json` with every optional field, and without them, and the
    /// headers of `checks.tsv` and `signals.tsv`, follow `dna.eval_run/v1`.
    #[test]
    fn the_run_header_and_tables_follow_the_contract() -> crate::contract::Checked {
        let header = fixture();
        crate::contract::check_document("run.json", &serde_json::to_value(&header)?)?;
        // `diff` reads back what `run` writes.
        assert_eq!(serde_json::from_str::<RunHeader>(&header.text()?)?, header);
        let bare = RunHeader {
            source: Source {
                revision: None,
                dirty: None,
            },
            identities: Vec::new(),
            current_since: None,
            run_signal: None,
            run_signal_batch_reads: None,
            forced_windows: None,
            recheck: None,
            outputs: BTreeMap::new(),
            forms_cache: None,
            workload: WorkloadHeader {
                all: header.workload.all,
                current: None,
            },
            ..header
        };
        crate::contract::check_document("run.json", &serde_json::to_value(&bare)?)?;
        assert_eq!(serde_json::from_str::<RunHeader>(&bare.text()?)?, bare);
        // A run of an earlier revision, without output digests, reads back.
        let mut earlier = serde_json::to_value(&bare)?;
        earlier
            .as_object_mut()
            .map(|header| header.remove("outputs"));
        assert_eq!(serde_json::from_value::<RunHeader>(earlier)?, bare);
        crate::contract::check_table("checks.tsv", CHECKS_HEADER)?;
        crate::contract::check_table("notes.tsv", NOTES_HEADER)?;
        crate::contract::check_table("run_context.tsv", RUN_CONTEXT_HEADER)?;
        crate::contract::check_table("signals.tsv", SIGNALS_HEADER)
    }

    /// Every file but `run.json` is hashed; a directory is refused.
    #[test]
    fn hashes_every_file_of_a_run_but_its_header()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        let directory = tempfile::tempdir()?;
        write(directory.path(), RUN_JSON, "{}\n")?;
        write(directory.path(), SCORES, "abc")?;
        write(directory.path(), CHECKS, "")?;
        assert_eq!(
            digests(directory.path())?,
            BTreeMap::from([
                (
                    CHECKS.to_owned(),
                    "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855".to_owned()
                ),
                (
                    SCORES.to_owned(),
                    "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad".to_owned()
                ),
            ])
        );
        fs::create_dir(directory.path().join("nested"))?;
        assert!(matches!(
            digests(directory.path()),
            Err(DiffError::Run { message, .. }) if message == "not a file of a run"
        ));
        Ok(())
    }

    #[test]
    fn reads_the_columns_asked_for_and_refuses_a_malformed_table()
    -> std::result::Result<(), DiffError> {
        let path = Path::new("t.tsv");
        let lines = parse_table(path, "# c\na\tb\tc\n1\t2\t3\n4\t5\t6\n", &["c", "a"])?;
        assert_eq!(
            lines
                .iter()
                .map(|line| (line.number, line.cells.clone()))
                .collect::<Vec<_>>(),
            [
                (3, vec!["3".to_owned(), "1".to_owned()]),
                (4, vec!["6".to_owned(), "4".to_owned()])
            ]
        );
        for (text, message) in [
            ("a\tb\n1\t2\n", "no column c"),
            ("a\tc\ta\n1\t2\t3\n", "the column a appears twice"),
            ("a\tc\n1\n", "line 2 has 1 cells, the header 2"),
            ("", "no header"),
        ] {
            assert!(
                matches!(
                    parse_table(path, text, &["a", "c"]),
                    Err(DiffError::Run { message: found, .. }) if found == message
                ),
                "{text:?}"
            );
        }
        Ok(())
    }
}
