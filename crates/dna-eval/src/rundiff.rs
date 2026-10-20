//! `dna-eval diff`: what a candidate run changes against a base run of the
//! same panel, read from the two run directories alone, with a verdict.
//!
//! Both runs must measure the same thing: the same manifest body, truth
//! digest, re-check list, current-data month, review window, and zone
//! widths, each from a clean build; otherwise the comparison is refused.
//! The passes are counted from `scores.tsv` and `unreviewed.tsv` as the
//! run's own summary counts them: a sample, or a zone of it, passes when no
//! check item touches it, and a passed unit counts under its most severe
//! unreviewed block (non-het, then truth `N`, then het, then re-check).
//!
//! The gate fails when, in any period or source, more passed samples or
//! more passed zones of any zone carry a non-het unreviewed block, or more
//! samples carry any unreviewed block other than one the re-check list or a
//! truth `N` explains.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use serde_json::Value;

use crate::baseline::{self, Class as StatusClass};
use crate::compare::Measure;
use crate::error::DiffError;
use crate::manifest;

/// Rows of a sample-change listing shown in the summary; the TSV has all.
const SHOWN: usize = 40;

fn failure(path: &Path, message: impl std::fmt::Display) -> DiffError {
    DiffError::Run {
        path: path.to_path_buf(),
        message: message.to_string(),
    }
}

/// One unreviewed block of a sample, as `unreviewed.tsv` writes it.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct Escape {
    block: String,
    class: Severity,
    zone: Option<String>,
    zone_checked: bool,
}

/// The pass measures' classes, from the least to the most severe.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Severity {
    Recheck,
    Het,
    TruthN,
    NonHet,
}

impl Severity {
    fn parse(label: &str) -> Option<Self> {
        Some(match label {
            "recheck" => Self::Recheck,
            "het" => Self::Het,
            "truth_n" => Self::TruthN,
            "non_het" => Self::NonHet,
            _ => return None,
        })
    }

    const fn label(self) -> &'static str {
        match self {
            Self::Recheck => "recheck",
            Self::Het => "het",
            Self::TruthN => "truth_n",
            Self::NonHet => "non_het",
        }
    }
}

/// One sample of a run.
#[derive(Debug, Clone, Default)]
struct Sample {
    /// Counted by the pass measures: compared under `haplotype_lod`.
    counted: bool,
    clear: bool,
    positions: i64,
    zones: BTreeSet<String>,
    escapes: Vec<Escape>,
    statuses: BTreeMap<Measure, String>,
    uncovered: BTreeSet<String>,
}

/// One run directory, read.
struct Run {
    directory: PathBuf,
    header: Value,
    samples: BTreeMap<String, Sample>,
}

/// Reads a TSV with a header row into rows of named cells.
fn table(path: &Path) -> Result<Vec<BTreeMap<String, String>>, DiffError> {
    let text = std::fs::read_to_string(path).map_err(|error| failure(path, error))?;
    let mut lines = text.lines().filter(|line| !line.starts_with("# "));
    let columns = lines
        .next()
        .ok_or_else(|| failure(path, "no header"))?
        .split('\t')
        .map(str::to_owned)
        .collect::<Vec<_>>();
    lines
        .map(|line| {
            let cells = line.split('\t').collect::<Vec<_>>();
            if cells.len() == columns.len() {
                Ok(columns
                    .iter()
                    .cloned()
                    .zip(cells.into_iter().map(str::to_owned))
                    .collect())
            } else {
                Err(failure(path, "a row does not match the header"))
            }
        })
        .collect()
}

impl Run {
    fn read(directory: &Path) -> Result<Self, DiffError> {
        let path = directory.join("run.json");
        let header = serde_json::from_slice::<Value>(
            &std::fs::read(&path).map_err(|error| failure(&path, error))?,
        )
        .map_err(|error| failure(&path, error))?;
        let scores_path = directory.join(baseline::SCORES);
        let statuses = baseline::parse(
            &std::fs::read(&scores_path).map_err(|error| failure(&scores_path, error))?,
        )
        .map_err(|error| failure(&scores_path, error))?;
        let mut samples = BTreeMap::new();
        for row in table(&scores_path)? {
            let cell = |name: &str| row.get(name).map_or("", String::as_str);
            let id = cell("sample_id").to_owned();
            let lod = cell("haplotype_lod_status");
            let sample = Sample {
                counted: lod == "concordant" || lod == "discordant",
                clear: cell("review_clear") == "Y",
                positions: cell("review_positions").parse().unwrap_or(0),
                zones: cell("review_zones")
                    .split_whitespace()
                    .map(str::to_owned)
                    .collect(),
                escapes: Vec::new(),
                statuses: statuses.get(&id).cloned().unwrap_or_default(),
                uncovered: cell("uncovered")
                    .split_whitespace()
                    .map(str::to_owned)
                    .collect(),
            };
            samples.insert(id, sample);
        }
        let unreviewed_path = directory.join("unreviewed.tsv");
        for row in table(&unreviewed_path)? {
            let cell = |name: &str| row.get(name).map_or("", String::as_str);
            let class = Severity::parse(cell("class"))
                .ok_or_else(|| failure(&unreviewed_path, "an unknown class"))?;
            let zone = Some(cell("zone"))
                .filter(|zone| *zone != ".")
                .map(str::to_owned);
            let escape = Escape {
                block: cell("block").to_owned(),
                class,
                zone,
                zone_checked: cell("zone_checked") == "yes",
            };
            samples
                .get_mut(cell("sample_id"))
                .ok_or_else(|| failure(&unreviewed_path, "a sample scores.tsv does not list"))?
                .escapes
                .push(escape);
        }
        Ok(Self {
            directory: directory.to_path_buf(),
            header,
            samples,
        })
    }

    fn field(&self, pointer: &str) -> &Value {
        self.header.pointer(pointer).unwrap_or(&Value::Null)
    }
}

/// What makes two runs measure the same thing, as each records it.
fn mismatches(base: &Run, run: &Run) -> Vec<String> {
    let mut found = Vec::new();
    for (name, pointer) in [
        ("manifest body", "/manifest/body_sha256"),
        ("truth digest", "/truth/digest"),
        ("re-check list", "/recheck/sha256"),
        ("current-data month", "/current_since"),
        ("review window", "/review_window"),
        ("zone widths", "/zones"),
    ] {
        if base.field(pointer) != run.field(pointer) {
            found.push(format!(
                "{name} differs ({} against {})",
                base.field(pointer),
                run.field(pointer)
            ));
        }
    }
    for (label, run) in [("base", base), ("run", run)] {
        if run.field("/source/dirty") != &Value::Bool(false) {
            found.push(format!(
                "the {label} run was not built from a clean commit ({})",
                run.directory.display()
            ));
        }
    }
    found
}

/// The groups a sample counts in: all, its period, and its source.
fn groups(
    id: &str,
    periods: &BTreeMap<String, String>,
    sources: &BTreeMap<String, String>,
) -> Vec<String> {
    let mut groups = vec!["all".to_owned()];
    groups.extend(periods.get(id).map(|period| format!("period {period}")));
    groups.extend(sources.get(id).map(|source| format!("source {source}")));
    groups
}

/// One group's sample-level counts.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
struct Counts {
    samples: usize,
    clear: usize,
    /// Samples with an unreviewed block the re-check list or a truth `N`
    /// does not explain.
    unreviewed: usize,
    /// Passed samples by their most severe unreviewed block.
    passed_non_het: usize,
    passed_truth_n: usize,
    passed_het_only: usize,
    passed_recheck_only: usize,
    positions: i64,
}

/// One zone's counts in one group: units needing a check, and passed units
/// with a non-het unreviewed block.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
struct ZoneCounts {
    checked: usize,
    passed_non_het: usize,
}

/// What a run gives per group, and per group and zone.
struct Tally {
    counts: BTreeMap<String, Counts>,
    zones: BTreeMap<(String, String), ZoneCounts>,
}

fn tally(
    run: &Run,
    periods: &BTreeMap<String, String>,
    sources: &BTreeMap<String, String>,
) -> Tally {
    let mut counts = BTreeMap::<String, Counts>::new();
    let mut zones = BTreeMap::<(String, String), ZoneCounts>::new();
    for (id, sample) in &run.samples {
        if !sample.counted {
            continue;
        }
        let worst = sample.escapes.iter().map(|escape| escape.class).max();
        let open = sample
            .escapes
            .iter()
            .any(|escape| matches!(escape.class, Severity::NonHet | Severity::Het));
        let mut passed_zones = BTreeMap::<&str, Severity>::new();
        for escape in &sample.escapes {
            if let Some(zone) = escape.zone.as_deref().filter(|_| !escape.zone_checked) {
                let entry = passed_zones.entry(zone).or_insert(escape.class);
                *entry = (*entry).max(escape.class);
            }
        }
        for group in groups(id, periods, sources) {
            let count = counts.entry(group.clone()).or_default();
            count.samples += 1;
            count.positions += sample.positions;
            count.unreviewed += usize::from(open);
            if sample.clear {
                count.clear += 1;
                match worst {
                    Some(Severity::NonHet) => count.passed_non_het += 1,
                    Some(Severity::TruthN) => count.passed_truth_n += 1,
                    Some(Severity::Het) => count.passed_het_only += 1,
                    Some(Severity::Recheck) => count.passed_recheck_only += 1,
                    None => {}
                }
            }
            for zone in &sample.zones {
                zones
                    .entry((group.clone(), zone.clone()))
                    .or_default()
                    .checked += 1;
            }
            for (zone, severity) in &passed_zones {
                let entry = zones
                    .entry((group.clone(), (*zone).to_owned()))
                    .or_default();
                entry.passed_non_het += usize::from(*severity == Severity::NonHet);
            }
        }
    }
    Tally { counts, zones }
}

/// Checks that the files count the passes the run's own `run.json` records
/// over all samples: the auto-passed samples by their most severe block, and
/// every zone's checked units and passed units with a non-het block.
fn consistent(run: &Run, tally: &Tally) -> Result<(), DiffError> {
    let all = tally.counts.get("all").copied().unwrap_or_default();
    let workload = run.field("/workload/all");
    let number = |value: &Value| value.as_u64().and_then(|value| usize::try_from(value).ok());
    let mut found = Vec::new();
    for (name, recorded, counted) in [
        (
            "auto-pass samples",
            &workload["auto_pass"]["samples"],
            all.clear,
        ),
        (
            "auto-pass non-het",
            &workload["auto_pass"]["unreviewed"]["non_het"],
            all.passed_non_het,
        ),
        (
            "auto-pass het only",
            &workload["auto_pass"]["unreviewed"]["het_only"],
            all.passed_het_only,
        ),
    ] {
        if number(recorded) != Some(counted) {
            found.push(format!("{name}: {recorded} recorded, {counted} counted"));
        }
    }
    for zone in workload["zones"]["by_zone"]
        .as_array()
        .into_iter()
        .flatten()
    {
        let name = zone["name"].as_str().unwrap_or_default().to_owned();
        let counted = tally
            .zones
            .get(&("all".to_owned(), name.clone()))
            .copied()
            .unwrap_or_default();
        if number(&zone["checked"]) != Some(counted.checked)
            || number(&zone["unreviewed"]["non_het"]) != Some(counted.passed_non_het)
        {
            found.push(format!("zone {name}"));
        }
    }
    if found.is_empty() {
        Ok(())
    } else {
        Err(DiffError::Inconsistent(
            run.directory.clone(),
            found.join("; "),
        ))
    }
}

/// Each sample's period, by its manifest `batch_date` against the run's
/// current-data month, and its source, when the manifest has the column.
/// Each sample's label in one split.
type Labels = BTreeMap<String, String>;

fn split(run: &Run) -> Result<(Labels, Labels), DiffError> {
    let path = PathBuf::from(run.field("/manifest/path").as_str().unwrap_or_default());
    let bytes = std::fs::read(&path).map_err(|error| failure(&path, error))?;
    let manifest = manifest::parse(&bytes).map_err(|error| failure(&path, error))?;
    let since = run.field("/current_since").as_str().map(str::to_owned);
    let source = manifest.columns.iter().position(|name| name == "source");
    let mut periods = BTreeMap::new();
    let mut sources = BTreeMap::new();
    for sample in &manifest.samples {
        if let (Some(since), Some(date)) = (&since, &sample.batch_date) {
            let period = if date.as_str() >= since.as_str() {
                format!("current (≥ {since})")
            } else {
                format!("earlier (< {since})")
            };
            periods.insert(sample.id.clone(), period);
        }
        if let Some(value) = source.and_then(|index| sample.fields.get(index)) {
            sources.insert(sample.id.clone(), value.clone());
        }
    }
    Ok((periods, sources))
}

/// The comparison: the report, the per-sample change table, the status
/// change table, and whether the gate fails.
pub(crate) struct Comparison {
    pub(crate) report: String,
    pub(crate) samples: String,
    pub(crate) statuses: String,
    pub(crate) gate_failed: bool,
}

/// Compares the run at `run` against the base run at `base`.
///
/// # Errors
///
/// Returns [`DiffError`] when a run file cannot be read, or when the runs
/// measure different things and `force` is not set.
#[allow(clippy::too_many_lines)]
pub(crate) fn compare(base: &Path, run: &Path, force: bool) -> Result<Comparison, DiffError> {
    let base = Run::read(base)?;
    let run = Run::read(run)?;
    let mismatched = mismatches(&base, &run);
    if !mismatched.is_empty() && !force {
        return Err(DiffError::Mismatch(mismatched.join("; ")));
    }
    let (periods, sources) = split(&base)?;
    let before = tally(&base, &periods, &sources);
    let after = tally(&run, &periods, &sources);
    consistent(&base, &before)?;
    consistent(&run, &after)?;
    let mut report = String::from("# Run comparison\n\n");
    let revision = |run: &Run| {
        run.field("/source/revision")
            .as_str()
            .map_or("unknown", |revision| revision.get(..12).unwrap_or(revision))
            .to_owned()
    };
    let _ = writeln!(
        report,
        "Base `{}` ({}) against run `{}` ({}).\n",
        base.directory.display(),
        revision(&base),
        run.directory.display(),
        revision(&run)
    );
    for mismatch in &mismatched {
        let _ = writeln!(report, "**WARNING, not comparable:** {mismatch}.\n");
    }
    let mut failures = Vec::new();

    let _ = writeln!(
        report,
        "## Samples\n\nUnreviewed counts samples with a non-het or het unreviewed block; passed samples are counted under their most severe block.\n\n| Group | Samples | Clear | Unreviewed | Passed non-het | Passed truth-N | Passed het only | Passed re-check only | Positions to review |\n|---|---:|---:|---:|---:|---:|---:|---:|---:|"
    );
    let names = before
        .counts
        .keys()
        .chain(after.counts.keys())
        .cloned()
        .collect::<BTreeSet<_>>();
    for name in &names {
        let old = before.counts.get(name).copied().unwrap_or_default();
        let new = after.counts.get(name).copied().unwrap_or_default();
        let cell = |old: usize, new: usize| {
            if old == new {
                old.to_string()
            } else {
                format!("{old} → {new}")
            }
        };
        let _ = writeln!(
            report,
            "| {name} | {} | {} | {} | {} | {} | {} | {} | {} |",
            cell(old.samples, new.samples),
            cell(old.clear, new.clear),
            cell(old.unreviewed, new.unreviewed),
            cell(old.passed_non_het, new.passed_non_het),
            cell(old.passed_truth_n, new.passed_truth_n),
            cell(old.passed_het_only, new.passed_het_only),
            cell(old.passed_recheck_only, new.passed_recheck_only),
            if old.positions == new.positions {
                old.positions.to_string()
            } else {
                format!("{} → {}", old.positions, new.positions)
            },
        );
        if new.passed_non_het > old.passed_non_het {
            failures.push(format!(
                "passed samples with a non-het block rise in {name}: {} → {}",
                old.passed_non_het, new.passed_non_het
            ));
        }
        if new.unreviewed > old.unreviewed {
            failures.push(format!(
                "samples with an unreviewed block rise in {name}: {} → {}",
                old.unreviewed, new.unreviewed
            ));
        }
    }

    let _ = writeln!(
        report,
        "\n## Zones\n\nZones that changed: units needing a check, and passed units with a non-het unreviewed block.\n\n| Group | Zone | Needing a check | Passed non-het |\n|---|---|---:|---:|"
    );
    let keys = before
        .zones
        .keys()
        .chain(after.zones.keys())
        .cloned()
        .collect::<BTreeSet<_>>();
    for key in &keys {
        let old = before.zones.get(key).copied().unwrap_or_default();
        let new = after.zones.get(key).copied().unwrap_or_default();
        if old == new {
            continue;
        }
        let _ = writeln!(
            report,
            "| {} | {} | {} → {} | {} → {} |",
            key.0, key.1, old.checked, new.checked, old.passed_non_het, new.passed_non_het
        );
        if new.passed_non_het > old.passed_non_het {
            failures.push(format!(
                "passed {} zones with a non-het block rise in {}: {} → {}",
                key.1, key.0, old.passed_non_het, new.passed_non_het
            ));
        }
    }

    // Uncovered truth alleles that are new or gone, by allele.
    let mut uncovered = BTreeMap::<String, (usize, usize)>::new();
    for (id, sample) in &run.samples {
        let old = base
            .samples
            .get(id)
            .map(|sample| &sample.uncovered)
            .cloned()
            .unwrap_or_default();
        for allele in sample.uncovered.difference(&old) {
            uncovered.entry(allele.clone()).or_default().0 += 1;
        }
        for allele in old.difference(&sample.uncovered) {
            uncovered.entry(allele.clone()).or_default().1 += 1;
        }
    }
    let _ = writeln!(
        report,
        "\n## Coverage\n\nTruth alleles DNA's callable intervals newly leave out, or now cover, by allele.\n\n| Allele | Newly uncovered | Now covered |\n|---|---:|---:|"
    );
    let mut alleles = uncovered.into_iter().collect::<Vec<_>>();
    alleles.sort_by(|left, right| {
        (right.1.0 + right.1.1)
            .cmp(&(left.1.0 + left.1.1))
            .then(left.0.cmp(&right.0))
    });
    for (allele, (new, gone)) in alleles.iter().take(SHOWN) {
        let _ = writeln!(report, "| {allele} | {new} | {gone} |");
    }

    // Status changes per measure.
    let mut statuses = String::from("sample_id\tgroup\tmeasure\tbefore\tafter\n");
    let mut changes = BTreeMap::<(Measure, String, String), usize>::new();
    for (id, sample) in &run.samples {
        let Some(old) = base.samples.get(id) else {
            continue;
        };
        for measure in Measure::ALL {
            let (Some(before), Some(after)) =
                (old.statuses.get(&measure), sample.statuses.get(&measure))
            else {
                continue;
            };
            if StatusClass::of(before) != StatusClass::of(after) {
                *changes
                    .entry((measure, before.clone(), after.clone()))
                    .or_default() += 1;
                let _ = writeln!(
                    statuses,
                    "{id}\t{}\t{}\t{before}\t{after}",
                    periods.get(id).map_or("", String::as_str),
                    measure.label()
                );
            }
        }
    }
    let _ = writeln!(
        report,
        "\n## Status changes\n\n| Measure | Before | After | Samples |\n|---|---|---|---:|"
    );
    for ((measure, before, after), count) in &changes {
        let _ = writeln!(
            report,
            "| {} | {before} | {after} | {count} |",
            measure.label()
        );
    }

    // Per-sample changes of the unreviewed blocks or the clear state.
    let mut samples =
        String::from("sample_id\tperiod\tsource\tclear_before\tclear_after\tgone\tnew\n");
    let mut listed = Vec::new();
    for (id, sample) in &run.samples {
        let Some(old) = base.samples.get(id) else {
            continue;
        };
        let text = |escapes: &[Escape], other: &[Escape]| {
            escapes
                .iter()
                .filter(|escape| !other.contains(escape))
                .map(|escape| {
                    format!(
                        "{} {} {}",
                        escape.block,
                        escape.class.label(),
                        escape.zone.as_deref().unwrap_or(".")
                    )
                })
                .collect::<Vec<_>>()
                .join("; ")
        };
        let gone = text(&old.escapes, &sample.escapes);
        let new = text(&sample.escapes, &old.escapes);
        if gone.is_empty() && new.is_empty() && old.clear == sample.clear {
            continue;
        }
        let flag = |clear: bool| if clear { "Y" } else { "N" };
        let row = format!(
            "{id}\t{}\t{}\t{}\t{}\t{gone}\t{new}",
            periods.get(id).map_or("", String::as_str),
            sources.get(id).map_or("", String::as_str),
            flag(old.clear),
            flag(sample.clear),
        );
        let _ = writeln!(samples, "{row}");
        listed.push(row);
    }
    let _ = writeln!(
        report,
        "\n## Samples that changed\n\n{} samples changed their unreviewed blocks or their clear state; the first {} (all in `samples.tsv`):\n\n```",
        listed.len(),
        listed.len().min(SHOWN)
    );
    for row in listed.iter().take(SHOWN) {
        let _ = writeln!(report, "{row}");
    }
    report.push_str("```\n\n## Verdict\n\n");
    if failures.is_empty() {
        report.push_str("GATE PASS\n");
    } else {
        report.push_str("GATE FAIL\n\n");
        for failure in &failures {
            let _ = writeln!(report, "- {failure}");
        }
    }
    Ok(Comparison {
        report,
        samples,
        statuses,
        gate_failed: !failures.is_empty(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Writes a run directory: its manifest, its `run.json` with the passes
    /// `workload` records, its scores, and its unreviewed blocks.
    fn run(
        root: &Path,
        name: &str,
        dirty: bool,
        scores: &[&str],
        unreviewed: &[&str],
        workload: (usize, usize, usize),
    ) -> Result<PathBuf, Box<dyn std::error::Error>> {
        let directory = root.join(name);
        std::fs::create_dir_all(&directory)?;
        let manifest = root.join("panel.tsv");
        std::fs::write(
            &manifest,
            manifest::tests::manifest(&[
                "S1\tb1\tclean\t1\t1\t0\t/t/a.ab1",
                "S2\tb1\tclean\t1\t1\t0\t/t/b.ab1",
            ]),
        )?;
        let header = serde_json::json!({
            "source": {"revision": "0".repeat(40), "dirty": dirty},
            "manifest": {"path": manifest, "body_sha256": "b"},
            "truth": {"digest": "t"},
            "review_window": 10,
            "zones": {"after": 20, "edge": 40, "overrun": 2},
            "workload": {"all": {
                "auto_pass": {"samples": workload.0, "unreviewed": {"non_het": workload.1, "het_only": workload.2}},
                "zones": {"by_zone": []}
            }}
        });
        std::fs::write(directory.join("run.json"), serde_json::to_vec(&header)?)?;
        let mut text = String::from(
            "sample_id\tstatus\thaplotype_lod_status\treview_clear\treview_positions\treview_zones\tuncovered\n",
        );
        for row in scores {
            text.push_str(row);
            text.push('\n');
        }
        std::fs::write(directory.join("scores.tsv"), text)?;
        let mut text = String::from(crate::report::UNREVIEWED_HEADER);
        for row in unreviewed {
            text.push_str(row);
            text.push('\n');
        }
        std::fs::write(directory.join("unreviewed.tsv"), text)?;
        Ok(directory)
    }

    #[test]
    fn gates_a_new_non_het_escape() -> Result<(), Box<dyn std::error::Error>> {
        let root = tempfile::tempdir()?;
        let base = run(
            root.path(),
            "base",
            false,
            &[
                "S1\tconcordant\tconcordant\tY\t0\t\t",
                "S2\tdiscordant\tdiscordant\tY\t0\t\t",
            ],
            &["S2\t461\t461\t461-461:C/S\thet\twindow@HV2\tno"],
            (2, 0, 1),
        )?;
        let candidate = run(
            root.path(),
            "run",
            false,
            &[
                "S1\tdiscordant\tdiscordant\tY\t0\t\t73G",
                "S2\tdiscordant\tdiscordant\tN\t10\tstart@HV2\t",
            ],
            &[
                "S1\t309\t309\t309-309+:CC/C\tnon_het\twindow@HV2\tno",
                "S2\t461\t461\t461-461:C/S\thet\twindow@HV2\tno",
            ],
            (1, 1, 0),
        )?;
        let comparison = compare(&base, &candidate, false)?;
        assert!(comparison.gate_failed);
        assert!(comparison.report.contains("GATE FAIL"));
        assert!(
            comparison
                .report
                .contains("| all | 2 | 2 → 1 | 1 → 2 | 0 → 1 |")
        );
        assert!(
            comparison
                .report
                .contains("| all | window@HV2 | 0 → 0 | 0 → 1 |")
        );
        assert!(comparison.report.contains("| 73G | 1 | 0 |"));
        assert!(
            comparison
                .statuses
                .contains("S1\t\tprimary\tconcordant\tdiscordant")
        );
        assert!(
            comparison
                .samples
                .contains("S1\t\t\tY\tY\t\t309-309+:CC/C non_het window@HV2")
        );
        let same = compare(&base, &base, false)?;
        assert!(!same.gate_failed);
        assert!(same.report.contains("GATE PASS"));
        Ok(())
    }

    #[test]
    fn refuses_runs_that_measure_different_things() -> Result<(), Box<dyn std::error::Error>> {
        let root = tempfile::tempdir()?;
        let rows = ["S1\tconcordant\tconcordant\tY\t0\t\t"];
        let base = run(root.path(), "base", false, &rows, &[], (1, 0, 0))?;
        let dirty = run(root.path(), "dirty", true, &rows, &[], (1, 0, 0))?;
        assert!(matches!(
            compare(&base, &dirty, false),
            Err(DiffError::Mismatch(_))
        ));
        assert!(
            compare(&base, &dirty, true)?
                .report
                .contains("WARNING, not comparable")
        );
        let miscounted = run(root.path(), "miscounted", false, &rows, &[], (0, 0, 0))?;
        assert!(matches!(
            compare(&base, &miscounted, false),
            Err(DiffError::Inconsistent(..))
        ));
        Ok(())
    }
}
