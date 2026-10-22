//! `dna-eval diff`: what a candidate run changes against a base run of the
//! same panel, read from the two run directories alone, with a verdict.
//!
//! Both runs must measure the same thing: the same manifest body, truth
//! digest, re-check list, current-data month, review window, zone widths,
//! reference, profile, run signal grouping and batch reads, and run
//! contexts, each from a clean build; otherwise the comparison is refused.
//! The manifest the base run names is read again for the periods and
//! sources, and must still hash to the body the run recorded.
//! The run directories are read strictly ([`rundir`]): a table without a
//! column `diff` needs, or with a sample twice, is refused. The passes are
//! counted from `scores.tsv` and `unreviewed.tsv` with the review workload's
//! own counting ([`crate::workload::PassSample`]): a sample, or a zone of it,
//! passes when no check item touches it, and a passed unit counts under its
//! most severe unreviewed block (non-het, then truth `N`, then het, then
//! re-check, then length note). Before comparing, every count of a run is
//! checked against the passes its `run.json` records, over all samples and
//! over the current data.
//!
//! The gate reads one group: the current period, the samples whose manifest
//! `batch_date` is in or after the runs' current-data month, or all samples
//! when the runs record no such month. It fails when, in that group, more
//! passed samples or more passed zones of any zone carry a non-het
//! unreviewed block, or more samples carry a non-het unreviewed block, or
//! when the two runs failed to score different samples of it. The earlier
//! period and the sources, which mix periods, are reported in full and never
//! gated (user decision 2026-10-21: the LAB keeps improving, so recent data
//! decides). Heteroplasmy is optional (user decision): het-only, truth `N`,
//! and re-check moves are reported, never gated.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use serde_json::Value;

use crate::classes::StatusClass;
use crate::compare::Measure;
use crate::error::DiffError;
use crate::manifest::{self, Month};
use crate::rundir::{self, RunHeader};
use crate::workload::{AutoPass, Class, Escaped, PassSample, Passes, Unit};

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
    class: Class,
    zone: Option<String>,
    zone_checked: bool,
}

/// One sample of a run.
#[derive(Debug, Clone)]
struct Sample {
    /// Counted by the pass measures: compared under `haplotype_lod`.
    counted: bool,
    /// The harness could not score it (`status` `failed`).
    failed: bool,
    positions: i64,
    /// The sample as the pass measures count it, its zones by label.
    pass: PassSample<String>,
    escapes: Vec<Escape>,
    statuses: BTreeMap<Measure, String>,
    uncovered: BTreeSet<String>,
}

/// One run directory, read.
struct Run {
    directory: PathBuf,
    header: RunHeader,
    samples: BTreeMap<String, Sample>,
}

/// The `scores.tsv` columns `diff` reads besides the measures' statuses.
const SCORE_COLUMNS: [&str; 5] = [
    "sample_id",
    "review_clear",
    "review_positions",
    "review_zones",
    "uncovered",
];

/// The `unreviewed.tsv` columns `diff` reads.
const UNREVIEWED_COLUMNS: [&str; 5] = ["sample_id", "block", "class", "zone", "zone_checked"];

impl Run {
    fn read(directory: &Path) -> Result<Self, DiffError> {
        let header = RunHeader::read(directory)?;
        let scores_path = directory.join(rundir::SCORES);
        let columns = SCORE_COLUMNS
            .into_iter()
            .chain(Measure::ALL.map(Measure::column))
            .collect::<Vec<_>>();
        let mut samples = BTreeMap::new();
        for line in rundir::read_table(directory, rundir::SCORES, &columns)? {
            let invalid =
                |message: &str| failure(&scores_path, format!("line {}: {message}", line.number));
            let [id, clear, positions, zones, uncovered] =
                [0, 1, 2, 3, 4].map(|at| &line.cells[at]);
            let statuses = Measure::ALL
                .into_iter()
                .zip(&line.cells[SCORE_COLUMNS.len()..])
                .map(|(measure, status)| (measure, status.clone()))
                .collect::<BTreeMap<_, _>>();
            let lod = statuses[&Measure::HaplotypeLod].as_str();
            let counted = lod == "concordant" || lod == "discordant";
            // A compared sample has a review workload; a failed one has
            // none, and the pass measures count only a counted one.
            let clear = match clear.as_str() {
                "Y" => true,
                "N" => false,
                "" if !counted => false,
                _ => return Err(invalid("review_clear is neither Y nor N")),
            };
            let positions = match positions.as_str() {
                "" if !counted => 0,
                text => text
                    .parse()
                    .map_err(|_| invalid("review_positions is not a number"))?,
            };
            let sample = Sample {
                counted,
                failed: statuses[&Measure::Primary] == "failed",
                positions,
                pass: PassSample {
                    clear,
                    zones: zones.split_whitespace().map(str::to_owned).collect(),
                    escapes: Vec::new(),
                },
                escapes: Vec::new(),
                statuses,
                uncovered: uncovered.split_whitespace().map(str::to_owned).collect(),
            };
            if samples.insert(id.clone(), sample).is_some() {
                return Err(invalid(&format!("sample {id} appears twice")));
            }
        }
        let unreviewed_path = directory.join(rundir::UNREVIEWED);
        for line in rundir::read_table(directory, rundir::UNREVIEWED, &UNREVIEWED_COLUMNS)? {
            let invalid = |message: &str| {
                failure(&unreviewed_path, format!("line {}: {message}", line.number))
            };
            let [id, block, class, zone, zone_checked] = [0, 1, 2, 3, 4].map(|at| &line.cells[at]);
            let class = Class::parse(class).ok_or_else(|| invalid("an unknown class"))?;
            let zone = Some(zone.as_str())
                .filter(|zone| *zone != ".")
                .map(str::to_owned);
            let zone_checked = match zone_checked.as_str() {
                "yes" => true,
                "no" => false,
                _ => return Err(invalid("zone_checked is neither yes nor no")),
            };
            let sample = samples
                .get_mut(id)
                .ok_or_else(|| invalid("a sample scores.tsv does not list"))?;
            if zone
                .as_ref()
                .is_some_and(|zone| sample.pass.zones.contains(zone))
                != zone_checked
            {
                return Err(invalid(
                    "zone_checked disagrees with the sample's review_zones",
                ));
            }
            sample.pass.escapes.push((zone.clone(), class));
            sample.escapes.push(Escape {
                block: block.clone(),
                class,
                zone,
                zone_checked,
            });
        }
        Ok(Self {
            directory: directory.to_path_buf(),
            header,
            samples,
        })
    }

    /// The samples the harness could not score.
    fn failed(&self) -> BTreeSet<&str> {
        self.samples
            .iter()
            .filter(|(_, sample)| sample.failed)
            .map(|(id, _)| id.as_str())
            .collect()
    }

    /// Every zone the run's files or header name.
    fn zone_names(&self) -> BTreeSet<String> {
        let recorded = [
            Some(&self.header.workload.all),
            self.header.workload.current.as_ref(),
        ]
        .into_iter()
        .flatten()
        .flat_map(|passes| &passes.zones.by_zone)
        .map(|unit| unit.name.clone());
        let seen = self.samples.values().flat_map(|sample| {
            sample.pass.zones.iter().cloned().chain(
                sample
                    .pass
                    .escapes
                    .iter()
                    .filter_map(|(zone, _)| zone.clone()),
            )
        });
        recorded.chain(seen).collect()
    }
}

/// What makes two runs measure the same thing, as each records it.
fn mismatches(base: &Run, run: &Run) -> Vec<String> {
    let value = |value: Result<Value, serde_json::Error>| value.unwrap_or(Value::Null);
    let fields = |run: &Run| {
        let header = &run.header;
        [
            (
                "manifest body",
                Value::from(header.manifest.body_sha256.clone()),
            ),
            ("truth digest", Value::from(header.truth.digest.clone())),
            (
                "re-check list",
                header
                    .recheck
                    .as_ref()
                    .map_or(Value::Null, |recheck| Value::from(recheck.sha256.clone())),
            ),
            (
                "current-data month",
                value(serde_json::to_value(header.current_since)),
            ),
            ("review window", Value::from(header.review_window)),
            ("zone widths", value(serde_json::to_value(header.zones))),
            (
                "reference",
                header
                    .identities
                    .iter()
                    .map(|identity| identity.reference_sha256.clone())
                    .collect::<BTreeSet<_>>()
                    .into_iter()
                    .collect(),
            ),
            (
                "profile",
                header
                    .identities
                    .iter()
                    .map(|identity| identity.profile_sha256.clone())
                    .collect::<BTreeSet<_>>()
                    .into_iter()
                    .collect(),
            ),
            ("run signal", value(serde_json::to_value(header.run_signal))),
            (
                "run signal batch reads",
                value(serde_json::to_value(header.run_signal_batch_reads)),
            ),
            (
                "run contexts",
                header
                    .run_context
                    .as_ref()
                    .map_or(Value::Null, |context| Value::from(context.sha256.clone())),
            ),
        ]
    };
    let mut found = Vec::new();
    for ((name, before), (_, after)) in fields(base).into_iter().zip(fields(run)) {
        if before != after {
            found.push(format!("{name} differs ({before} against {after})"));
        }
    }
    for (label, run) in [("base", base), ("run", run)] {
        if run.header.source.dirty != Some(false) {
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
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Counts {
    samples: usize,
    /// Samples with a non-het unreviewed block.
    unreviewed_non_het: usize,
    /// Samples whose only unreviewed blocks the re-check list or a truth
    /// `N` does not explain are het.
    unreviewed_het_only: usize,
    /// The samples without a check item, and those of them with an
    /// unreviewed block by its most severe class.
    auto_pass: AutoPass,
    positions: i64,
}

impl Default for Counts {
    fn default() -> Self {
        Self {
            samples: 0,
            unreviewed_non_het: 0,
            unreviewed_het_only: 0,
            auto_pass: AutoPass {
                samples: 0,
                unreviewed: Escaped::default(),
            },
            positions: 0,
        }
    }
}

/// What a run gives per group, and per group and zone.
struct Tally {
    counts: BTreeMap<String, Counts>,
    zones: BTreeMap<(String, String), Unit>,
}

/// Counts the passes of `run` per group, and per group and zone over
/// `names`.
fn tally(
    run: &Run,
    periods: &BTreeMap<String, String>,
    sources: &BTreeMap<String, String>,
    names: &BTreeSet<String>,
) -> Tally {
    let mut counts = BTreeMap::<String, Counts>::new();
    let mut zones = BTreeMap::<(String, String), Unit>::new();
    for (id, sample) in &run.samples {
        if !sample.counted {
            continue;
        }
        let non_het = sample
            .escapes
            .iter()
            .any(|escape| escape.class == Class::NonHet);
        let het_only = !non_het
            && sample
                .escapes
                .iter()
                .any(|escape| escape.class == Class::Het);
        for group in groups(id, periods, sources) {
            let count = counts.entry(group.clone()).or_default();
            count.samples += 1;
            count.positions += sample.positions;
            count.unreviewed_non_het += usize::from(non_het);
            count.unreviewed_het_only += usize::from(het_only);
            count.auto_pass.add(&sample.pass);
            for name in names {
                zones
                    .entry((group.clone(), name.clone()))
                    .or_insert_with(|| Unit::new(name.clone()))
                    .add(&sample.pass, name);
            }
        }
    }
    Tally { counts, zones }
}

/// Checks that the files count the passes the run's own `run.json` records:
/// over all samples, and over the current data in `current` when the run
/// records it, the counted samples, the auto-passed samples by their most
/// severe block, and every zone's units.
fn consistent(run: &Run, tally: &Tally, current: &str) -> Result<(), DiffError> {
    let mut found = Vec::new();
    let mut check = |group: &str, passes: &Passes| {
        let counted = tally.counts.get(group).copied().unwrap_or_default();
        if passes.samples != counted.samples {
            found.push(format!(
                "{group}: {} samples recorded, {} counted",
                passes.samples, counted.samples
            ));
        }
        if passes.auto_pass != counted.auto_pass {
            found.push(format!("{group}: the auto-pass samples"));
        }
        for unit in &passes.zones.by_zone {
            let key = (group.to_owned(), unit.name.clone());
            if tally.zones.get(&key) != Some(unit) {
                found.push(format!("{group}: zone {}", unit.name));
            }
        }
    };
    check("all", &run.header.workload.all);
    if let Some(passes) = &run.header.workload.current {
        check(current, passes);
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

/// Each sample's label in one split.
type Labels = BTreeMap<String, String>;

/// The label of the current period from `since`.
fn current(since: Month) -> String {
    format!("current (≥ {since})")
}

/// The group the gate reads: the current period when the run records a
/// current-data month, else all samples.
fn gated(run: &Run) -> String {
    run.header.current_since.map_or_else(
        || "all".to_owned(),
        |since| format!("period {}", current(since)),
    )
}

/// Each sample's period, by its manifest `batch_date` against the run's
/// current-data month, and its source, when the manifest has the column.
fn split(run: &Run) -> Result<(Labels, Labels), DiffError> {
    let path = &run.header.manifest.path;
    let bytes = std::fs::read(path).map_err(|error| failure(path, error))?;
    let manifest = manifest::parse(&bytes).map_err(|error| failure(path, error))?;
    if manifest.body_sha256 != run.header.manifest.body_sha256 {
        return Err(failure(
            path,
            format!(
                "its body has sha256 {}, the run recorded {}",
                manifest.body_sha256, run.header.manifest.body_sha256
            ),
        ));
    }
    let since = run.header.current_since;
    let source = manifest.columns.iter().position(|name| name == "source");
    let mut periods = BTreeMap::new();
    let mut sources = BTreeMap::new();
    for sample in &manifest.samples {
        if let (Some(since), Some(month)) = (since, sample.batch_month) {
            let period = if month >= since {
                current(since)
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
    let gate = gated(&base);
    let in_gate = |id: &str| groups(id, &periods, &sources).contains(&gate);
    let names = base
        .zone_names()
        .into_iter()
        .chain(run.zone_names())
        .collect::<BTreeSet<_>>();
    let before = tally(&base, &periods, &sources, &names);
    let after = tally(&run, &periods, &sources, &names);
    consistent(&base, &before, &gate)?;
    consistent(&run, &after, &gate)?;
    let mut report = String::from("# Run comparison\n\n");
    let revision = |run: &Run| {
        run.header
            .source
            .revision
            .as_deref()
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
    // A sample that failed on one side only is compared against nothing: the
    // two runs no longer cover the same samples, whether the failure is the
    // change's or a transient read error. Only the gated group's samples gate.
    let (failed_before, failed_after) = (base.failed(), run.failed());
    if failed_before != failed_after {
        let one_side = failed_before
            .symmetric_difference(&failed_after)
            .collect::<Vec<_>>();
        let differing = one_side.iter().filter(|id| in_gate(id)).count();
        let message = format!(
            "failed samples differ: {} in the base, {} in the run, {} in one only ({differing} in {gate}); rerun a side that failed transiently, or see the run's failures",
            failed_before.len(),
            failed_after.len(),
            one_side.len()
        );
        let _ = writeln!(report, "**WARNING:** {message}.\n");
        if differing > 0 {
            failures.push(message);
        }
    }

    let _ = writeln!(
        report,
        "## Samples\n\nUnreviewed counts samples with a non-het unreviewed block, and those whose only unreviewed blocks are het; passed samples are counted under their most severe block. Only the non-het counts of {gate} gate; the other groups are reported.\n\n| Group | Samples | Clear | Unreviewed non-het | Unreviewed het only | Passed non-het | Passed truth-N | Passed het only | Passed re-check only | Passed length-note only | Positions to review |\n|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"
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
        let (passed_old, passed_new) = (old.auto_pass.unreviewed, new.auto_pass.unreviewed);
        let cell = |old: usize, new: usize| {
            if old == new {
                old.to_string()
            } else {
                format!("{old} → {new}")
            }
        };
        let _ = writeln!(
            report,
            "| {name} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |",
            cell(old.samples, new.samples),
            cell(old.auto_pass.samples, new.auto_pass.samples),
            cell(old.unreviewed_non_het, new.unreviewed_non_het),
            cell(old.unreviewed_het_only, new.unreviewed_het_only),
            cell(passed_old.non_het, passed_new.non_het),
            cell(passed_old.truth_n_only, passed_new.truth_n_only),
            cell(passed_old.het_only, passed_new.het_only),
            cell(passed_old.recheck_only, passed_new.recheck_only),
            cell(passed_old.length_note_only, passed_new.length_note_only),
            if old.positions == new.positions {
                old.positions.to_string()
            } else {
                format!("{} → {}", old.positions, new.positions)
            },
        );
        if *name != gate {
            continue;
        }
        if passed_new.non_het > passed_old.non_het {
            failures.push(format!(
                "passed samples with a non-het block rise in {name}: {} → {}",
                passed_old.non_het, passed_new.non_het
            ));
        }
        if new.unreviewed_non_het > old.unreviewed_non_het {
            failures.push(format!(
                "samples with a non-het unreviewed block rise in {name}: {} → {}",
                old.unreviewed_non_het, new.unreviewed_non_het
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
        let counts = |tally: &Tally| {
            tally
                .zones
                .get(key)
                .map_or((0, 0), |unit| (unit.checked, unit.unreviewed.non_het))
        };
        let ((old_checked, old_non_het), (new_checked, new_non_het)) =
            (counts(&before), counts(&after));
        if (old_checked, old_non_het) == (new_checked, new_non_het) {
            continue;
        }
        let _ = writeln!(
            report,
            "| {} | {} | {old_checked} → {new_checked} | {old_non_het} → {new_non_het} |",
            key.0, key.1
        );
        if key.0 == gate && new_non_het > old_non_het {
            failures.push(format!(
                "passed {} zones with a non-het block rise in {}: {old_non_het} → {new_non_het}",
                key.1, key.0
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
        if gone.is_empty() && new.is_empty() && old.pass.clear == sample.pass.clear {
            continue;
        }
        let flag = |clear: bool| if clear { "Y" } else { "N" };
        let row = format!(
            "{id}\t{}\t{}\t{}\t{}\t{gone}\t{new}",
            periods.get(id).map_or("", String::as_str),
            sources.get(id).map_or("", String::as_str),
            flag(old.pass.clear),
            flag(sample.pass.clear),
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
    let _ = write!(report, "```\n\n## Verdict\n\nGated on {gate}.\n\n");
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

    /// Writes a run directory: its manifest (S1 from an earlier batch, S2
    /// from a current one), its `run.json` with the passes `workload`
    /// records, its scores, and its unreviewed blocks.
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
        let body = "S1\tb1\t2025-06-01\tclean\t1\t1\t0\t/t/a.ab1\nS2\tb2\t2026-02-01\tclean\t1\t1\t0\t/t/b.ab1\n";
        std::fs::write(
            &manifest,
            format!(
                "# body_sha256={:x}\nsample_id\tbatch\tbatch_date\tstratum\tweight\tweight_corpus\tanchor_candidate\ttrace_paths\n{body}",
                <sha2::Sha256 as sha2::Digest>::digest(body.as_bytes())
            ),
        )?;
        let counted = scores
            .iter()
            .filter(|row| {
                let lod = row.split('\t').nth(2).unwrap_or_default();
                lod == "concordant" || lod == "discordant"
            })
            .count();
        let mut header = rundir::tests::fixture();
        header.source.dirty = Some(dirty);
        header.manifest.path = manifest;
        header.manifest.body_sha256 = format!(
            "{:x}",
            <sha2::Sha256 as sha2::Digest>::digest(body.as_bytes())
        );
        header.truth.digest = "t".to_owned();
        header.recheck = None;
        header.current_since = None;
        header.review_window = 10;
        header.workload.current = None;
        let all = &mut header.workload.all;
        all.samples = counted;
        all.zones.by_zone.clear();
        all.auto_pass = AutoPass {
            samples: workload.0,
            unreviewed: Escaped {
                total: workload.1 + workload.2,
                non_het: workload.1,
                het_only: workload.2,
                ..Escaped::default()
            },
        };
        std::fs::write(directory.join("run.json"), header.text()?)?;
        // The scores give the primary status, `haplotype_lod`'s, and the
        // review columns; the other measures take the primary status.
        let mut text = String::from("sample_id");
        for measure in Measure::ALL {
            text.push('\t');
            text.push_str(measure.column());
        }
        text.push_str("\treview_clear\treview_positions\treview_zones\tuncovered\n");
        for row in scores {
            let cells = row.split('\t').collect::<Vec<_>>();
            text.push_str(cells[0]);
            for measure in Measure::ALL {
                text.push('\t');
                text.push_str(if measure == Measure::HaplotypeLod {
                    cells[2]
                } else {
                    cells[1]
                });
            }
            for cell in &cells[3..] {
                text.push('\t');
                text.push_str(cell);
            }
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

    /// Records `since` as the run's current-data month.
    fn dated(directory: &Path, since: &str) -> Result<(), Box<dyn std::error::Error>> {
        let mut header = RunHeader::read(directory)?;
        header.current_since = Some(since.parse()?);
        std::fs::write(directory.join("run.json"), header.text()?)?;
        Ok(())
    }

    // With a current-data month, only the current period gates: a non-het
    // rise or a lost sample in the earlier period is reported and passes,
    // a non-het rise in the current period fails.
    #[test]
    fn gates_the_current_period_only() -> Result<(), Box<dyn std::error::Error>> {
        let root = tempfile::tempdir()?;
        let clear = [
            "S1\tconcordant\tconcordant\tY\t0\t\t",
            "S2\tconcordant\tconcordant\tY\t0\t\t",
        ];
        let base = run(root.path(), "base", false, &clear, &[], (2, 0, 0))?;
        dated(&base, "2026-01")?;
        let earlier = run(
            root.path(),
            "earlier",
            false,
            &[
                "S1\tdiscordant\tdiscordant\tY\t0\t\t",
                "S2\tconcordant\tconcordant\tY\t0\t\t",
            ],
            &["S1\t309\t309\t309-309+:CC/C\tnon_het\twindow@HV2\tno"],
            (2, 1, 0),
        )?;
        dated(&earlier, "2026-01")?;
        let comparison = compare(&base, &earlier, false)?;
        assert!(!comparison.gate_failed);
        assert!(
            comparison
                .report
                .contains("Gated on period current (≥ 2026-01).\n\nGATE PASS")
        );
        assert!(
            comparison
                .report
                .contains("| period earlier (< 2026-01) | 1 | 1 | 0 → 1 | 0 | 0 → 1 |")
        );
        assert!(
            comparison
                .report
                .contains("| all | window@HV2 | 0 → 0 | 0 → 1 |")
        );
        let current = run(
            root.path(),
            "current",
            false,
            &[
                "S1\tconcordant\tconcordant\tY\t0\t\t",
                "S2\tdiscordant\tdiscordant\tY\t0\t\t",
            ],
            &["S2\t309\t309\t309-309+:CC/C\tnon_het\twindow@HV2\tno"],
            (2, 1, 0),
        )?;
        dated(&current, "2026-01")?;
        let comparison = compare(&base, &current, false)?;
        assert!(comparison.gate_failed);
        assert!(comparison.report.contains(
            "- passed samples with a non-het block rise in period current (≥ 2026-01): 0 → 1"
        ));
        assert!(!comparison.report.contains("rise in all"));
        let lost = run(
            root.path(),
            "lost",
            false,
            &[
                "S1\tfailed\tfailed\tN\t0\t\t",
                "S2\tconcordant\tconcordant\tY\t0\t\t",
            ],
            &[],
            (1, 0, 0),
        )?;
        dated(&lost, "2026-01")?;
        let comparison = compare(&base, &lost, false)?;
        assert!(!comparison.gate_failed);
        assert!(comparison.report.contains(
            "failed samples differ: 0 in the base, 1 in the run, 1 in one only (0 in period current"
        ));
        Ok(())
    }

    // Without a current-data month, every sample gates.
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
                .contains("| all | 2 | 2 → 1 | 0 → 1 | 1 | 0 → 1 |")
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

    // A het block the candidate leaves unreviewed is reported, not gated:
    // heteroplasmy is optional.
    #[test]
    fn reports_a_new_het_escape_without_gating() -> Result<(), Box<dyn std::error::Error>> {
        let root = tempfile::tempdir()?;
        let base = run(
            root.path(),
            "base",
            false,
            &["S1\tconcordant\tconcordant\tY\t0\t\t"],
            &[],
            (1, 0, 0),
        )?;
        let candidate = run(
            root.path(),
            "run",
            false,
            &["S1\tdiscordant\tdiscordant\tY\t0\t\t"],
            &["S1\t461\t461\t461-461:C/S\thet\twindow@HV2\tno"],
            (1, 0, 1),
        )?;
        let comparison = compare(&base, &candidate, false)?;
        assert!(!comparison.gate_failed);
        assert!(comparison.report.contains("GATE PASS"));
        assert!(
            comparison
                .report
                .contains("| all | 1 | 1 | 0 | 0 → 1 | 0 | 0 | 0 → 1 |")
        );
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
        // Another profile, or a manifest edited since the run, is refused.
        let profile = run(root.path(), "profile", false, &rows, &[], (1, 0, 0))?;
        let mut header = RunHeader::read(&profile)?;
        header.identities[0].profile_sha256 = "another".to_owned();
        std::fs::write(profile.join(rundir::RUN_JSON), header.text()?)?;
        assert!(matches!(
            compare(&base, &profile, false),
            Err(DiffError::Mismatch(message)) if message.starts_with("profile differs")
        ));
        let edited = run(root.path(), "edited", false, &rows, &[], (1, 0, 0))?;
        let mut header = RunHeader::read(&edited)?;
        header.manifest.body_sha256 = "c".repeat(64);
        std::fs::write(edited.join(rundir::RUN_JSON), header.text()?)?;
        assert!(matches!(
            compare(&edited, &edited, false),
            Err(DiffError::Run { message, .. }) if message.contains("the run recorded cccc")
        ));
        let miscounted = run(root.path(), "miscounted", false, &rows, &[], (0, 0, 0))?;
        assert!(matches!(
            compare(&base, &miscounted, false),
            Err(DiffError::Inconsistent(..))
        ));
        Ok(())
    }
    #[test]
    fn gates_samples_that_failed_on_one_side_only() -> Result<(), Box<dyn std::error::Error>> {
        let root = tempfile::tempdir()?;
        let scored = [
            "S1\tconcordant\tconcordant\tY\t0\t\t",
            "S2\tconcordant\tconcordant\tY\t0\t\t",
        ];
        let base = run(root.path(), "base", false, &scored, &[], (2, 0, 0))?;
        let lost = [
            "S1\tconcordant\tconcordant\tY\t0\t\t",
            "S2\tfailed\tfailed\tN\t0\t\t",
        ];
        let flaky = run(root.path(), "flaky", false, &lost, &[], (1, 0, 0))?;
        let comparison = compare(&base, &flaky, false)?;
        assert!(comparison.gate_failed);
        assert!(
            comparison
                .report
                .contains("failed samples differ: 0 in the base, 1 in the run")
        );
        assert!(!compare(&flaky, &flaky, false)?.gate_failed);
        Ok(())
    }

    /// A run written by the harness's own writers reads back through
    /// `diff`: its files count the passes its `run.json` records, over all
    /// samples and the current data, and the run compared with itself passes.
    #[test]
    fn a_written_run_reads_back_consistently() -> Result<(), Box<dyn std::error::Error>> {
        use crate::cause::forms::Forms;
        use crate::compare::Allele;
        use crate::report::Report;
        use crate::runner::{Called, Outcome};
        use crate::truth::{TruthCheck, TruthDatabase, TruthSample};

        let root = tempfile::tempdir()?;
        // Four samples, two earlier and two current; each a different case.
        let body = "E1\tb1\t2025-06-01\tclean\t1\t1\t0\t/t/a.ab1\n\
                    E2\tb1\t2025-06-01\tclean\t1\t1\t0\t/t/b.ab1\n\
                    C1\tb2\t2026-02-01\tclean\t1\t1\t0\t/t/c.ab1\n\
                    C2\tb2\t2026-02-01\tclean\t1\t1\t0\t/t/d.ab1\n";
        let path = root.path().join("panel.tsv");
        std::fs::write(
            &path,
            format!(
                "# panel=dev\n# body_sha256={:x}\nsample_id\tbatch\tbatch_date\tstratum\tweight\tweight_corpus\tanchor_candidate\ttrace_paths\n{body}",
                <sha2::Sha256 as sha2::Digest>::digest(body.as_bytes())
            ),
        )?;
        let manifest = manifest::parse(&std::fs::read(&path)?)?;
        let truth = manifest
            .samples
            .iter()
            .map(|sample| {
                (
                    sample.id.clone(),
                    TruthSample {
                        batch: sample.batch.clone(),
                        status: None,
                        alleles: vec![Allele::new("263", "A", "G")],
                        intervals: vec![(73, 340)],
                    },
                )
            })
            .collect::<BTreeMap<_, _>>();
        // E1 and C1 call a DNA-only 150T, unchecked; E2 calls it under a
        // check item; C2 agrees with the truth.
        let outcomes = manifest
            .samples
            .iter()
            .map(|sample| {
                let mut alleles = vec![Allele::new("263", "A", "G")];
                if sample.id != "C2" {
                    alleles.push(Allele::new("150", "C", "T"));
                }
                let checks = if sample.id == "E2" {
                    vec![(145, 155)]
                } else {
                    Vec::new()
                };
                Outcome {
                    traces: Some(1),
                    result: Ok(Called {
                        rejected: 0,
                        alleles,
                        not_recorded: BTreeSet::new(),
                        intervals: vec![(73, 340)],
                        signals: crate::lod::Signals::default(),
                        loci: crate::cause::loci::Loci::default(),
                        identity: rundir::tests::fixture().identities[0].clone(),
                        empop: None,
                        check_labels: checks
                            .iter()
                            .map(|&(first, last)| (first, last, "x".to_owned()))
                            .collect(),
                        checks,
                        notes: Vec::new(),
                        components: Vec::new(),
                        tracts: Vec::new(),
                        length_context: None,
                    }),
                    millis: 0,
                }
            })
            .collect::<Vec<_>>();
        let rows = crate::report::rows(
            &manifest.samples,
            &outcomes,
            &truth,
            crate::compare::tests::rcrs(),
            (crate::geometry::tests::rcrs(), &[]),
        );
        let database = TruthDatabase {
            info: BTreeMap::new(),
            samples: truth.clone(),
        };
        let check = TruthCheck::new(&manifest, &database);
        let zones = crate::zones::Zones::new(
            &crate::geometry::tests::spans(&[(73, 340)]),
            &[(303, 315)],
            crate::zones::Widths {
                after: 20,
                edge: 40,
                overrun: 2,
            },
        );
        let report = Report::new(
            &rows,
            (crate::geometry::tests::rcrs(), &check),
            (&Forms::default(), ""),
        )
        .with_review(10, None, Some(zones), None, Vec::new())
        .with_current_since("2026-01".parse().ok());
        let directory = root.path().join("run");
        std::fs::create_dir_all(&directory)?;
        for (name, contents) in report.files(&[]) {
            rundir::write(&directory, name, &contents)?;
        }
        let (all, current) = report.passes();
        let mut header = rundir::tests::fixture();
        header.manifest.path = path;
        header
            .manifest
            .body_sha256
            .clone_from(&manifest.body_sha256);
        header.recheck = None;
        header.current_since = "2026-01".parse().ok();
        header.workload = rundir::WorkloadHeader { all, current };
        std::fs::write(directory.join(rundir::RUN_JSON), header.text()?)?;

        let comparison = compare(&directory, &directory, false)?;
        assert!(!comparison.gate_failed, "{}", comparison.report);
        assert!(
            comparison
                .report
                .contains("| period current (≥ 2026-01) | 2 | 2 | 1 | 0 | 1 | 0 | 0 | 0 | 0 | 0 |"),
            "{}",
            comparison.report
        );
        // A run.json that records other passes is refused.
        let mut miscounted = header.clone();
        if let Some(current) = &mut miscounted.workload.current {
            current.zones.by_zone[0].checked += 1;
        }
        std::fs::write(directory.join(rundir::RUN_JSON), miscounted.text()?)?;
        assert!(matches!(
            compare(&directory, &directory, false),
            Err(DiffError::Inconsistent(..))
        ));
        Ok(())
    }

    #[test]
    fn refuses_a_sample_scored_twice_or_a_missing_column() -> Result<(), Box<dyn std::error::Error>>
    {
        let root = tempfile::tempdir()?;
        let rows = [
            "S1\tconcordant\tconcordant\tY\t0\t\t",
            "S1\tconcordant\tconcordant\tY\t0\t\t",
        ];
        let twice = run(root.path(), "twice", false, &rows, &[], (2, 0, 0))?;
        assert!(matches!(
            compare(&twice, &twice, false),
            Err(DiffError::Run { message, .. }) if message.contains("sample S1 appears twice")
        ));
        let once = run(root.path(), "once", false, &rows[..1], &[], (1, 0, 0))?;
        let scores = once.join(rundir::SCORES);
        let text = std::fs::read_to_string(&scores)?.replace("review_zones", "zones");
        std::fs::write(&scores, text)?;
        assert!(matches!(
            compare(&once, &once, false),
            Err(DiffError::Run { message, .. }) if message == "no column review_zones"
        ));
        Ok(())
    }
}
