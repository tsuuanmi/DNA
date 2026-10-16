//! The run's outputs: per-sample comparisons, estimates, anchors, the diff
//! against a baseline, and a Markdown summary.
//!
//! Every sample is compared under each [`Measure`]. Every output but the run
//! header (`run.json`, which records timings) is a function of the panel, the
//! truth, and DNA's calls only, so two runs of the same inputs give
//! byte-identical files whatever the thread count.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

use crate::baseline::{Class, Statuses};
use crate::compare::{self, Comparison, Concordance, Failure, Measure, Reason, Region, Side};
use crate::estimate::{self, Estimate, Observation};
use crate::haplotype;
use crate::manifest::Sample;
use crate::runner::Outcome;
use crate::truth::{TruthCheck, TruthSample, TruthState};

/// One compared panel sample.
pub(crate) struct Row<'a> {
    pub(crate) sample: &'a Sample,
    pub(crate) outcome: &'a Outcome,
    /// The sample's truth record, when the truth has one.
    pub(crate) truth: Option<&'a TruthSample>,
    /// One comparison per measure, in [`Measure::ALL`] order.
    pub(crate) comparisons: Result<[Comparison; Measure::COUNT], Failure>,
}

/// Compares every sample's outcome, given in sample order, with its truth
/// under every measure, on `reference` (the reference sequence DNA ran on),
/// and sorts the rows by sample identifier.
pub(crate) fn rows<'a>(
    samples: &'a [Sample],
    outcomes: &'a [Outcome],
    truth: &'a BTreeMap<String, TruthSample>,
    reference: &str,
) -> Vec<Row<'a>> {
    let mut rows = samples
        .iter()
        .zip(outcomes)
        .map(|(sample, outcome)| {
            let record = truth.get(&sample.id);
            Row {
                sample,
                outcome,
                truth: record,
                comparisons: match (&outcome.result, record) {
                    (Err(failure), _) => Err(failure.clone()),
                    (Ok(_), None) => Err(Failure::TruthMissing),
                    (Ok(called), Some(record)) => {
                        let dna = Side {
                            alleles: &called.alleles,
                            intervals: &called.intervals,
                        };
                        let truth = Side {
                            alleles: &record.alleles,
                            intervals: &record.intervals,
                        };
                        let haplotypes = compare::haplotypes(dna, truth, reference);
                        Ok(Measure::ALL
                            .map(|measure| compare::compare(dna, truth, measure, &haplotypes)))
                    }
                },
            }
        })
        .collect::<Vec<_>>();
    rows.sort_by(|left, right| left.sample.id.cmp(&right.sample.id));
    rows
}

/// Outcome of one sample under one measure.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Status {
    /// Concordance `Y`.
    Concordant,
    /// Concordance `N`.
    Discordant,
    /// Concordance `N/A`: one side has no allele.
    NotApplicable,
    /// DNA or the truth gave nothing to compare.
    Failed,
}

impl Status {
    const ALL: [Self; 4] = [
        Self::Concordant,
        Self::Discordant,
        Self::NotApplicable,
        Self::Failed,
    ];

    pub(crate) const fn label(self) -> &'static str {
        match self {
            Self::Concordant => "concordant",
            Self::Discordant => "discordant",
            Self::NotApplicable => "not_applicable",
            Self::Failed => "failed",
        }
    }

    /// The concordance code of the lab's comparison.
    const fn code(self) -> &'static str {
        match self {
            Self::Concordant => "Y",
            Self::Discordant => "N",
            Self::NotApplicable => "N/A",
            Self::Failed => "",
        }
    }
}

impl Row<'_> {
    /// The comparison under `measure`, or the failure.
    pub(crate) fn comparison(&self, measure: Measure) -> Result<&Comparison, &Failure> {
        self.comparisons
            .as_ref()
            .map(|comparisons| &comparisons[measure.index()])
    }

    pub(crate) fn status(&self, measure: Measure) -> Status {
        match self.comparison(measure) {
            Ok(comparison) => match comparison.concordance {
                Concordance::Concordant => Status::Concordant,
                Concordance::Discordant => Status::Discordant,
                Concordance::NotApplicable => Status::NotApplicable,
            },
            Err(_) => Status::Failed,
        }
    }
}

/// How samples are weighted for an estimate. Retired samples weigh nothing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Weighting {
    /// `weight_corpus`: the panel scaled to the whole corpus.
    Corpus,
    /// `weight`: the panel scaled to its sampling pool.
    Pool,
}

impl Weighting {
    const ALL: [Self; 2] = [Self::Corpus, Self::Pool];

    const fn label(self) -> &'static str {
        match self {
            Self::Corpus => "corpus",
            Self::Pool => "pool",
        }
    }

    pub(crate) fn of(self, sample: &Sample) -> f64 {
        if sample.retired {
            return 0.0;
        }
        match self {
            Self::Corpus => sample.weight_corpus.value,
            Self::Pool => sample.weight.value,
        }
    }
}

/// The event an estimate counts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Metric {
    /// Anything but `Y` under the measure: discordant, not applicable, or
    /// failed.
    NotConcordant(Measure),
    /// Concordance `N` under the measure.
    Discordant(Measure),
    /// The sample failed, whatever the measure.
    Failed,
}

impl Metric {
    fn all() -> Vec<Self> {
        Measure::ALL
            .into_iter()
            .flat_map(|measure| [Self::NotConcordant(measure), Self::Discordant(measure)])
            .chain([Self::Failed])
            .collect()
    }

    const fn measure(self) -> &'static str {
        match self {
            Self::NotConcordant(measure) | Self::Discordant(measure) => measure.label(),
            Self::Failed => "all",
        }
    }

    const fn label(self) -> &'static str {
        match self {
            Self::NotConcordant(_) => "not_concordant",
            Self::Discordant(_) => "discordant",
            Self::Failed => "failed",
        }
    }

    fn of(self, row: &Row<'_>) -> bool {
        match self {
            Self::NotConcordant(measure) => row.status(measure) != Status::Concordant,
            Self::Discordant(measure) => row.status(measure) == Status::Discordant,
            Self::Failed => row.comparisons.is_err(),
        }
    }
}

/// How a sample's outcome moved against the baseline under one measure.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Change {
    /// Exact or concordant in the baseline, not now.
    Regression,
    /// Not exact or concordant in the baseline, exact or concordant now.
    Fix,
    /// Another change of status class.
    Changed,
    /// Absent from the baseline.
    Added,
    /// In the baseline, absent from this panel.
    Removed,
    /// The same status class.
    Same,
}

impl Change {
    const ALL: [Self; 5] = [
        Self::Regression,
        Self::Fix,
        Self::Changed,
        Self::Added,
        Self::Removed,
    ];

    const fn label(self) -> &'static str {
        match self {
            Self::Regression => "regression",
            Self::Fix => "fix",
            Self::Changed => "changed",
            Self::Added => "added",
            Self::Removed => "removed",
            Self::Same => "same",
        }
    }
}

/// Every derived output of one run.
pub(crate) struct Report<'a> {
    rows: &'a [Row<'a>],
    baseline: Option<&'a Statuses>,
    check: &'a TruthCheck,
    estimates: Vec<(Weighting, Metric, Estimate)>,
}

impl<'a> Report<'a> {
    /// Derives the estimates of `rows`, which are sorted by sample identifier.
    pub(crate) fn new(
        rows: &'a [Row<'a>],
        baseline: Option<&'a Statuses>,
        check: &'a TruthCheck,
    ) -> Self {
        let estimates = Weighting::ALL
            .into_iter()
            .flat_map(|weighting| {
                Metric::all()
                    .into_iter()
                    .map(move |metric| (weighting, metric))
            })
            .filter_map(|(weighting, metric)| {
                let observations = rows
                    .iter()
                    .map(|row| Observation {
                        stratum: &row.sample.stratum,
                        weight: weighting.of(row.sample),
                        event: metric.of(row),
                    })
                    .collect::<Vec<_>>();
                estimate::estimate(&observations).map(|estimate| (weighting, metric, estimate))
            })
            .collect();
        Self {
            rows,
            baseline,
            check,
            estimates,
        }
    }

    /// The change under `measure`, or `None` without a baseline or when the
    /// baseline does not report the measure.
    fn change(&self, row: &Row<'_>, measure: Measure) -> Option<Change> {
        let baseline = self.baseline.filter(|baseline| baseline.reports(measure))?;
        let Some(before) = baseline.get(&row.sample.id) else {
            return Some(Change::Added);
        };
        let before = Class::of(before.get(&measure)?);
        let now = Class::of(row.status(measure).label());
        Some(if before == Class::Good && now != Class::Good {
            Change::Regression
        } else if before != Class::Good && now == Class::Good {
            Change::Fix
        } else if before != now {
            Change::Changed
        } else {
            Change::Same
        })
    }

    /// Baseline samples absent from this panel, with their statuses.
    fn removed(&self) -> Vec<(&'a str, &'a BTreeMap<Measure, String>)> {
        let present = self
            .rows
            .iter()
            .map(|row| row.sample.id.as_str())
            .collect::<BTreeSet<_>>();
        self.baseline
            .into_iter()
            .flat_map(|baseline| &baseline.samples)
            .filter(|(sample, _)| !present.contains(sample.as_str()))
            .map(|(sample, statuses)| (sample.as_str(), statuses))
            .collect()
    }

    /// Every deterministic output file, by name; `diff.tsv` only with a
    /// baseline, `truth_changes.tsv` whenever the panel's truth changed.
    /// `context` heads the summary.
    pub(crate) fn files(&self, context: &[(&str, String)]) -> Vec<(&'static str, String)> {
        let mut files = vec![
            ("scores.tsv", self.scores()),
            ("estimates.tsv", self.estimates()),
            ("anchors.tsv", self.anchors()),
        ];
        if let Some(diff) = self.diff() {
            files.push(("diff.tsv", diff));
        }
        if self.check.state() == TruthState::Changed {
            files.push(("truth_changes.tsv", self.check.changes()));
        }
        files.push(("summary.md", self.summary(context)));
        files
    }

    /// `scores.tsv`: one row per sample.
    fn scores(&self) -> String {
        let mut text = String::from(
            "sample_id\tbatch\ttruth_batch\tstratum\tweight\tweight_corpus\tanchor_candidate\t\
             retired\tstatus\tconcordant\tdna_alleles\ttruth_alleles\tdna_unique\ttruth_unique\t\
             reasons\tpositions\tiupac_only\tregions\tstrict_status\tstrict_concordant\t\
             uncovered\tiupac_compatible_status\tiupac_compatible_concordant\t\
             haplotype_status\thaplotype_concordant\thaplotype_blocks\t\
             haplotype_iupac_compatible_status\thaplotype_iupac_compatible_concordant\t\
             haplotype_iupac_compatible_blocks\ttraces\trejected_traces\tfailure\n",
        );
        for row in self.rows {
            let sample = row.sample;
            let primary = row.status(Measure::Primary);
            let strict = row.status(Measure::Strict);
            let compatible = row.status(Measure::IupacCompatible);
            let haplotype = [Measure::Haplotype, Measure::HaplotypeIupacCompatible]
                .map(|measure| {
                    let status = row.status(measure);
                    let blocks = row
                        .comparison(measure)
                        .map(|comparison| {
                            comparison
                                .blocks
                                .iter()
                                .map(haplotype::Block::text)
                                .collect::<Vec<_>>()
                                .join(" ")
                        })
                        .unwrap_or_default();
                    format!("{}\t{}\t{blocks}", status.label(), status.code())
                })
                .join("\t");
            let detail = row.comparison(Measure::Primary).map_or_else(
                |_| "\t\t\t\t\t\t\t".to_owned(),
                |comparison| {
                    format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        comparison.dna_alleles,
                        comparison.truth_alleles,
                        compare::texts(&comparison.dna_unique),
                        compare::texts(&comparison.truth_unique),
                        labels(&comparison.reasons),
                        compare::sorted_positions(&comparison.positions).join(" "),
                        u8::from(comparison.iupac_only),
                        regions(comparison)
                            .iter()
                            .map(|region| region.label())
                            .collect::<Vec<_>>()
                            .join(","),
                    )
                },
            );
            let uncovered = row
                .comparison(Measure::Primary)
                .map(|comparison| compare::texts(&comparison.uncovered))
                .unwrap_or_default();
            let traces = row
                .outcome
                .traces
                .map(|traces| traces.to_string())
                .unwrap_or_default();
            let rejected = row
                .outcome
                .result
                .as_ref()
                .map(|called| called.rejected.to_string())
                .unwrap_or_default();
            let failure = row
                .comparisons
                .as_ref()
                .err()
                .map_or_else(String::new, |failure| {
                    clean(&format!("{}: {failure}", failure.kind()))
                });
            let _ = writeln!(
                text,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{detail}\t{}\t{}\t{uncovered}\t{}\t{}\t\
                 {haplotype}\t{traces}\t{rejected}\t{failure}",
                sample.id,
                sample.batch,
                row.truth.map_or("", |truth| truth.batch.as_str()),
                sample.stratum,
                sample.weight.text,
                sample.weight_corpus.text,
                u8::from(sample.anchor),
                u8::from(sample.retired),
                primary.label(),
                primary.code(),
                strict.label(),
                strict.code(),
                compatible.label(),
                compatible.code(),
            );
        }
        text
    }

    /// `estimates.tsv`: per weighting, measure, and metric, each stratum then
    /// `all`.
    fn estimates(&self) -> String {
        let mut text = String::from(
            "weighting\tmeasure\tmetric\tstratum\tsamples\tpopulation\tevents\trate\tci95\tunweighted\n",
        );
        for (weighting, metric, estimate) in &self.estimates {
            let prefix = format!(
                "{}\t{}\t{}",
                weighting.label(),
                metric.measure(),
                metric.label()
            );
            for stratum in &estimate.strata {
                let _ = writeln!(
                    text,
                    "{prefix}\t{}\t{}\t{:.2}\t{}\t{:.4}\t\t",
                    stratum.name, stratum.samples, stratum.population, stratum.events, stratum.rate
                );
            }
            let _ = writeln!(
                text,
                "{prefix}\tall\t{}\t{:.2}\t{}\t{:.4}\t{:.4}\t{:.4}",
                estimate.samples,
                estimate.population,
                estimate.events,
                estimate.rate,
                estimate.ci95,
                estimate.unweighted
            );
        }
        text
    }

    /// `anchors.tsv`: every anchor candidate under each measure, with its
    /// baseline status and primary change.
    fn anchors(&self) -> String {
        let mut text = String::from("sample_id\tstratum");
        for measure in Measure::ALL {
            let _ = write!(text, "\t{}", measure.column());
        }
        text.push_str("\tbaseline_status\tchange\n");
        for row in self.rows.iter().filter(|row| row.sample.anchor) {
            let before = self
                .baseline
                .and_then(|baseline| baseline.get(&row.sample.id))
                .and_then(|statuses| statuses.get(&Measure::Primary))
                .map_or("", String::as_str);
            let change = self.change(row, Measure::Primary).map_or("", Change::label);
            let _ = write!(text, "{}\t{}", row.sample.id, row.sample.stratum);
            for measure in Measure::ALL {
                let _ = write!(text, "\t{}", row.status(measure).label());
            }
            let _ = writeln!(text, "\t{before}\t{change}");
        }
        text
    }

    /// `diff.tsv`: per measure the baseline reports, the samples whose status
    /// class changed.
    fn diff(&self) -> Option<String> {
        let baseline = self.baseline?;
        let mut lines = BTreeMap::new();
        for row in self.rows {
            for measure in Measure::ALL {
                let Some(change) = self
                    .change(row, measure)
                    .filter(|change| *change != Change::Same)
                else {
                    continue;
                };
                let before = baseline
                    .get(&row.sample.id)
                    .and_then(|statuses| statuses.get(&measure))
                    .map_or("", String::as_str);
                lines.insert(
                    (row.sample.id.as_str(), measure),
                    format!(
                        "{}\t{}\t{}\t{before}\t{}\t{}",
                        row.sample.id,
                        row.sample.stratum,
                        measure.label(),
                        row.status(measure).label(),
                        change.label()
                    ),
                );
            }
        }
        for (sample, statuses) in self.removed() {
            for (measure, before) in statuses {
                lines.insert(
                    (sample, *measure),
                    format!(
                        "{sample}\t\t{}\t{before}\t\t{}",
                        measure.label(),
                        Change::Removed.label()
                    ),
                );
            }
        }
        let mut text =
            String::from("sample_id\tstratum\tmeasure\tbaseline_status\tstatus\tchange\n");
        for line in lines.into_values() {
            text.push_str(&line);
            text.push('\n');
        }
        Some(text)
    }

    /// `summary.md`: the run in aggregate.
    fn summary(&self, context: &[(&str, String)]) -> String {
        let mut text = String::from("# Evaluation run\n\n");
        if self.check.state() == TruthState::Changed {
            if self.check.changed.is_empty() {
                text.push_str(
                    "**Truth changed:** the panel truth digest differs from the panel's \
                     record although no sample's recorded digest does, so the manifest's \
                     header and rows disagree (`truth_changes.tsv` lists no sample); derive a \
                     new panel version.\n\n",
                );
            } else {
                let _ = writeln!(
                    text,
                    "**Truth changed:** the truth of {} panel samples differs from the \
                     panel's record (`truth_changes.tsv`); derive a new panel version.\n",
                    self.check.changed.len()
                );
            }
        }
        text.push_str("| Item | Value |\n|---|---|\n");
        for (item, value) in context {
            let _ = writeln!(text, "| {item} | {value} |");
        }
        self.panel_section(&mut text);
        self.outcome_section(&mut text);
        self.estimate_section(&mut text);
        self.discordance_section(&mut text);
        self.secondary_section(&mut text);
        self.anchor_section(&mut text);
        self.baseline_section(&mut text);

        let mut failures = BTreeMap::<&str, usize>::new();
        for row in self.rows {
            if let Err(failure) = &row.comparisons {
                *failures.entry(failure.kind()).or_default() += 1;
            }
        }
        if !failures.is_empty() {
            text.push_str("\n## Failures\n\n| Kind | Samples |\n|---|---:|\n");
            for (kind, count) in failures {
                let _ = writeln!(text, "| {kind} | {count} |");
            }
        }
        text
    }

    fn panel_section(&self, text: &mut String) {
        let count =
            |keep: &dyn Fn(&Row<'_>) -> bool| self.rows.iter().filter(|row| keep(row)).count();
        let _ = write!(
            text,
            "\n## Panel against the truth\n\nPanel truth: {} (digest {}, recorded {}).\n\n\
             | Check | Samples |\n|---|---:|\n\
             | panel samples | {} |\n| truth changed since the panel was built | {} |\n\
             | retired (left out of the estimates) | {} |\n\
             | without a truth record | {} |\n\
             | truth batch differs from the panel batch | {} |\n\
             | truth status not `OK` in its batch | {} |\n",
            self.check.state().label(),
            self.check.current,
            self.check.recorded.as_deref().unwrap_or("-"),
            self.rows.len(),
            self.check.changed.len(),
            count(&|row| row.sample.retired),
            count(&|row| row.truth.is_none()),
            count(&|row| row
                .truth
                .is_some_and(|truth| truth.batch != row.sample.batch)),
            count(&|row| {
                row.truth
                    .is_some_and(|truth| truth.status.as_deref() != Some("OK"))
            }),
        );
    }

    fn outcome_section(&self, text: &mut String) {
        text.push_str(
            "\n## Outcomes\n\nThe primary measure is the lab's concordance, DNA as tool A and the \
             truth as tool B. The strict measure also counts truth alleles where DNA has no \
             callable base; the IUPAC-compatible measure lets DNA's base (its call, or the \
             reference where it calls nothing) match a truth IUPAC code that admits it. The \
             haplotype measures apply the primary and IUPAC-compatible rules only to alleles \
             that do not describe the same sequence as the other side's. All panel samples, \
             retired and zero-weight samples included.\n\n| Status | Code |",
        );
        table_header(text, "|---|---|");
        for status in Status::ALL {
            let _ = write!(text, "| {} | {} |", status.label(), status.code());
            for measure in Measure::ALL {
                let count = self
                    .rows
                    .iter()
                    .filter(|row| row.status(measure) == status)
                    .count();
                let _ = write!(text, " {count} |");
            }
            text.push('\n');
        }
    }

    fn estimate_section(&self, text: &mut String) {
        text.push_str(
            "\n## Estimates\n\nRate = Σ w·y / Σ w with the 95 % interval of the stratified \
             variance; retired and zero-weight samples are left out.\n\n\
             | Weighting | Measure | Metric | Samples | Rate | 95 % CI | Unweighted |\n\
             |---|---|---|---:|---:|---:|---:|\n",
        );
        for (weighting, metric, estimate) in &self.estimates {
            let _ = writeln!(
                text,
                "| {} | {} | {} | {} | {:.4} | ± {:.4} | {:.4} |",
                weighting.label(),
                metric.measure(),
                metric.label(),
                estimate.samples,
                estimate.rate,
                estimate.ci95,
                estimate.unweighted
            );
        }
        let corpus = |measure: Measure| {
            self.estimates
                .iter()
                .find(|(weighting, metric, _)| {
                    *weighting == Weighting::Corpus && *metric == Metric::NotConcordant(measure)
                })
                .map(|(_, _, estimate)| estimate)
        };
        let Some(estimates) = Measure::ALL
            .into_iter()
            .map(corpus)
            .collect::<Option<Vec<_>>>()
        else {
            return;
        };
        text.push_str(
            "\n### Not concordant by stratum (corpus weighting)\n\n\
             | Stratum | Population | Samples |",
        );
        table_header(text, "|---|---:|---:|");
        for stratum in &estimates[0].strata {
            let _ = write!(
                text,
                "| {} | {:.0} | {} |",
                stratum.name, stratum.population, stratum.samples
            );
            for estimate in &estimates {
                let rate = estimate
                    .strata
                    .iter()
                    .find(|other| other.name == stratum.name)
                    .map_or(0.0, |other| other.rate);
                let _ = write!(text, " {rate:.3} |");
            }
            text.push('\n');
        }
    }

    fn discordance_section(&self, text: &mut String) {
        let discordant = self
            .rows
            .iter()
            .filter_map(|row| row.comparison(Measure::Primary).ok())
            .filter(|comparison| comparison.concordance == Concordance::Discordant)
            .collect::<Vec<_>>();
        let _ = write!(
            text,
            "\n## Discordance (primary)\n\n{} discordant samples; {} of them differ only at \
             positions where the truth has an IUPAC allele.\n\n| Rule | Samples |\n|---|---:|\n",
            discordant.len(),
            discordant
                .iter()
                .filter(|comparison| comparison.iupac_only)
                .count()
        );
        for reason in Reason::ALL
            .into_iter()
            .filter(|reason| *reason != Reason::Uncovered)
        {
            let count = discordant
                .iter()
                .filter(|comparison| comparison.reasons.contains(&reason))
                .count();
            let _ = writeln!(text, "| {} | {count} |", reason.label());
        }
        text.push_str(
            "\nBy the window of a discordant position:\n\n| Region | Samples |\n|---|---:|\n",
        );
        for region in Region::ALL {
            let count = discordant
                .iter()
                .filter(|comparison| regions(comparison).contains(&region))
                .count();
            let _ = writeln!(text, "| {} | {count} |", region.label());
        }
    }

    fn secondary_section(&self, text: &mut String) {
        let compared = self
            .rows
            .iter()
            .filter(|row| row.comparisons.is_ok())
            .collect::<Vec<_>>();
        let uncovered = compared
            .iter()
            .filter_map(|row| row.comparison(Measure::Primary).ok())
            .map(|comparison| comparison.uncovered.len())
            .filter(|count| *count > 0)
            .collect::<Vec<_>>();
        let moved = |measure: Measure, from: Status, to: Status| {
            compared
                .iter()
                .filter(|row| row.status(Measure::Primary) == from && row.status(measure) == to)
                .count()
        };
        let _ = write!(
            text,
            "\n## Coverage gaps (strict)\n\n{} samples have {} truth alleles in an HV region \
             inside the truth's analysed intervals where DNA has no callable base. {} samples \
             concordant on the primary measure are discordant on the strict one.\n\
             \n## IUPAC alleles (IUPAC-compatible)\n\n{} samples discordant on the primary measure \
             are concordant when DNA's base may match a truth IUPAC code that admits it; {} stay \
             discordant.\n",
            uncovered.len(),
            uncovered.iter().sum::<usize>(),
            moved(Measure::Strict, Status::Concordant, Status::Discordant),
            moved(
                Measure::IupacCompatible,
                Status::Discordant,
                Status::Concordant
            ),
            moved(
                Measure::IupacCompatible,
                Status::Discordant,
                Status::Discordant
            ),
        );
        self.haplotype_section(text);
    }

    fn haplotype_section(&self, text: &mut String) {
        let pairs = [
            (Measure::Primary, Measure::Haplotype),
            (Measure::IupacCompatible, Measure::HaplotypeIupacCompatible),
        ];
        text.push_str(
            "\n## Representation (haplotype)\n\nDiscordant samples on a base measure that the \
             haplotype measure finds concordant differ only in how the same sequence is \
             written, by the window of their base-measure discordant positions; the rest \
             keep blocks where the sequences differ (open: reaching the end of a comparable \
             interval with unequal lengths).\n\n\
             | Base measure | Haplotype measure | Discordant | Same sequence | Still discordant | \
             Blocks | Open blocks |",
        );
        for region in Region::ALL {
            let _ = write!(text, " Same: {} |", region.label());
        }
        text.push_str("\n|---|---|---:|---:|---:|---:|---:|");
        text.push_str(&"---:|".repeat(Region::ALL.len()));
        text.push('\n');
        for (base, measure) in pairs {
            let discordant = self
                .rows
                .iter()
                .filter(|row| row.status(base) == Status::Discordant)
                .collect::<Vec<_>>();
            let same = discordant
                .iter()
                .filter(|row| row.status(measure) == Status::Concordant)
                .filter_map(|row| row.comparison(base).ok())
                .collect::<Vec<_>>();
            let blocks = discordant
                .iter()
                .filter_map(|row| row.comparison(measure).ok())
                .flat_map(|comparison| &comparison.blocks)
                .collect::<Vec<_>>();
            let _ = write!(
                text,
                "| {} | {} | {} | {} | {} | {} | {} |",
                base.label(),
                measure.label(),
                discordant.len(),
                same.len(),
                discordant
                    .iter()
                    .filter(|row| row.status(measure) == Status::Discordant)
                    .count(),
                blocks.len(),
                blocks.iter().filter(|block| block.open).count(),
            );
            for region in Region::ALL {
                let count = same
                    .iter()
                    .filter(|comparison| regions(comparison).contains(&region))
                    .count();
                let _ = write!(text, " {count} |");
            }
            text.push('\n');
        }
    }

    fn anchor_section(&self, text: &mut String) {
        let anchors = self
            .rows
            .iter()
            .filter(|row| row.sample.anchor)
            .collect::<Vec<_>>();
        let concordant = |measure: Measure| {
            anchors
                .iter()
                .filter(|row| row.status(measure) == Status::Concordant)
                .count()
        };
        let _ = writeln!(
            text,
            "\n## Anchors\n\nOf {} anchor candidates, concordant: {}.",
            anchors.len(),
            Measure::ALL
                .map(|measure| format!("{} {}", measure.label(), concordant(measure)))
                .join(", ")
        );
        if let Some(baseline) = self.baseline {
            let set = anchors
                .iter()
                .filter(|row| {
                    baseline
                        .get(&row.sample.id)
                        .and_then(|statuses| statuses.get(&Measure::Primary))
                        .is_some_and(|status| Class::of(status) == Class::Good)
                })
                .collect::<Vec<_>>();
            let regressions = set
                .iter()
                .filter(|row| row.status(Measure::Primary) != Status::Concordant)
                .count();
            let _ = writeln!(
                text,
                "{} candidates were exact or concordant in the baseline (the anchor set); \
                 {regressions} of them are not concordant now.",
                set.len()
            );
        }
    }

    fn baseline_section(&self, text: &mut String) {
        let Some(baseline) = self.baseline else {
            return;
        };
        let mut transitions = BTreeMap::<(&str, &str), usize>::new();
        let mut changes = BTreeMap::<(Measure, Change), usize>::new();
        for row in self.rows {
            let before = baseline
                .get(&row.sample.id)
                .and_then(|statuses| statuses.get(&Measure::Primary))
                .map_or("(absent)", String::as_str);
            *transitions
                .entry((before, row.status(Measure::Primary).label()))
                .or_default() += 1;
            for measure in Measure::ALL {
                if let Some(change) = self.change(row, measure) {
                    *changes.entry((measure, change)).or_default() += 1;
                }
            }
        }
        for (_, statuses) in self.removed() {
            for measure in statuses.keys() {
                *changes.entry((*measure, Change::Removed)).or_default() += 1;
            }
        }
        text.push_str(
            "\n## Against the baseline\n\nPrimary status against the baseline's `status`:\n\n\
             | Baseline status | Status | Samples |\n|---|---|---:|\n",
        );
        for ((before, now), count) in transitions {
            let _ = writeln!(text, "| {before} | {now} | {count} |");
        }
        text.push_str("\nChanges of status class, per measure the baseline reports:\n\n| Change |");
        table_header(text, "|---|");
        for change in Change::ALL {
            let _ = write!(text, "| {} |", change.label());
            for measure in Measure::ALL {
                let count = changes.get(&(measure, change)).copied().unwrap_or_default();
                let _ = write!(text, " {count} |");
            }
            text.push('\n');
        }
    }
}

/// Completes a Markdown table header with one numeric column per measure:
/// the titles after the leading columns already written, then `separator`
/// (the leading columns' separator) and the measures' separators.
fn table_header(text: &mut String, separator: &str) {
    for measure in Measure::ALL {
        let _ = write!(text, " {} |", measure.title());
    }
    text.push('\n');
    text.push_str(separator);
    text.push_str(&"---:|".repeat(Measure::COUNT));
    text.push('\n');
}

/// Reason labels, comma-separated.
fn labels(reasons: &BTreeSet<Reason>) -> String {
    reasons
        .iter()
        .map(|reason| reason.label())
        .collect::<Vec<_>>()
        .join(",")
}

/// The windows of a comparison's discordant positions.
fn regions(comparison: &Comparison) -> BTreeSet<Region> {
    comparison
        .positions
        .iter()
        .map(|position| Region::of(position))
        .collect()
}

/// One TSV field: control characters become spaces.
fn clean(text: &str) -> String {
    text.chars()
        .map(|symbol| if symbol.is_control() { ' ' } else { symbol })
        .collect()
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use super::*;
    use crate::compare::Allele;
    use crate::runner::{Called, Identity, Threads};

    fn identity() -> Identity {
        Identity {
            configuration_sha256: "c".into(),
            reference_name: "r".into(),
            reference_sha256: "s".into(),
            profile_id: "p".into(),
            profile_sha256: "q".into(),
        }
    }

    /// From manifest and truth database to every output file, with a
    /// synthetic DNA: the files are byte-identical for one and for eight
    /// threads, and each status, change, and check appears.
    #[test]
    fn outputs_do_not_depend_on_the_thread_count() -> Result<(), Box<dyn std::error::Error>> {
        let directory = tempfile::tempdir()?;
        let first = directory.path().join("first.ab1");
        let second = directory.path().join("second.ab1");
        std::fs::write(&first, b"first trace")?;
        std::fs::write(&second, b"second trace")?;
        let traces = format!("{0};{1};{0}", first.display(), second.display());
        let lines = (0..30)
            .map(|index| {
                let (stratum, weight) = if index % 3 == 0 {
                    ("hard", "4.5")
                } else {
                    ("clean", "10")
                };
                format!(
                    "S{index:02}\tb1\t{stratum}\t{weight}\t{weight}\t{}\t{traces}",
                    u8::from(index % 5 == 0)
                )
            })
            .collect::<Vec<_>>();
        let text =
            crate::manifest::tests::manifest(&lines.iter().map(String::as_str).collect::<Vec<_>>());
        let manifest = crate::manifest::parse(text.as_bytes())?;

        // Every sample but S09 has a truth record; S08 sits in batch b2.
        let mut sql = String::from(crate::truth::tests::FIXTURE);
        sql.push_str(
            "DELETE FROM samples; DELETE FROM variants; DELETE FROM sample_status; \
             DELETE FROM intervals;",
        );
        for index in (0..30).filter(|index| *index != 9) {
            let batch = if index == 8 { "b2" } else { "b1" };
            let _ = write!(
                sql,
                "INSERT INTO samples VALUES ('S{index:02}', '{batch}', 2, 1, 0, 0, 0);
                 INSERT INTO sample_status VALUES ('{batch}', 'S{index:02}', 'x', 'OK', 1, 1, '', 1, 1);
                 INSERT INTO intervals VALUES ('S{index:02}', 'HV1', 16024, 16365),
                    ('S{index:02}', 'HV2', 73, 340);
                 INSERT INTO variants VALUES ('S{index:02}', 'snp', '263', 263, 'A', 'G'),
                    ('S{index:02}', 'insertion', '309.1', 309, '-', 'C'),
                    ('S{index:02}', 'snp', '16189', 16189, 'T', 'Y');"
            );
        }
        let path = crate::truth::tests::database(directory.path(), &sql)?;
        let wanted = manifest
            .samples
            .iter()
            .map(|sample| sample.id.as_str())
            .collect();
        let truth = crate::truth::load(&path, &wanted)?;

        let called = |sample: &Sample, traces: &[PathBuf]| {
            let index = sample.id[1..]
                .parse::<usize>()
                .map_err(|_| Failure::Panic)?;
            let alleles = match index % 5 {
                // Concordant except for the IUPAC site, which is outside DNA's intervals.
                0 => vec![Allele::new("263", "A", "G"), Allele::new("309.1", "-", "C")],
                // The major allele at the IUPAC site: discordant, IUPAC only.
                1 => vec![
                    Allele::new("263", "A", "G"),
                    Allele::new("309.1", "-", "C"),
                    Allele::new("16189", "T", "C"),
                ],
                2 => return Err(Failure::Dna(format!("{} traces", traces.len()))),
                // Nothing called: not applicable.
                3 => Vec::new(),
                // A DNA-only allele in HV2.
                _ => vec![
                    Allele::new("263", "A", "G"),
                    Allele::new("309.1", "-", "C"),
                    Allele::new("150", "C", "T"),
                ],
            };
            Ok(Called {
                rejected: index % 2,
                alleles,
                intervals: if index % 5 == 0 {
                    vec![(73, 340)]
                } else {
                    vec![(73, 340), (16024, 16365)]
                },
                identity: identity(),
            })
        };
        let call = |sample: &Sample, traces: &[PathBuf]| -> crate::runner::Ran {
            (called(sample, traces), None)
        };
        let baseline = crate::baseline::parse(
            b"sample_id\tstatus\tstrict_status\nS00\texact\tconcordant\nS01\texact\tdiscordant\n\
              S02\terror\terror\nGONE\texact\texact\n",
        )?;
        let check = TruthCheck::new(&manifest, &truth);
        assert_eq!(check.state(), TruthState::NotRecorded);
        let outputs = |threads: usize| -> Result<Vec<(&str, String)>, Box<dyn std::error::Error>> {
            let outcomes =
                crate::runner::run_with(&manifest.samples, Threads::new(threads)?, call)?;
            let rows = rows(
                &manifest.samples,
                &outcomes,
                &truth.samples,
                compare::tests::rcrs(),
            );
            Ok(Report::new(&rows, Some(&baseline), &check)
                .files(&[("panel", "synthetic".to_owned())]))
        };
        let one = outputs(1)?;
        assert_eq!(one, outputs(8)?);
        let names = one.iter().map(|(name, _)| *name).collect::<Vec<_>>();
        assert_eq!(
            names,
            [
                "scores.tsv",
                "estimates.tsv",
                "anchors.tsv",
                "diff.tsv",
                "summary.md"
            ]
        );
        let scores = &one[0].1;
        let expected = [
            // The IUPAC site lies outside DNA's callable intervals but inside
            // the truth's: concordant, but missed on the strict measure.
            "S00\tb1\tb1\thard\t4.5\t4.5\t1\t0\tconcordant\tY\t2\t3\t\t16189Y\t\t\t0\t\t\
             discordant\tN\t16189Y\tconcordant\tY\tconcordant\tY\t\tconcordant\tY\t\t2\t0\t",
            // The major allele at the IUPAC site: compatible with Y.
            "S01\tb1\tb1\tclean\t10\t10\t0\t0\tdiscordant\tN\t3\t3\t16189C\t16189Y\t\
             shared_position,dna_only,truth_only\t16189\t1\thv1_polyc\tdiscordant\tN\t\t\
             concordant\tY\tdiscordant\tN\t16189-16189:C/Y\tconcordant\tY\t\t2\t1\t",
            // Duplicate trace content is removed before DNA runs.
            "S02\tb1\tb1\tclean\t10\t10\t0\t0\tfailed\t\t\t\t\t\t\t\t\t\tfailed\t\t\tfailed\t\tfailed\t\t\tfailed\t\t\t\
             2\t\tdna: 2 traces",
            "S03\tb1\tb1\thard\t4.5\t4.5\t0\t0\tnot_applicable\tN/A\t0\t3\t\t263G 16189Y\t\t\t0\t\t\
             not_applicable\tN/A\t\tnot_applicable\tN/A\tnot_applicable\tN/A\t\t\
             not_applicable\tN/A\t\t2\t1\t",
            "S04\tb1\tb1\tclean\t10\t10\t0\t0\tdiscordant\tN\t3\t3\t150T\t16189Y\t\
             dna_only,truth_only\t150 16189\t0\thv1_polyc,other\tdiscordant\tN\t\tdiscordant\tN\t\
             discordant\tN\t150-150:T/C 16189-16189:T/Y\tdiscordant\tN\t150-150:T/C\t2\t0\t",
            "S08\tb1\tb2\tclean\t10\t10\t0\t0\tnot_applicable\tN/A\t0\t3\t\t263G 16189Y\t\t\t0\t\t\
             not_applicable\tN/A\t\tnot_applicable\tN/A\tnot_applicable\tN/A\t\t\
             not_applicable\tN/A\t\t2\t0\t",
            "S09\tb1\t\thard\t4.5\t4.5\t0\t0\tfailed\t\t\t\t\t\t\t\t\t\tfailed\t\t\tfailed\t\tfailed\t\t\tfailed\t\t\t\
             2\t1\ttruth_missing: the truth has no record of the sample",
        ];
        for line in expected {
            assert!(scores.contains(&format!("{line}\n")), "{line}");
        }

        let diff = &one[3].1;
        for line in [
            "S00\thard\tstrict\tconcordant\tdiscordant\tregression",
            "S01\tclean\tprimary\texact\tdiscordant\tregression",
            "S02\tclean\tprimary\terror\tfailed\tchanged",
            "S02\tclean\tstrict\terror\tfailed\tchanged",
            "S03\thard\tprimary\t\tnot_applicable\tadded",
            "GONE\t\tprimary\texact\t\tremoved",
            "GONE\t\tstrict\texact\t\tremoved",
        ] {
            assert!(diff.contains(&format!("{line}\n")), "{line}");
        }
        assert!(!diff.contains("S00\thard\tprimary"));
        assert!(!diff.contains("S01\tclean\tstrict"));
        // Added and removed samples count under every measure the baseline
        // reports, and under no other.
        assert!(diff.contains("S03\thard\tstrict\t\tnot_applicable\tadded\n"));
        assert!(!diff.contains("\tiupac_compatible\t"));
        let summary_counts = |summary: &str, change: &str| {
            summary
                .lines()
                .find(|line| line.starts_with(&format!("| {change} |")))
                .map(str::to_owned)
        };
        for (change, line) in [
            ("added", "| added | 27 | 27 | 0 | 0 | 0 |"),
            ("removed", "| removed | 1 | 1 | 0 | 0 | 0 |"),
            ("regression", "| regression | 1 | 1 | 0 | 0 | 0 |"),
        ] {
            assert_eq!(summary_counts(&one[4].1, change).as_deref(), Some(line));
            for (measure, column) in [("primary", 2), ("strict", 3), ("iupac_compatible", 4)] {
                let rows = diff
                    .lines()
                    .filter(|row| row.contains(&format!("\t{measure}\t")) && row.ends_with(change))
                    .count();
                let summary = line.split(" | ").nth(column - 1).unwrap_or_default();
                assert_eq!(
                    rows.to_string(),
                    summary.trim_end_matches(" |"),
                    "{measure} {change}"
                );
            }
        }
        let anchors = &one[2].1;
        assert!(anchors.contains(
            "S00\thard\tconcordant\tdiscordant\tconcordant\tconcordant\tconcordant\texact\tsame\n"
        ));
        let estimates = &one[1].1;
        assert!(estimates.contains("corpus\tstrict\tnot_concordant\tall\t"));
        assert!(estimates.contains("pool\tall\tfailed\tall\t"));
        let summary = &one[4].1;
        for line in [
            "| without a truth record | 1 |",
            "| truth batch differs from the panel batch | 1 |",
            "| exact | discordant | 1 |",
            "| regression | 1 | 1 | 0 | 0 | 0 |",
        ] {
            assert!(summary.contains(line), "{line}");
        }
        Ok(())
    }

    #[test]
    fn a_retired_sample_weighs_nothing() -> Result<(), crate::error::ManifestError> {
        let manifest = crate::manifest::parse(
            crate::manifest::tests::manifest(&["S1\tb1\tclean\t3\t6\t0\t/t/a.ab1"]).as_bytes(),
        )?;
        let mut sample = manifest.samples[0].clone();
        assert_eq!(Weighting::Corpus.of(&sample), 6.0);
        sample.retired = true;
        assert_eq!(Weighting::Corpus.of(&sample), 0.0);
        assert_eq!(Weighting::Pool.of(&sample), 0.0);
        Ok(())
    }

    #[test]
    fn a_changed_truth_is_reported_and_listed() {
        let digest = |text: &str| text.repeat(64);
        let check = TruthCheck {
            recorded: Some(digest("a")),
            current: digest("b"),
            changed: vec![crate::truth::ChangedSample {
                sample: "S1".into(),
                stratum: "clean".into(),
                retired: false,
                recorded: digest("c"),
                current: digest("d"),
            }],
        };
        let files = Report::new(&[], None, &check).files(&[]);
        let names = files.iter().map(|(name, _)| *name).collect::<Vec<_>>();
        assert_eq!(
            names,
            [
                "scores.tsv",
                "estimates.tsv",
                "anchors.tsv",
                "truth_changes.tsv",
                "summary.md"
            ]
        );
        assert!(
            files[3]
                .1
                .contains(&format!("S1\tclean\t0\t{}\t{}\n", digest("c"), digest("d")))
        );
        let summary = &files[4].1;
        assert!(summary.starts_with("# Evaluation run\n\n**Truth changed:** the truth of 1 panel"));
        assert!(summary.contains("| truth changed since the panel was built | 1 |"));
        let unchanged = TruthCheck {
            recorded: Some(digest("b")),
            changed: Vec::new(),
            ..check
        };
        let files = Report::new(&[], None, &unchanged).files(&[]);
        assert!(files.iter().all(|(name, _)| *name != "truth_changes.tsv"));
        assert!(files[3].1.contains("Panel truth: unchanged"));
        // The header digest differs while no row does: the file is still
        // written, and the notice says why it lists nothing.
        let header_only = TruthCheck {
            recorded: Some(digest("a")),
            changed: Vec::new(),
            ..unchanged
        };
        let files = Report::new(&[], None, &header_only).files(&[]);
        let tsv = files
            .iter()
            .find(|(name, _)| *name == "truth_changes.tsv")
            .map(|(_, text)| text.as_str());
        assert_eq!(
            tsv,
            Some("sample_id\tstratum\tretired\trecorded_digest\tcurrent_digest\n")
        );
        let summary = &files[files.len() - 1].1;
        assert!(summary.contains("although no sample's recorded digest does"));
    }
}
