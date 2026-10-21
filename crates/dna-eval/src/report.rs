//! The run's outputs: per-sample comparisons, estimates, anchors, and a
//! Markdown summary.
//!
//! Every sample is compared under each [`Measure`]. Every output but the run
//! header (`run.json`, which records timings) is a function of the panel, the
//! truth, and DNA's calls only, so two runs of the same inputs give
//! byte-identical files whatever the thread count.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

use crate::cause::forms::Forms;
use crate::cause::{self, Cause, Finding, Scope};
use crate::classes::{Discrepancy, StatusClass};
use crate::compare::{
    self, Comparison, Concordance, Failure, Measure, Reason, Region, Side, Yield,
};
use crate::estimate::{self, Estimate, Observation};
use crate::haplotype;
use crate::lod;
use crate::manifest::{Month, Sample};
use crate::recheck::Recheck;
use crate::rundir;
use crate::runner::{Empop, Outcome};
use crate::truth::{TruthCheck, TruthSample, TruthState};
use crate::unrecorded;
use crate::workload::{self, Class as BlockClass, Escaped, Passes, Tract, Unit, Window, Workload};
use crate::zones::{Kind, Zones};
use dna_kernel::profile::Convention;

/// One compared panel sample.
pub(crate) struct Row<'a> {
    pub(crate) sample: &'a Sample,
    pub(crate) outcome: &'a Outcome,
    /// The sample's truth record, when the truth has one.
    pub(crate) truth: Option<&'a TruthSample>,
    /// One comparison per measure, in [`Measure::ALL`] order.
    pub(crate) comparisons: Result<[Comparison; Measure::COUNT], Failure>,
    /// Truth IUPAC positions below the limit of detection, which
    /// `haplotype_lod` compares IUPAC-compatibly.
    pub(crate) below_lod: usize,
    /// Per HV region, the truth's analysed positions DNA's callable
    /// intervals cover.
    pub(crate) callable_yield: Option<[Yield; 3]>,
}

/// Compares every sample's outcome, given in sample order, with its truth
/// under every measure, on `reference` (the reference sequence DNA ran on),
/// and sorts the rows by sample identifier.
pub(crate) fn rows<'a>(
    samples: &'a [Sample],
    outcomes: &'a [Outcome],
    truth: &'a BTreeMap<String, TruthSample>,
    reference: &str,
    conventions: &[Convention],
) -> Vec<Row<'a>> {
    let mut rows = samples
        .iter()
        .zip(outcomes)
        .map(|(sample, outcome)| {
            let record = truth.get(&sample.id);
            let mut below_lod = 0;
            let mut callable_yield = None;
            let comparisons = match (&outcome.result, record) {
                (Err(failure), _) => Err(failure.clone()),
                (Ok(_), None) => Err(Failure::TruthMissing),
                (Ok(called), Some(record)) => {
                    // Run lengths the reviewers do not record are not
                    // compared, on either side.
                    let dna_alleles = unrecorded::recorded(&called.alleles, reference, conventions);
                    let truth_alleles =
                        unrecorded::recorded(&record.alleles, reference, conventions);
                    let dna = Side {
                        alleles: &dna_alleles,
                        intervals: &called.intervals,
                    };
                    let truth = Side {
                        alleles: &truth_alleles,
                        intervals: &record.intervals,
                    };
                    callable_yield = Some(compare::callable_yield(
                        &called.intervals,
                        &record.intervals,
                    ));
                    let below = lod::below(&record.alleles, &called.signals);
                    below_lod = below.len();
                    let haplotypes = compare::haplotypes(dna, truth, reference, None);
                    let lod_haplotypes = compare::haplotypes(dna, truth, reference, Some(&below));
                    Ok(Measure::ALL.map(|measure| {
                        let haplotypes = if measure == Measure::HaplotypeLod {
                            &lod_haplotypes
                        } else {
                            &haplotypes
                        };
                        compare::compare(dna, truth, measure, haplotypes, &below)
                    }))
                }
            };
            Row {
                sample,
                outcome,
                truth: record,
                comparisons,
                below_lod,
                callable_yield,
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

    /// The sample's discrepancy class.
    pub(crate) fn discrepancy(&self) -> Discrepancy {
        Discrepancy::of(
            |measure| StatusClass::of(self.status(measure).label()),
            self.unanalysed(),
        )
    }

    /// Whether every primary difference is a DNA call outside the truth's
    /// analysed intervals.
    fn unanalysed(&self) -> bool {
        let (Ok(primary), Some(truth)) = (self.comparison(Measure::Primary), self.truth) else {
            return false;
        };
        Discrepancy::only_unanalysed(&primary.reasons, &primary.positions, &truth.intervals)
    }

    /// The sample's cause, scope, and findings, on the truth's written `forms` and
    /// `reference`; `None` for a sample without a discrepancy.
    pub(crate) fn cause(
        &self,
        forms: &Forms,
        reference: &str,
    ) -> Option<(Cause, Scope, Vec<Finding>)> {
        let discrepancy = self.discrepancy();
        let regions = self.discrepancy_regions(discrepancy);
        let called = self.outcome.result.as_ref().ok();
        cause::findings(
            &cause::Sample {
                discrepancy,
                failure: self.comparisons.as_ref().err(),
                comparison: discrepancy
                    .measure()
                    .and_then(|measure| self.comparison(measure).ok()),
                primary: self.comparison(Measure::Primary).ok(),
                regions: &regions,
                dna: called.map_or(&[][..], |called| called.alleles.as_slice()),
                not_recorded: called.map(|called| &called.not_recorded),
                truth: self.truth.map_or(&[][..], |truth| truth.alleles.as_slice()),
                loci: called.map(|called| &called.loci),
                signals: called.map(|called| &called.signals),
            },
            forms,
            reference,
        )
    }

    /// The windows of the class's discordant positions.
    fn discrepancy_regions(&self, discrepancy: Discrepancy) -> BTreeSet<Region> {
        discrepancy
            .measure()
            .and_then(|measure| self.comparison(measure).ok())
            .map(regions)
            .unwrap_or_default()
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

/// Every derived output of one run.
/// The header of `unreviewed.tsv`.
pub(crate) const UNREVIEWED_HEADER: &str =
    "sample_id\tfirst\tlast\tblock\tclass\tzone\tzone_checked\n";

pub(crate) struct Report<'a> {
    rows: &'a [Row<'a>],
    check: &'a TruthCheck,
    estimates: Vec<(Weighting, Metric, Estimate)>,
    /// Per row, its cause and findings.
    causes: Vec<Option<(Cause, Scope, Vec<Finding>)>>,
    /// The review workload of every row.
    review: Review,
    /// The first batch month (`YYYY-MM`) of the current data, whose slice
    /// the summary reports first.
    current_since: Option<Month>,
}

/// Default window of the review workload, in positions.
pub(crate) const REVIEW_WINDOW: i64 = 10;

/// The review workload of a run: its windows, its forced windows, its
/// re-check list, and each row's workload when the row was compared.
struct Review {
    window: i64,
    windows: Vec<Window>,
    forced: Option<Vec<(i64, i64)>>,
    zones: Option<Zones>,
    recheck: Option<Recheck>,
    /// The profile's length tracts, whose length differences are noted.
    tracts: Vec<Tract>,
    workloads: Vec<Option<Workload>>,
}

impl Review {
    fn of(
        rows: &[Row<'_>],
        window: i64,
        forced: Option<Vec<(i64, i64)>>,
        zones: Option<Zones>,
        recheck: Option<Recheck>,
        tracts: Vec<Tract>,
    ) -> Self {
        let windows = workload::windows(window);
        let workloads = rows
            .iter()
            .map(|row| {
                let called = row.outcome.result.as_ref().ok()?;
                let comparison = row.comparison(Measure::HaplotypeLod).ok()?;
                Some(workload::of(
                    &called.checks,
                    &comparison.blocks,
                    &windows,
                    forced.as_deref().unwrap_or_default(),
                    zones.as_ref(),
                    (
                        recheck
                            .as_ref()
                            .map_or(&[], |recheck| recheck.spans(&row.sample.id)),
                        &tracts,
                    ),
                ))
            })
            .collect();
        Self {
            window,
            windows,
            forced,
            zones,
            recheck,
            tracts,
            workloads,
        }
    }

    /// The rows the workload counts: compared under `haplotype_lod` with a
    /// concordant or discordant status, with their workloads.
    fn counted<'r, 'a>(&'r self, rows: &'r [Row<'a>]) -> Vec<(&'r Row<'a>, &'r Workload)> {
        rows.iter()
            .zip(&self.workloads)
            .filter(|(row, _)| {
                matches!(
                    row.status(Measure::HaplotypeLod),
                    Status::Concordant | Status::Discordant
                )
            })
            .filter_map(|(row, workload)| Some((row, workload.as_ref()?)))
            .collect()
    }
}

impl<'a> Report<'a> {
    /// Derives the estimates and causes of `rows`, which are sorted by sample
    /// identifier, the causes on the truth's written `forms` and `reference`.
    pub(crate) fn new(
        rows: &'a [Row<'a>],
        check: &'a TruthCheck,
        (forms, reference): (&Forms, &str),
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
        let causes = rows.iter().map(|row| row.cause(forms, reference)).collect();
        Self {
            rows,
            check,
            estimates,
            causes,
            review: Review::of(rows, REVIEW_WINDOW, None, None, None, Vec::new()),
            current_since: None,
        }
    }

    /// The compared rows, sorted by sample identifier.
    pub(crate) const fn rows(&self) -> &'a [Row<'a>] {
        self.rows
    }

    /// The check of the panel against the truth.
    pub(crate) const fn check(&self) -> &'a TruthCheck {
        self.check
    }

    /// The same report with the samples whose batch date is in or after
    /// `since` (`YYYY-MM`) reported first as the current data.
    pub(crate) fn with_current_since(mut self, since: Option<Month>) -> Self {
        self.current_since = since;
        self
    }

    /// Whether a row belongs to the current data.
    fn current(&self, row: &Row<'_>) -> bool {
        self.current_since
            .is_some_and(|since| row.sample.batch_month.is_some_and(|month| month >= since))
    }

    /// The same report with the review workload over windows of `window`
    /// positions, `forced` windows counted as reviewed when given, the check
    /// items placed in `zones` when given, and the blocks of the `recheck`
    /// list in their own class when given.
    pub(crate) fn with_review(
        mut self,
        window: i64,
        forced: Option<Vec<(i64, i64)>>,
        zones: Option<Zones>,
        recheck: Option<Recheck>,
        tracts: Vec<Tract>,
    ) -> Self {
        self.review = Review::of(self.rows, window, forced, zones, recheck, tracts);
        self
    }

    /// Every deterministic output file, by name; `truth_changes.tsv` only
    /// when the panel's truth changed.
    /// `context` heads the summary.
    pub(crate) fn files(&self, context: &[(&str, String)]) -> Vec<(&'static str, String)> {
        let mut files = vec![
            (rundir::SCORES, self.scores()),
            ("estimates.tsv", self.estimates()),
            ("anchors.tsv", self.anchors()),
            ("windows.tsv", self.windows()),
        ];
        if self.check.state() == TruthState::Changed {
            files.push(("truth_changes.tsv", self.check.changes()));
        }
        files.push((rundir::UNREVIEWED, self.unreviewed()));
        files.push(("summary.md", self.summary(context)));
        files
    }

    /// `unreviewed.tsv`: one row per sample and unreviewed difference block,
    /// in sample order then block order: its anchors, its text, its class as
    /// the pass measures split it, the zone that holds its first anchor, and
    /// whether a check item touches that zone, so that a later comparison of
    /// two runs (`dna-eval diff`) counts the same passes from the files.
    fn unreviewed(&self) -> String {
        let mut text = String::from(UNREVIEWED_HEADER);
        for (row, workload) in self.rows.iter().zip(&self.review.workloads) {
            let Some(workload) = workload else {
                continue;
            };
            for (escape, block) in workload.escapes.iter().zip(&workload.unreviewed) {
                let (zone, checked) = match (escape.zone, self.review.zones.as_ref()) {
                    (Some(instance), Some(zones)) => (
                        zones.label(instance),
                        if workload.zones.contains(&instance) {
                            "yes"
                        } else {
                            "no"
                        },
                    ),
                    _ => (".".to_owned(), "."),
                };
                let _ = writeln!(
                    text,
                    "{}\t{}\t{}\t{block}\t{}\t{zone}\t{checked}",
                    row.sample.id,
                    escape.first,
                    escape.last,
                    escape.class.label(),
                );
            }
        }
        text
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
             haplotype_iupac_compatible_blocks\thaplotype_lod_status\thaplotype_lod_concordant\t\
             haplotype_lod_blocks\tbelow_lod\tclass\tclass_regions\tcause\tcause_scope\t\
             cause_evidence\tempop_ranges\tempop_mitotype\tempop_undecided\tempop_ties\t\
             callable_yield\t\
             traces\trejected_traces\tfailure\t\
             review_checks\treview_positions\treview_clear\tunreviewed\tunreviewed_forced\t\
             review_zones\treview_zone_count\n",
        );
        for ((row, cause), workload) in self
            .rows
            .iter()
            .zip(&self.causes)
            .zip(&self.review.workloads)
        {
            let sample = row.sample;
            let review = workload.as_ref().map_or_else(
                || "\t\t\t\t\t\t".to_owned(),
                |workload| {
                    let (zones, count) = self.review.zones.as_ref().map_or_else(
                        || (String::new(), String::new()),
                        |zones| {
                            (
                                workload
                                    .zones
                                    .iter()
                                    .map(|&instance| zones.label(instance))
                                    .collect::<Vec<_>>()
                                    .join(" "),
                                workload.zones.len().to_string(),
                            )
                        },
                    );
                    format!(
                        "{}\t{}\t{}\t{}\t{}\t{zones}\t{count}",
                        workload.checks,
                        workload.positions,
                        if workload.checks == 0 { "Y" } else { "N" },
                        workload.unreviewed.join(" "),
                        if self.review.forced.is_some() {
                            workload.unreviewed_forced.join(" ")
                        } else {
                            String::new()
                        },
                    )
                },
            );
            let (cause, cause_scope, cause_evidence) = cause.as_ref().map_or_else(
                || (String::new(), String::new(), String::new()),
                |(cause, scope, findings)| {
                    (
                        cause.label().to_owned(),
                        scope.label().to_owned(),
                        findings
                            .iter()
                            .map(Finding::text)
                            .collect::<Vec<_>>()
                            .join(" "),
                    )
                },
            );
            let primary = row.status(Measure::Primary);
            let strict = row.status(Measure::Strict);
            let compatible = row.status(Measure::IupacCompatible);
            let haplotype = [
                Measure::Haplotype,
                Measure::HaplotypeIupacCompatible,
                Measure::HaplotypeLod,
            ]
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
            let empop = row
                .outcome
                .result
                .as_ref()
                .ok()
                .and_then(|called| called.empop.as_ref())
                .map_or_else(
                    || "\t\t\t".to_owned(),
                    |empop| {
                        format!(
                            "{}\t{}\t{}\t{}",
                            empop.ranges, empop.mitotype, empop.undecided, empop.ties
                        )
                    },
                );
            let callable_yield = row.callable_yield.map_or_else(String::new, |yields| {
                yields.map(|value| value.text()).join(" ")
            });
            let below_lod = if row.comparisons.is_ok() {
                row.below_lod.to_string()
            } else {
                String::new()
            };
            let discrepancy = row.discrepancy();
            let class_regions = row
                .discrepancy_regions(discrepancy)
                .iter()
                .map(|region| region.label())
                .collect::<Vec<_>>()
                .join(",");
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
                 {haplotype}\t{below_lod}\t{}\t{class_regions}\t{cause}\t{cause_scope}\t\
                 {cause_evidence}\t{empop}\t{callable_yield}\t\
                 {traces}\t{rejected}\t{failure}\t{review}",
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
                discrepancy.label(),
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

    /// Per window of the review workload: the counted samples whose
    /// analysed intervals reach it, those with a check item in it, and those
    /// with an unreviewed difference in it, as a share of the ones without a
    /// check item there.
    fn window_counts(&self) -> Vec<(usize, usize, usize, usize)> {
        let counted = self.review.counted(self.rows);
        self.review
            .windows
            .iter()
            .enumerate()
            .map(|(index, window)| {
                let mut samples = 0;
                let mut checked = 0;
                let mut open = 0;
                let mut escaped = 0;
                for (row, workload) in &counted {
                    let analysed = row.truth.is_some_and(|truth| {
                        truth
                            .intervals
                            .iter()
                            .any(|&(first, last)| window.overlaps(first, last))
                    });
                    if !analysed {
                        continue;
                    }
                    samples += 1;
                    if workload.checked[index] {
                        checked += 1;
                    } else {
                        open += 1;
                        escaped += usize::from(workload.escaped[index]);
                    }
                }
                (samples, checked, open, escaped)
            })
            .collect()
    }

    /// `windows.tsv`: one row per review window, in region order.
    fn windows(&self) -> String {
        let mut text = String::from(
            "region\tstart\tend\tsamples\tchecked\tunreviewed\tunreviewed_share\t\
             unreviewed_upper95\n",
        );
        for (window, (samples, checked, open, escaped)) in
            self.review.windows.iter().zip(self.window_counts())
        {
            let share = if open == 0 {
                String::new()
            } else {
                format!("{:.4}", escaped as f64 / open as f64)
            };
            let _ = writeln!(
                text,
                "{}\t{}\t{}\t{samples}\t{checked}\t{escaped}\t{share}\t{:.4}",
                window.region,
                window.start,
                window.end,
                workload::wilson_upper(escaped, open),
            );
        }
        text
    }

    /// The review workload: positions to review, samples with an
    /// unreviewed difference overall, by region, and by batch month, clear
    /// samples, the windows with the most unreviewed differences, and the
    /// forced windows' effect.
    fn review_section(&self, text: &mut String) {
        let counted = self.review.counted(self.rows);
        if counted.is_empty() {
            return;
        }
        if let Some(since) = &self.current_since {
            let current = counted
                .iter()
                .filter(|(row, _)| self.current(row))
                .collect::<Vec<_>>();
            if !current.is_empty() {
                let span = self.review.windows.iter().map(Window::length).sum::<i64>();
                let total = span * i64::try_from(current.len()).unwrap_or(i64::MAX);
                let positions = current
                    .iter()
                    .map(|(_, workload)| workload.positions)
                    .sum::<i64>();
                let unreviewed = current
                    .iter()
                    .filter(|(_, workload)| !workload.unreviewed.is_empty())
                    .count();
                let clear = current
                    .iter()
                    .filter(|(_, workload)| workload.checks == 0)
                    .count();
                let share = |part: f64, whole: f64| {
                    if whole > 0.0 {
                        100.0 * part / whole
                    } else {
                        0.0
                    }
                };
                let _ = write!(
                    text,
                    "\n## Review workload on current data\n\nSamples with a batch date in or after \
                     {since}: {}.\n\n| Measure | Value |\n|---|---:|\n\
                     | Positions to review | {positions} ({:.2} %) |\n\
                     | Samples with an unreviewed difference | {unreviewed} ({:.2} %) |\n\
                     | Clear samples (no check item) | {clear} ({:.2} %) |\n",
                    current.len(),
                    share(positions as f64, total as f64),
                    share(unreviewed as f64, current.len() as f64),
                    share(clear as f64, current.len() as f64),
                );
                let current = current.into_iter().copied().collect::<Vec<_>>();
                self.pass_section(text, &current);
                self.zone_section(text, &current);
                self.note_section(text, &current);
            }
        }
        let samples = counted.len();
        let span = self.review.windows.iter().map(Window::length).sum::<i64>();
        let total = span * i64::try_from(samples).unwrap_or(i64::MAX);
        let percent = |part: f64, whole: f64| {
            if whole > 0.0 {
                100.0 * part / whole
            } else {
                0.0
            }
        };
        let positions = counted
            .iter()
            .map(|(_, workload)| workload.positions)
            .sum::<i64>();
        let unreviewed = counted
            .iter()
            .filter(|(_, workload)| !workload.unreviewed.is_empty())
            .count();
        let clear = counted
            .iter()
            .filter(|(_, workload)| workload.checks == 0)
            .count();
        let _ = write!(
            text,
            "\n## Review workload\n\nOver {samples} samples compared under `haplotype_lod`, \
             in windows of {} positions over the HV regions. A reviewer looks at every window \
             a check item touches; a difference block no check item overlaps is unreviewed.\n\n\
             | Measure | Value |\n|---|---:|\n\
             | Positions to review | {positions} ({:.2} %) |\n\
             | Samples with an unreviewed difference | {unreviewed} ({:.2} %) |\n\
             | Clear samples (no check item) | {clear} ({:.2} %) |\n",
            self.review.window,
            percent(positions as f64, total as f64),
            percent(unreviewed as f64, samples as f64),
            percent(clear as f64, samples as f64),
        );
        if self.review.forced.is_some() {
            let positions = counted
                .iter()
                .map(|(_, workload)| workload.positions_forced)
                .sum::<i64>();
            let unreviewed = counted
                .iter()
                .filter(|(_, workload)| !workload.unreviewed_forced.is_empty())
                .count();
            let _ = write!(
                text,
                "| With forced windows: positions to review | {positions} ({:.2} %) |\n\
                 | With forced windows: samples with an unreviewed difference | {unreviewed} \
                 ({:.2} %) |\n",
                percent(positions as f64, total as f64),
                percent(unreviewed as f64, samples as f64),
            );
        }
        if let Some(recheck) = &self.review.recheck {
            let matched = counted
                .iter()
                .filter_map(|(row, _)| {
                    let blocks = &row.comparison(Measure::HaplotypeLod).ok()?.blocks;
                    Some(
                        recheck
                            .spans(&row.sample.id)
                            .iter()
                            .filter(|&&(first, last)| {
                                blocks.iter().any(|block| block.overlaps(first, last))
                            })
                            .count(),
                    )
                })
                .sum::<usize>();
            let _ = writeln!(
                text,
                "| Re-check list rows overlapping no difference block of these samples | {} \
                 of {} |",
                recheck.rows() - matched,
                recheck.rows()
            );
        }
        self.pass_section(text, &counted);
        self.zone_section(text, &counted);
        self.note_section(text, &counted);
        text.push_str("\nSamples with an unreviewed difference, by region:\n\n| Region | Samples |\n|---|---:|\n");
        let mut regions = BTreeMap::<&str, usize>::new();
        for (_, workload) in &counted {
            for region in workload.regions.iter().collect::<BTreeSet<_>>() {
                *regions.entry(region).or_default() += 1;
            }
        }
        for (region, count) in regions {
            let _ = writeln!(text, "| {region} | {count} |");
        }
        text.push_str(
            "\nBy batch month:\n\n| Month | Samples | With an unreviewed difference |\n|---|---:|---:|\n",
        );
        let mut months = BTreeMap::<String, (usize, usize)>::new();
        for (row, workload) in &counted {
            let month = row
                .sample
                .batch_month
                .map_or_else(|| "unknown".to_owned(), |month| month.to_string());
            let entry = months.entry(month).or_default();
            entry.0 += 1;
            entry.1 += usize::from(!workload.unreviewed.is_empty());
        }
        for (month, (samples, unreviewed)) in months {
            let _ = writeln!(text, "| {month} | {samples} | {unreviewed} |");
        }
        let mut ranked = self
            .review
            .windows
            .iter()
            .zip(self.window_counts())
            .filter(|(_, (_, _, _, escaped))| *escaped > 0)
            .collect::<Vec<_>>();
        ranked.sort_by_key(|(window, (_, _, _, escaped))| {
            (std::cmp::Reverse(*escaped), window.start)
        });
        text.push_str(
            "\nWindows with the most unreviewed differences (forcing a window is a policy \
             decision, learned on one run and judged on another):\n\n\
             | Window | Samples | Checked | Unreviewed | Share of unchecked |\n|---|---:|---:|---:|---:|\n",
        );
        for (window, (samples, checked, open, escaped)) in ranked.into_iter().take(20) {
            let _ = writeln!(
                text,
                "| {} {}-{} | {samples} | {checked} | {escaped} | {:.2} % |",
                window.region,
                window.start,
                window.end,
                percent(escaped as f64, open as f64),
            );
        }
    }

    /// What passes without review over every counted sample and, with a
    /// current slice, over the current data: the measures `run.json`
    /// records.
    pub(crate) fn passes(&self) -> (Passes, Option<Passes>) {
        let counted = self.review.counted(self.rows);
        let zones = self.review.zones.as_ref();
        let all = workload::passes(counted.iter().map(|(_, workload)| *workload), zones);
        let current = self.current_since.as_ref().map(|_| {
            workload::passes(
                counted
                    .iter()
                    .filter(|(row, _)| self.current(row))
                    .map(|(_, workload)| *workload),
                zones,
            )
        });
        (all, current)
    }

    /// What passes without review of `samples`: the samples without a check
    /// item, the regions and zones no check item touches, overall, per
    /// region, per zone, and per zone kind, and the passed samples, regions,
    /// and zones that carry an unreviewed difference, by the class of their
    /// blocks; with a re-check list, the unreviewed blocks on it.
    fn pass_section(&self, text: &mut String, samples: &[(&Row<'_>, &Workload)]) {
        let Some(geometry) = &self.review.zones else {
            return;
        };
        if samples.is_empty() {
            return;
        }
        let passes = workload::passes(
            samples.iter().map(|(_, workload)| *workload),
            Some(geometry),
        );
        let percent = |part: usize, whole: usize| {
            if whole > 0 {
                100.0 * part as f64 / whole as f64
            } else {
                0.0
            }
        };
        let split = |escaped: &Escaped| {
            format!(
                "non-het {}, truth N only {}, het only {}, re-check only {}, length note only {}",
                escaped.non_het,
                escaped.truth_n_only,
                escaped.het_only,
                escaped.recheck_only,
                escaped.length_note_only
            )
        };
        let (auto_pass, regions, zones) = (&passes.auto_pass, &passes.regions, &passes.zones);
        let unchecked = regions.total - regions.checked;
        let open = zones.total - zones.checked;
        let _ = write!(
            text,
            "\nPasses (a sample, region, or zone that no check item touches passes unreviewed; \
             a difference block counts for the region of its first position. A block is \
             re-check when the re-check list names it, else truth N when the truth writes `N` \
             in it, else length note when it differs only in the length of a profile length \
             tract, else het when the truth writes an IUPAC code of two or three bases in it, \
             else non-het; a passed unit counts under its most severe block, in the order \
             non-het, truth N, het, re-check, length note; length notes are never gated):\n\n\
             | Measure | Value |\n|---|---:|\n\
             | Auto-pass samples (no check item) | {} of {} ({:.2} %) |\n\
             | Auto-pass samples with an unreviewed difference | {} ({:.2} % of them): {} |\n\
             | Regions needing a check | {} of {} ({:.2} %) |\n\
             | Passed regions | {unchecked} ({:.2} %) |\n\
             | Passed regions with an unreviewed difference | {} ({:.2} % of them): {} |\n\
             | Zones needing a check | {} of {} ({:.2} %) |\n\
             | Passed zones | {open} ({:.2} %) |\n\
             | Passed zones with an unreviewed difference | {} ({:.2} % of them): {} |\n",
            auto_pass.samples,
            passes.samples,
            percent(auto_pass.samples, passes.samples),
            auto_pass.unreviewed.total,
            percent(auto_pass.unreviewed.total, auto_pass.samples),
            split(&auto_pass.unreviewed),
            regions.checked,
            regions.total,
            percent(regions.checked, regions.total),
            percent(unchecked, regions.total),
            regions.unreviewed.total,
            percent(regions.unreviewed.total, unchecked),
            split(&regions.unreviewed),
            zones.checked,
            zones.total,
            percent(zones.checked, zones.total),
            percent(open, zones.total),
            zones.unreviewed.total,
            percent(zones.unreviewed.total, open),
            split(&zones.unreviewed),
        );
        if self.review.recheck.is_some() {
            let listed = samples
                .iter()
                .flat_map(|(_, workload)| &workload.escapes)
                .filter(|escape| escape.class == BlockClass::Recheck)
                .count();
            let _ = writeln!(
                text,
                "| Unreviewed blocks on the re-check list | {listed} |"
            );
        }
        text.push_str(
            "\n| Region | Needing a check | Passed | Passed with an unreviewed difference \
             | Non-het | Truth N only | Het only | Re-check only | Length note only |\n\
             |---|---:|---:|---:|---:|---:|---:|---:|---:|\n",
        );
        for region in &regions.by_region {
            let _ = writeln!(
                text,
                "| {} | {} ({:.2} %) | {:.2} % | {} | {} | {} | {} | {} | {} |",
                region.name,
                region.checked,
                percent(region.checked, region.total),
                percent(region.total - region.checked, region.total),
                region.unreviewed.total,
                region.unreviewed.non_het,
                region.unreviewed.truth_n_only,
                region.unreviewed.het_only,
                region.unreviewed.recheck_only,
                region.unreviewed.length_note_only,
            );
        }
        let table = |text: &mut String, label: &str, units: &[Unit]| {
            let _ = write!(
                text,
                "\n| {label} | Zones | Needing a check | Passed | Passed with an unreviewed \
                 difference | Non-het | Truth N only | Het only | Re-check only | Length note \
                 only |\n|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n"
            );
            for unit in units {
                let _ = writeln!(
                    text,
                    "| {} | {} | {} ({:.2} %) | {:.2} % | {} | {} | {} | {} | {} | {} |",
                    unit.name,
                    unit.total,
                    unit.checked,
                    percent(unit.checked, unit.total),
                    percent(unit.total - unit.checked, unit.total),
                    unit.unreviewed.total,
                    unit.unreviewed.non_het,
                    unit.unreviewed.truth_n_only,
                    unit.unreviewed.het_only,
                    unit.unreviewed.recheck_only,
                    unit.unreviewed.length_note_only,
                );
            }
        };
        table(text, "Zone", &zones.by_zone);
        let instances = geometry.instances();
        let kinds = Kind::ALL
            .into_iter()
            .map(|kind| {
                let units = instances
                    .iter()
                    .zip(&zones.by_zone)
                    .filter(|((_, found), _)| *found == kind)
                    .map(|(_, unit)| unit)
                    .collect::<Vec<_>>();
                Unit {
                    name: kind.label().to_owned(),
                    total: units.iter().map(|unit| unit.total).sum(),
                    checked: units.iter().map(|unit| unit.checked).sum(),
                    unreviewed: Escaped::sum(units.iter().map(|unit| &unit.unreviewed)),
                }
            })
            .filter(|unit| unit.total > 0)
            .collect::<Vec<_>>();
        table(text, "Zone kind", &kinds);
    }

    /// The noted length tracts of `samples`: per profile length tract, the
    /// items the review notes (`tract_length`) and their samples, and the
    /// unreviewed blocks that differ only in the tract's length
    /// (`length_note`, never gated), split by the truth's length: the
    /// note's convention reading, DNA's length less one, or another.
    fn note_section(&self, text: &mut String, samples: &[(&Row<'_>, &Workload)]) {
        if self.review.tracts.is_empty() || samples.is_empty() {
            return;
        }
        text.push_str(
            "\nLength notes (a run of a profile length tract whose dominant length the reads \
             determine, or that DNA does not report, is noted, not checked; an unreviewed \
             block that differs only in a tract's length is a length note, reported and never \
             gated, split by the truth's length: the note's convention reading, DNA's length \
             less one, or another):\n\n\
             | Tract | Noted items | Samples | Length-note blocks | Convention | Mode - 1 | \
             Other |\n|---|---:|---:|---:|---:|---:|---:|\n",
        );
        for tract in &self.review.tracts {
            let (mut items, mut samples_noted, mut blocks) = (0, 0, [0_usize; 3]);
            for (row, workload) in samples {
                let Ok(called) = &row.outcome.result else {
                    continue;
                };
                let tract_notes = called
                    .notes
                    .iter()
                    .filter(|note| note.tract == tract.name)
                    .collect::<Vec<_>>();
                items += tract_notes.len();
                samples_noted += usize::from(!tract_notes.is_empty());
                let Ok(comparison) = row.comparison(Measure::HaplotypeLod) else {
                    continue;
                };
                for escape in workload
                    .escapes
                    .iter()
                    .filter(|escape| escape.class == BlockClass::LengthNote)
                {
                    let Some(block) = comparison.blocks.iter().find(|block| {
                        block.start == escape.first
                            && block.end == escape.last
                            && tract.length_only(block)
                    }) else {
                        continue;
                    };
                    let excess = tract.excess(block);
                    let convention = tract_notes.iter().flat_map(|note| &note.runs).any(|run| {
                        run.called
                            .zip(run.convention)
                            .is_some_and(|(called, convention)| {
                                i64::try_from(called).ok().map(|called| called - excess)
                                    == i64::try_from(convention).ok()
                            })
                    });
                    blocks[if convention {
                        0
                    } else if excess == 1 {
                        1
                    } else {
                        2
                    }] += 1;
                }
            }
            let _ = writeln!(
                text,
                "| {} | {items} | {samples_noted} | {} | {} | {} | {} |",
                tract.name,
                blocks.iter().sum::<usize>(),
                blocks[0],
                blocks[1],
                blocks[2],
            );
        }
    }

    /// The zones to check of `samples`: per sample (mean, median, 90th
    /// percentile, none), the items that touch more than one zone, the
    /// samples with a check in each zone, the unreviewed blocks by the zone
    /// of their first position.
    fn zone_section(&self, text: &mut String, samples: &[(&Row<'_>, &Workload)]) {
        let Some(zones) = &self.review.zones else {
            return;
        };
        if samples.is_empty() {
            return;
        }
        let mut counts = samples
            .iter()
            .map(|(_, workload)| workload.zones.len())
            .collect::<Vec<_>>();
        counts.sort_unstable();
        let rank = |quantile: f64| {
            // Nearest rank.
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let index = ((quantile * counts.len() as f64).ceil() as usize).max(1) - 1;
            counts[index.min(counts.len() - 1)]
        };
        let n = samples.len() as f64;
        let mean = counts.iter().sum::<usize>() as f64 / n;
        let none = counts.iter().filter(|&&count| count == 0).count();
        let items = samples
            .iter()
            .map(|(_, workload)| workload.checks)
            .sum::<usize>();
        let spanning = samples
            .iter()
            .map(|(_, workload)| workload.spanning)
            .sum::<usize>();
        let _ = write!(
            text,
            "\nZones to check (a check item counts for every zone it touches):\n\n\
             | Measure | Value |\n|---|---:|\n\
             | Zones per sample, mean | {mean:.2} |\n\
             | Zones per sample, median / 90th percentile | {} / {} |\n\
             | Samples with no zone | {none} ({:.2} %) |\n\
             | Check items touching more than one zone | {spanning} of {items} ({:.2} %) |\n",
            rank(0.5),
            rank(0.9),
            100.0 * none as f64 / n,
            if items > 0 {
                100.0 * spanning as f64 / items as f64
            } else {
                0.0
            },
        );
        text.push_str("\n| Zone | Samples with a check | Share |\n|---|---:|---:|\n");
        for instance in zones.instances() {
            let count = samples
                .iter()
                .filter(|(_, workload)| workload.zones.contains(&instance))
                .count();
            let _ = writeln!(
                text,
                "| {} | {count} | {:.2} % |",
                zones.label(instance),
                100.0 * count as f64 / n
            );
        }
        text.push_str(
            "\n| Zone kind | Samples with a check | Share | Unreviewed blocks |\n\
             |---|---:|---:|---:|\n",
        );
        for kind in Kind::ALL {
            let checked = samples
                .iter()
                .filter(|(_, workload)| workload.zones.iter().any(|&(_, found)| found == kind))
                .count();
            let blocks = samples
                .iter()
                .flat_map(|(_, workload)| &workload.escapes)
                .filter(|escape| escape.zone.is_some_and(|(_, found)| found == kind))
                .count();
            let _ = writeln!(
                text,
                "| {} | {checked} | {:.2} % | {blocks} |",
                kind.label(),
                100.0 * checked as f64 / n
            );
        }
    }

    /// `anchors.tsv`: every anchor candidate under each measure.
    fn anchors(&self) -> String {
        let mut text = String::from("sample_id\tstratum");
        for measure in Measure::ALL {
            let _ = write!(text, "\t{}", measure.column());
        }
        text.push('\n');
        for row in self.rows.iter().filter(|row| row.sample.anchor) {
            let _ = write!(text, "{}\t{}", row.sample.id, row.sample.stratum);
            for measure in Measure::ALL {
                let _ = write!(text, "\t{}", row.status(measure).label());
            }
            text.push('\n');
        }
        text
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
        self.class_section(&mut text);
        self.secondary_section(&mut text);
        self.anchor_section(&mut text);
        self.review_section(&mut text);

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

    /// Samples per discrepancy class and window, retired ones left out.
    fn class_section(&self, text: &mut String) {
        let current = self.rows.iter().filter(|row| !row.sample.retired);
        let mut classes = BTreeMap::<Discrepancy, (usize, BTreeMap<Region, usize>)>::new();
        for row in current.clone() {
            let discrepancy = row.discrepancy();
            let (samples, windows) = classes.entry(discrepancy).or_default();
            *samples += 1;
            for region in row.discrepancy_regions(discrepancy) {
                *windows.entry(region).or_default() += 1;
            }
        }
        text.push_str(
            "\n## Discrepancy classes\n\nEach sample not retired that failed or has a side \
             without an allele is `failed` or `not_comparable`; any other takes the first class \
             whose measure is not concordant: sequence (haplotype IUPAC-compatible), heteroplasmy \
             (haplotype), unanalysed or representation (primary), coverage (strict). A window column counts \
             the samples with a discordant position of their class in it.\n\n\
             | Class | Owner | Samples |",
        );
        for region in Region::ALL {
            let _ = write!(text, " {} |", region.label());
        }
        text.push_str("\n|---|---|---:|");
        text.push_str(&"---:|".repeat(Region::ALL.len()));
        text.push('\n');
        for discrepancy in Discrepancy::ALL {
            let (samples, windows) = classes.get(&discrepancy).cloned().unwrap_or_default();
            let _ = write!(
                text,
                "| {} | {} | {samples} |",
                discrepancy.label(),
                discrepancy.owner()
            );
            for region in Region::ALL {
                let _ = write!(
                    text,
                    " {} |",
                    windows.get(&region).copied().unwrap_or_default()
                );
            }
            text.push('\n');
        }
        self.cause_tables(text);
        self.empop_section(text);
        self.yield_section(text);
    }

    /// DNA's callable yield: per HV region, the truth's analysed positions
    /// DNA's callable intervals cover, retired samples left out.
    fn yield_section(&self, text: &mut String) {
        let yields = self
            .rows
            .iter()
            .filter(|row| !row.sample.retired)
            .filter_map(|row| row.callable_yield)
            .collect::<Vec<_>>();
        if yields.is_empty() {
            return;
        }
        text.push_str(
            "\nCallable yield: per HV region, the truth's analysed positions that DNA's \
             callable intervals cover.\n\n\
             | Region | Samples | Analysed | Callable | Share | Samples short |\n\
             |---|---:|---:|---:|---:|---:|\n",
        );
        for index in 0..3 {
            let region = yields
                .iter()
                .map(|values| values[index])
                .filter(|value| value.analysed > 0)
                .collect::<Vec<_>>();
            let analysed = region.iter().map(|value| value.analysed).sum::<u64>();
            let callable = region.iter().map(|value| value.callable).sum::<u64>();
            #[allow(clippy::cast_precision_loss)]
            let share = if analysed == 0 {
                0.0
            } else {
                callable as f64 / analysed as f64
            };
            let _ = writeln!(
                text,
                "| {} | {} | {analysed} | {callable} | {share:.4} | {} |",
                yields[0][index].region,
                region.len(),
                region
                    .iter()
                    .filter(|value| value.callable < value.analysed)
                    .count(),
            );
        }
    }

    /// The EMPOP queries (PROP-0007): how many samples have one, and the
    /// consensus decisions their ranges rest on.
    fn empop_section(&self, text: &mut String) {
        let queries = self
            .rows
            .iter()
            .filter(|row| !row.sample.retired)
            .filter_map(|row| row.outcome.result.as_ref().ok()?.empop.as_ref())
            .collect::<Vec<_>>();
        if queries.is_empty() {
            return;
        }
        let total =
            |count: fn(&Empop) -> usize| queries.iter().map(|empop| count(empop)).sum::<usize>();
        let samples =
            |count: fn(&Empop) -> usize| queries.iter().filter(|empop| count(empop) > 0).count();
        let _ = write!(
            text,
            "\nEMPOP queries: a range must cover only positions the consensus decided; a \
             tie the consensus resolved is exported as its decision.\n\n\
             | Item | Positions | Samples |\n|---|---:|---:|\n\
             | query | | {} |\n\
             | undecided position inside a range | {} | {} |\n\
             | resolved tie inside a range | {} | {} |\n",
            queries.len(),
            total(|empop| empop.undecided),
            samples(|empop| empop.undecided),
            total(|empop| empop.ties),
            samples(|empop| empop.ties),
        );
    }

    /// Samples per class and cause, retired ones left out, and per cause the
    /// samples with any finding of it: who has work to do.
    fn cause_tables(&self, text: &mut String) {
        let current = self
            .rows
            .iter()
            .zip(&self.causes)
            .filter(|(row, _)| !row.sample.retired);
        // Defects first, then deferred items, each in cause order; the
        // technology is always deferred.
        let columns = [
            (Cause::Dna, Scope::Defect),
            (Cause::Lab, Scope::Defect),
            (Cause::Undetermined, Scope::Defect),
            (Cause::Dna, Scope::Deferred),
            (Cause::Lab, Scope::Deferred),
            (Cause::Technology, Scope::Deferred),
        ];
        let mut headline = BTreeMap::<(Discrepancy, Cause, Scope), usize>::new();
        let mut routed = BTreeMap::<(Cause, Scope, &str), usize>::new();
        for (row, cause) in current {
            let Some((cause, scope, findings)) = cause else {
                continue;
            };
            *headline
                .entry((row.discrepancy(), *cause, *scope))
                .or_default() += 1;
            let evidence = findings
                .iter()
                .map(|finding| (finding.cause, finding.scope, finding.evidence))
                .collect::<BTreeSet<_>>();
            for key in evidence {
                *routed.entry(key).or_default() += 1;
            }
        }
        text.push_str(
            "\nCause of each sample's discrepancy, the first of dna, lab, technology, and \
             undetermined among its findings, a defect before an item deferred by design:\n\n\
             | Class | dna | lab | undetermined | dna (deferred) | lab (deferred) | technology (deferred) |\n\
             |---|---:|---:|---:|---:|---:|---:|\n",
        );
        for discrepancy in Discrepancy::ALL
            .into_iter()
            .filter(|discrepancy| *discrepancy != Discrepancy::None)
        {
            let _ = write!(text, "| {} |", discrepancy.label());
            for (cause, scope) in columns {
                let count = headline
                    .get(&(discrepancy, cause, scope))
                    .copied()
                    .unwrap_or_default();
                let _ = write!(text, " {count} |");
            }
            text.push('\n');
        }
        text.push_str(
            "\nSamples with a finding of each cause and evidence (a sample can need more \
             than one of them):\n\n| Cause | Scope | Evidence | Samples |\n|---|---|---|---:|\n",
        );
        for ((cause, scope, evidence), count) in routed {
            let _ = writeln!(
                text,
                "| {} | {} | {evidence} | {count} |",
                cause.label(),
                scope.label()
            );
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
                components: Vec::new(),
                tracts: Vec::new(),
                length_context: None,
                rejected: index % 2,
                alleles,
                not_recorded: BTreeSet::new(),
                empop: None,
                checks: Vec::new(),
                check_labels: Vec::new(),
                notes: Vec::new(),
                signals: lod::Signals::default(),
                loci: cause::loci::Loci::default(),
                intervals: if index % 5 == 0 {
                    vec![(73, 340)]
                } else {
                    vec![(73, 340), (16024, 16365)]
                },
                identity: identity(),
            })
        };
        let call =
            |sample: &Sample, traces: &[PathBuf]| -> crate::runner::Ran { called(sample, traces) };
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
                &[],
            );
            Ok(Report::new(&rows, &check, (&Forms::default(), ""))
                .files(&[("panel", "synthetic".to_owned())]))
        };
        let one = outputs(1)?;
        assert_eq!(one, outputs(8)?);
        crate::contract::check_files(one.iter().map(|(name, text)| (*name, text.as_str())))?;
        let names = one.iter().map(|(name, _)| *name).collect::<Vec<_>>();
        assert_eq!(
            names,
            [
                "scores.tsv",
                "estimates.tsv",
                "anchors.tsv",
                "windows.tsv",
                "unreviewed.tsv",
                "summary.md"
            ]
        );
        let scores = &one[0].1;
        let expected = [
            // The IUPAC site lies outside DNA's callable intervals but inside
            // the truth's: concordant, but missed on the strict measure.
            "S00\tb1\tb1\thard\t4.5\t4.5\t1\t0\tconcordant\tY\t2\t3\t\t16189Y\t\t\t0\t\t\
             discordant\tN\t16189Y\tconcordant\tY\tconcordant\tY\t\tconcordant\tY\t\tconcordant\tY\t\t1\tcoverage\thv1_polyc\tlab\tdefect\tlab:no_read@16189\t\t\t\t\tHV1=0/342 HV2=268/268 HV3=0/0\t2\t0\t",
            // The major allele at the IUPAC site: compatible with Y.
            "S01\tb1\tb1\tclean\t10\t10\t0\t0\tdiscordant\tN\t3\t3\t16189C\t16189Y\t\
             shared_position,dna_only,truth_only\t16189\t1\thv1_polyc\tdiscordant\tN\t\t\
             concordant\tY\tdiscordant\tN\t16189-16189:C/Y\tconcordant\tY\t\tconcordant\tY\t\t1\theteroplasmy\thv1_polyc\tundetermined\tdefect\tundetermined:one_strand@16189\t\t\t\t\tHV1=342/342 HV2=268/268 HV3=0/0\t2\t1\t",
            // Duplicate trace content is removed before DNA runs.
            "S02\tb1\tb1\tclean\t10\t10\t0\t0\tfailed\t\t\t\t\t\t\t\t\t\tfailed\t\t\tfailed\t\tfailed\t\t\tfailed\t\t\tfailed\t\t\t\tfailed\t\tdna\tdefect\tdna:dna_failure\t\t\t\t\t\t\
             2\t\tdna: 2 traces",
            "S03\tb1\tb1\thard\t4.5\t4.5\t0\t0\tnot_applicable\tN/A\t0\t3\t\t263G 16189Y\t\t\t0\t\t\
             not_applicable\tN/A\t\tnot_applicable\tN/A\tnot_applicable\tN/A\t\t\
             not_applicable\tN/A\t\tnot_applicable\tN/A\t\t1\tnot_comparable\t\tundetermined\tdefect\tundetermined:not_comparable\t\t\t\t\tHV1=342/342 HV2=268/268 HV3=0/0\t2\t1\t",
            "S04\tb1\tb1\tclean\t10\t10\t0\t0\tdiscordant\tN\t3\t3\t150T\t16189Y\t\
             dna_only,truth_only\t150 16189\t0\thv1_polyc,other\tdiscordant\tN\t\tdiscordant\tN\t\
             discordant\tN\t150-150:T/C 16189-16189:T/Y\tdiscordant\tN\t150-150:T/C\tdiscordant\tN\t150-150:T/C\t1\tsequence\tother\tundetermined\tdefect\tundetermined:adjudicate\t\t\t\t\tHV1=342/342 HV2=268/268 HV3=0/0\t2\t0\t",
            "S08\tb1\tb2\tclean\t10\t10\t0\t0\tnot_applicable\tN/A\t0\t3\t\t263G 16189Y\t\t\t0\t\t\
             not_applicable\tN/A\t\tnot_applicable\tN/A\tnot_applicable\tN/A\t\t\
             not_applicable\tN/A\t\tnot_applicable\tN/A\t\t1\tnot_comparable\t\tundetermined\tdefect\tundetermined:not_comparable\t\t\t\t\tHV1=342/342 HV2=268/268 HV3=0/0\t2\t0\t",
            "S09\tb1\t\thard\t4.5\t4.5\t0\t0\tfailed\t\t\t\t\t\t\t\t\t\tfailed\t\t\tfailed\t\tfailed\t\t\tfailed\t\t\tfailed\t\t\t\tfailed\t\tlab\tdefect\tlab:truth_missing\t\t\t\t\t\t\
             2\t1\ttruth_missing: the truth has no record of the sample",
        ];
        for line in expected {
            // The review workload columns follow the failure column.
            assert!(scores.contains(&format!("{line}\t")), "{line}");
        }

        let anchors = &one[2].1;
        assert!(anchors.contains(
            "S00\thard\tconcordant\tdiscordant\tconcordant\tconcordant\tconcordant\tconcordant\n"
        ));
        let estimates = &one[1].1;
        assert!(estimates.contains("corpus\tstrict\tnot_concordant\tall\t"));
        assert!(estimates.contains("pool\tall\tfailed\tall\t"));
        let summary = &one[5].1;
        for line in [
            "| without a truth record | 1 |",
            "| truth batch differs from the panel batch | 1 |",
            "| coverage | core callability, or the lab's assay | 6 | 0 | 6 | 0 | 0 | 0 |",
            "| sequence | core, or the truth (adjudication) | 5 | 0 | 0 | 0 | 0 | 5 |",
            "| not_comparable | the input, or coverage | 6 | 0 | 0 | 0 | 0 | 0 |",
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
    fn a_changed_truth_is_reported_and_listed() -> crate::contract::Checked {
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
        let files = Report::new(&[], &check, (&Forms::default(), "")).files(&[]);
        let names = files.iter().map(|(name, _)| *name).collect::<Vec<_>>();
        assert_eq!(
            names,
            [
                "scores.tsv",
                "estimates.tsv",
                "anchors.tsv",
                "windows.tsv",
                "truth_changes.tsv",
                "unreviewed.tsv",
                "summary.md"
            ]
        );
        assert!(
            files[4]
                .1
                .contains(&format!("S1\tclean\t0\t{}\t{}\n", digest("c"), digest("d")))
        );
        crate::contract::check_table("truth_changes.tsv", &files[4].1)?;
        let summary = &files[6].1;
        assert!(summary.starts_with("# Evaluation run\n\n**Truth changed:** the truth of 1 panel"));
        assert!(summary.contains("| truth changed since the panel was built | 1 |"));
        let unchanged = TruthCheck {
            recorded: Some(digest("b")),
            changed: Vec::new(),
            ..check
        };
        let files = Report::new(&[], &unchanged, (&Forms::default(), "")).files(&[]);
        assert!(files.iter().all(|(name, _)| *name != "truth_changes.tsv"));
        assert!(files[5].1.contains("Panel truth: unchanged"));
        // The header digest differs while no row does: the file is still
        // written, and the notice says why it lists nothing.
        let header_only = TruthCheck {
            recorded: Some(digest("a")),
            changed: Vec::new(),
            ..unchanged
        };
        let files = Report::new(&[], &header_only, (&Forms::default(), "")).files(&[]);
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
        Ok(())
    }
}
