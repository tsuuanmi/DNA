//! Scores the draft triage against the truth (PROP-0004 phase 0, stage C):
//! A/B/C at sample and region level, leaks of passed samples and regions
//! under every measure, catches and burden per detector, and a sweep of
//! policies as a risk–coverage table.
//!
//! Retired samples are left out. A sample DNA failed on is its own class; a
//! sample without a comparison (no truth) counts in A/B/C but not in leaks
//! or catches.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

use super::detect::{
    Detector, Footprints, Outcome, Policy, Proxy, SampleDecision, TraceFinding, decide,
};
use crate::compare::{self, Allele, Comparison, Measure, Region as Window};
use crate::estimate::{self, Estimate, Observation};
use crate::evidence::{Cause, Evidence, TraceState};
use crate::report::{Row, Status, Weighting};

/// One-sided level of the Clopper–Pearson upper bound: a 95 % bound.
const ALPHA: f64 = 0.05;

/// The triage class of a sample.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Class {
    Decided(Outcome),
    /// DNA failed on the sample.
    Failed,
}

impl Class {
    const ALL: [Self; 4] = [
        Self::Decided(Outcome::Pass),
        Self::Decided(Outcome::Review),
        Self::Decided(Outcome::Rerun),
        Self::Failed,
    ];

    const fn label(self) -> &'static str {
        match self {
            Self::Decided(outcome) => outcome.label(),
            Self::Failed => "failed",
        }
    }
}

/// The policies of the sweep: the draft first, then one change each.
pub(crate) fn sweep() -> Vec<Policy> {
    let draft = Policy::draft();
    let with = |name: &str, change: &dyn Fn(&mut Policy)| {
        let mut policy = Policy {
            name: name.to_owned(),
            ..draft.clone()
        };
        change(&mut policy);
        policy
    };
    vec![
        draft.clone(),
        with("short_off", &|policy| policy.short_calls = 0),
        with("short_50", &|policy| policy.short_calls = 50),
        with("short_150", &|policy| policy.short_calls = 150),
        with("short_200", &|policy| policy.short_calls = 200),
        with("short_300", &|policy| policy.short_calls = 300),
        with("proxy_any", &|policy| policy.proxy = Proxy::Any),
        with("one_strand_off", &|policy| policy.one_strand = false),
        with("masked_off", &|policy| {
            policy.causes.remove(&Cause::Masked);
        }),
        with("masked_read_on", &|policy| {
            policy.causes.insert(Cause::MaskedRead);
        }),
        with("n_contested_only", &|policy| {
            policy.causes = [Cause::N, Cause::Contested].into();
        }),
        with("unreliable_off", &|policy| policy.causes.clear()),
        with("one_strand_masked_off", &|policy| {
            policy.one_strand = false;
            policy.causes.remove(&Cause::Masked);
        }),
        with("gap_only", &|policy| {
            policy.causes.clear();
            policy.one_strand = false;
        }),
    ]
}

/// The triage of every panel sample under every policy of the sweep, measured
/// once: the reported (first) policy's decisions and discordant positions, and
/// each policy's A/B/C and leaks.
pub(crate) struct Triage<'a> {
    rows: &'a [Row<'a>],
    footprints: Footprints,
    policies: Vec<Policy>,
    /// Per row, the reported policy's decision; `None` without triage
    /// evidence.
    decisions: Vec<Option<SampleDecision>>,
    /// Per row, the reported policy's class.
    classes: Vec<Class>,
    /// Per policy, its measurement.
    measurements: Vec<Measurement>,
    /// Region labels, from the first sample with region evidence.
    labels: Vec<String>,
    /// Every discordant position of the scored, compared rows with region
    /// evidence, under every measure.
    cases: Vec<Case>,
    /// The indices into `cases` of each row.
    by_row: BTreeMap<usize, Vec<usize>>,
}

/// The triage evidence of a row, when the run took it.
fn evidence<'a>(row: &'a Row<'a>) -> Option<&'a Evidence> {
    row.outcome.triage.as_ref()
}

/// The class of a row: its decision, unless DNA failed on a sample that had
/// an admitted read. A sample without an admitted read is decided (`rerun`)
/// although DNA reports it as failed.
fn class(row: &Row<'_>, decision: Option<&SampleDecision>) -> Class {
    match decision {
        Some(decision) if row.outcome.result.is_ok() || !decision.admitted => {
            Class::Decided(decision.outcome)
        }
        _ => Class::Failed,
    }
}

impl<'a> Triage<'a> {
    /// Decides every row under every policy; the first policy is the one the
    /// per-sample outputs report. Primer footprints are learned from the
    /// active (not retired) samples.
    pub(crate) fn new(rows: &'a [Row<'a>], policies: Vec<Policy>) -> Self {
        let footprints = Footprints::learn(
            rows.iter()
                .filter(|row| !row.sample.retired)
                .filter_map(evidence),
        );
        let labels: Vec<String> = rows
            .iter()
            .filter_map(evidence)
            .find(|evidence| !evidence.regions.is_empty())
            .map(|evidence| {
                evidence
                    .regions
                    .iter()
                    .map(crate::evidence::Region::label)
                    .collect()
            })
            .unwrap_or_default();
        let mut reported = Vec::new();
        let mut measurements = Vec::new();
        for (index, policy) in policies.iter().enumerate() {
            let decisions = rows
                .iter()
                .map(|row| evidence(row).map(|evidence| decide(evidence, policy, &footprints)))
                .collect::<Vec<_>>();
            measurements.push(Measurement::of(rows, &decisions, labels.len()));
            if index == 0 {
                reported = decisions;
            }
        }
        let classes = rows
            .iter()
            .zip(&reported)
            .map(|(row, decision)| class(row, decision.as_ref()))
            .collect::<Vec<_>>();
        let cases = cases(rows, &reported, &classes);
        let mut by_row = BTreeMap::<usize, Vec<usize>>::new();
        for (index, case) in cases.iter().enumerate() {
            by_row.entry(case.row).or_default().push(index);
        }
        Self {
            rows,
            footprints,
            policies,
            decisions: reported,
            classes,
            measurements,
            labels,
            cases,
            by_row,
        }
    }

    /// The number of policies decided.
    pub(crate) fn policies(&self) -> usize {
        self.policies.len()
    }

    /// The policy the per-sample outputs report.
    pub(crate) fn policy(&self) -> Option<&Policy> {
        self.policies.first()
    }

    /// Every output file, by name.
    pub(crate) fn files(&self) -> Vec<(&'static str, String)> {
        vec![
            ("triage.tsv", self.samples()),
            ("triage_regions.tsv", self.regions()),
            ("triage_traces.tsv", self.traces()),
            ("triage_leaks.tsv", self.leak_cases()),
            ("triage_detectors.tsv", self.detectors()),
            ("triage_positions.tsv", self.positions()),
            ("triage_sweep.tsv", self.sweep_table()),
            ("triage.md", self.summary()),
        ]
    }

    /// The rows with their index, decision, and class under the reported
    /// policy.
    fn all(
        &self,
    ) -> impl Iterator<Item = (usize, &'a Row<'a>, Option<&SampleDecision>, Class)> + '_ {
        self.rows
            .iter()
            .zip(&self.decisions)
            .zip(&self.classes)
            .enumerate()
            .map(|(index, ((row, decision), class))| (index, row, decision.as_ref(), *class))
    }

    /// [`Triage::all`] without the retired rows.
    fn scored(
        &self,
    ) -> impl Iterator<Item = (usize, &'a Row<'a>, Option<&SampleDecision>, Class)> + '_ {
        self.all().filter(|(_, row, _, _)| !row.sample.retired)
    }

    /// The reported policy's measurement.
    fn reported(&self) -> Option<&Measurement> {
        self.measurements.first()
    }

    /// `triage.tsv`: one row per sample under the reported policy.
    fn samples(&self) -> String {
        let mut text = String::from(
            "sample_id\tstratum\tweight_corpus\tretired\toutcome\tregions\tfindings\t\
             poor_traces\trerun_traces\tstatus\tstrict_status\tiupac_compatible_status\n",
        );
        for (_, row, decision, class) in self.all() {
            let traces = evidence(row).map_or(&[][..], |evidence| &evidence.traces[..]);
            let (regions, findings, poor, rerun) =
                decision.map_or_else(Default::default, |decision| {
                    let regions = self
                        .labels
                        .iter()
                        .zip(&decision.regions)
                        .map(|(label, region)| format!("{label}={}", region.outcome.label()))
                        .collect::<Vec<_>>()
                        .join(",");
                    let findings = self
                        .labels
                        .iter()
                        .zip(&decision.regions)
                        .flat_map(|(label, region)| {
                            region.findings.iter().filter(|(_, count)| **count > 0).map(
                                move |(detector, count)| {
                                    format!("{label}:{}={count}", detector.label())
                                },
                            )
                        })
                        .collect::<Vec<_>>()
                        .join(",");
                    let poor = traces
                        .iter()
                        .zip(&decision.traces)
                        .filter(|(_, findings)| !findings.is_empty())
                        .map(|(trace, findings)| {
                            format!("{}={}", trace.name, finding_labels(findings))
                        })
                        .collect::<Vec<_>>()
                        .join(",");
                    let rerun = decision
                        .rerun
                        .iter()
                        .filter_map(|index| traces.get(*index).map(|trace| trace.name.as_str()))
                        .collect::<Vec<_>>()
                        .join(",");
                    (regions, findings, poor, rerun)
                });
            let _ = writeln!(
                text,
                "{}\t{}\t{}\t{}\t{}\t{regions}\t{findings}\t{poor}\t{rerun}\t{}\t{}\t{}",
                row.sample.id,
                row.sample.stratum,
                row.sample.weight_corpus.text,
                u8::from(row.sample.retired),
                class.label(),
                row.status(Measure::Primary).label(),
                row.status(Measure::Strict).label(),
                row.status(Measure::IupacCompatible).label(),
            );
        }
        text
    }

    /// `triage_regions.tsv`: one row per sample and region.
    fn regions(&self) -> String {
        let detectors = Detector::all();
        let mut text = String::from("sample_id\tstratum\tregion\toutcome");
        for detector in &detectors {
            let _ = write!(text, "\t{}", detector.label());
        }
        text.push_str(
            "\trerun_traces\tprimary_positions\tstrict_positions\tiupac_compatible_positions\n",
        );
        for (_, row, decision, _) in self.all() {
            let (Some(decision), Some(evidence)) = (decision, evidence(row)) else {
                continue;
            };
            for (index, (label, region)) in self.labels.iter().zip(&decision.regions).enumerate() {
                let _ = write!(
                    text,
                    "{}\t{}\t{label}\t{}",
                    row.sample.id,
                    row.sample.stratum,
                    region.outcome.label()
                );
                for detector in &detectors {
                    let _ = write!(
                        text,
                        "\t{}",
                        region.findings.get(detector).copied().unwrap_or(0)
                    );
                }
                let rerun = region
                    .rerun
                    .iter()
                    .filter_map(|trace| {
                        evidence.traces.get(*trace).map(|trace| trace.name.as_str())
                    })
                    .collect::<Vec<_>>()
                    .join(",");
                let _ = write!(text, "\t{rerun}");
                for measure in Measure::ALL {
                    let positions = row
                        .comparison(measure)
                        .map(|comparison| {
                            compare::sorted_positions(&comparison.positions)
                                .into_iter()
                                .filter(|position| {
                                    compare::anchor(position)
                                        .and_then(|anchor| evidence.region_of(anchor))
                                        == Some(index)
                                })
                                .collect::<Vec<_>>()
                                .join(" ")
                        })
                        .unwrap_or_default();
                    let _ = write!(text, "\t{positions}");
                }
                text.push('\n');
            }
        }
        text
    }

    /// `triage_traces.tsv`: one row per input trace.
    fn traces(&self) -> String {
        let mut text =
            String::from("sample_id\ttrace\tprimer\tstate\tcallable_span\tfindings\trerun_for\n");
        for (_, row, decision, _) in self.all() {
            let (Some(decision), Some(evidence)) = (decision, evidence(row)) else {
                continue;
            };
            for (index, (trace, findings)) in
                evidence.traces.iter().zip(&decision.traces).enumerate()
            {
                let state = match &trace.state {
                    TraceState::Admitted { strand, .. } => format!("admitted:{}", strand.label()),
                    TraceState::TooShort => "too_short".to_owned(),
                    TraceState::Weak => "weak".to_owned(),
                    TraceState::Unplaced(label) => format!("unplaced:{label}"),
                    TraceState::Rejected => "rejected".to_owned(),
                };
                let rerun_for = if decision.admitted {
                    self.labels
                        .iter()
                        .zip(&decision.regions)
                        .filter(|(_, region)| region.rerun.contains(&index))
                        .map(|(label, _)| label.as_str())
                        .collect::<Vec<_>>()
                        .join(",")
                } else if decision.rerun.contains(&index) {
                    NO_ADMITTED_READ.to_owned()
                } else {
                    String::new()
                };
                let _ = writeln!(
                    text,
                    "{}\t{}\t{}\t{state}\t{}\t{}\t{rerun_for}",
                    row.sample.id,
                    trace.name,
                    trace.primer.as_deref().unwrap_or(""),
                    trace.span.map(|span| span.to_string()).unwrap_or_default(),
                    finding_labels(findings),
                );
            }
        }
        text
    }

    /// `triage_leaks.tsv`: discordant positions inside passed regions.
    fn leak_cases(&self) -> String {
        let mut text = String::from(
            "sample_id\tstratum\tsample_outcome\tmeasure\tregion\tposition\twindow\tkind\t\
             class\tdna_alleles\ttruth_alleles\tmarked_by\n",
        );
        for case in self.cases.iter().filter(|case| case.leaks()) {
            let row = &self.rows[case.row];
            let _ = writeln!(
                text,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                row.sample.id,
                row.sample.stratum,
                case.sample.label(),
                case.measure.label(),
                self.region_label(case.region.map(|(index, _)| index)),
                case.position,
                Window::of(&case.position).label(),
                case.detail.kind,
                case.detail.class,
                compare::texts(&case.detail.dna),
                compare::texts(&case.detail.truth),
                case.flags
                    .iter()
                    .map(|detector| detector.label())
                    .collect::<Vec<_>>()
                    .join(","),
            );
        }
        text
    }

    /// The label of a region index; `-` outside every region.
    fn region_label(&self, index: Option<usize>) -> &str {
        index
            .and_then(|index| self.labels.get(index))
            .map_or("-", String::as_str)
    }

    /// Per region detector, under the reported policy's findings, over the
    /// scored, compared samples with a decision: samples flagged, catches,
    /// burden, and positional hits; for an applied detector the catches and
    /// burden no other applied detector shares; for a detector the policy does
    /// not apply, what applying it would add.
    fn detector_rows(&self) -> Vec<DetectorRow> {
        let Some(policy) = self.policy() else {
            return Vec::new();
        };
        let compared = self
            .scored()
            .filter_map(|(_, row, decision, _)| Some((row, decision?)))
            .filter(|(row, _)| row.comparisons.is_ok())
            .collect::<Vec<_>>();
        let flagged = |decision: &SampleDecision, detector: Detector| {
            decision.regions.iter().any(|region| {
                region
                    .findings
                    .get(&detector)
                    .is_some_and(|count| *count > 0)
            })
        };
        Detector::all()
            .into_iter()
            .map(|detector| {
                let applied = policy.applies(detector);
                let others = Detector::all()
                    .into_iter()
                    .filter(|other| *other != detector && policy.applies(*other))
                    .collect::<Vec<_>>();
                let mut row = DetectorRow::new(detector, applied);
                for (sample, decision) in &compared {
                    if !flagged(decision, detector) {
                        continue;
                    }
                    row.samples += 1;
                    for region in &decision.regions {
                        let count = region.findings.get(&detector).copied().unwrap_or(0);
                        row.regions += usize::from(count > 0);
                        row.positions += count;
                    }
                    // Flagged by no other applied detector: for an applied
                    // detector the samples only it holds back; for one the
                    // policy does not apply, the samples it would newly hold
                    // back.
                    let alone = !others.iter().any(|other| flagged(decision, *other));
                    for measure in Measure::ALL {
                        if sample.status(measure) != Status::Concordant {
                            row.catches[measure.index()] += 1;
                            row.alone[measure.index()] += usize::from(alone);
                        }
                    }
                    let concordant = sample.status(Measure::Primary) == Status::Concordant;
                    row.burden += usize::from(concordant);
                    row.alone_burden += usize::from(concordant && alone);
                }
                row.hits = self
                    .cases
                    .iter()
                    .filter(|case| {
                        case.measure == Measure::Primary && case.flags.contains(&detector)
                    })
                    .count();
                row
            })
            .collect()
    }

    /// `triage_positions.tsv`: per region position, the analysed samples
    /// each region detector marks it in, so findings can be located.
    fn positions(&self) -> String {
        let detectors = Detector::all();
        let mut text = String::from("region\tposition\tsamples");
        for detector in &detectors {
            let _ = write!(text, "\t{}", detector.label());
        }
        text.push('\n');
        let analysed = self
            .scored()
            .filter_map(|(_, row, _, _)| evidence(row))
            .filter(|evidence| !evidence.regions.is_empty())
            .collect::<Vec<_>>();
        let Some(first) = analysed.first() else {
            return text;
        };
        for (index, region) in first.regions.iter().enumerate() {
            for offset in 0..region.positions.len() {
                let mut counts = vec![0_usize; detectors.len()];
                let mut samples = 0;
                for evidence in &analysed {
                    let Some(at) = evidence
                        .regions
                        .get(index)
                        .and_then(|region| region.positions.get(offset))
                    else {
                        continue;
                    };
                    samples += 1;
                    for (count, detector) in counts.iter_mut().zip(&detectors) {
                        *count += usize::from(detector.marks(*at));
                    }
                }
                let _ = write!(
                    text,
                    "{}\t{}\t{samples}",
                    region.label(),
                    region.start + offset
                );
                for count in counts {
                    let _ = write!(text, "\t{count}");
                }
                text.push('\n');
            }
        }
        text
    }

    /// `triage_detectors.tsv`: `unique_*` and `sole_burden` for applied
    /// detectors, `if_applied_*` for the others; the other set is `-`.
    fn detectors(&self) -> String {
        let mut text = String::from(
            "detector\tapplied\tsamples\tregions\tpositions\tcatches_primary\tcatches_strict\t\
             catches_iupac_compatible\tburden\tunique_primary\tunique_strict\t\
             unique_iupac_compatible\tsole_burden\tif_applied_catches_primary\t\
             if_applied_catches_strict\tif_applied_catches_iupac_compatible\t\
             if_applied_burden\tposition_hits_primary\n",
        );
        for row in self.detector_rows() {
            let alone = format!(
                "{}\t{}\t{}\t{}",
                row.alone[0], row.alone[1], row.alone[2], row.alone_burden
            );
            let none = "-\t-\t-\t-";
            let (unique, added) = if row.applied {
                (alone.as_str(), none)
            } else {
                (none, alone.as_str())
            };
            let _ = writeln!(
                text,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{unique}\t{added}\t{}",
                row.detector.label(),
                u8::from(row.applied),
                row.samples,
                row.regions,
                row.positions,
                row.catches[0],
                row.catches[1],
                row.catches[2],
                row.burden,
                row.hits,
            );
        }
        text
    }

    /// `triage_sweep.tsv`: A/B/C and leaks per policy.
    fn sweep_table(&self) -> String {
        let mut text = String::from(
            "policy\tparameters\tsamples\tpass\treview\trerun\tfailed\tpass_rate\t\
             pass_corpus\tpass_corpus_ci95\treview_corpus\trerun_corpus\tfailed_corpus",
        );
        for measure in Measure::ALL {
            let label = measure.label();
            let _ = write!(
                text,
                "\tpassed_{label}\tleaks_{label}\tleak_rate_{label}\tleak_upper95_{label}\t\
                 leak_rate_corpus_{label}"
            );
        }
        text.push('\n');
        for (policy, measured) in self.policies.iter().zip(&self.measurements) {
            let corpus = |class: Class| {
                measured
                    .estimate(class)
                    .map_or((0.0, 0.0), |estimate| (estimate.rate, estimate.ci95))
            };
            let _ = write!(
                text,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
                policy.name,
                policy.parameters(),
                measured.samples,
                measured.count(Class::Decided(Outcome::Pass)),
                measured.count(Class::Decided(Outcome::Review)),
                measured.count(Class::Decided(Outcome::Rerun)),
                measured.count(Class::Failed),
                measured.share(Class::Decided(Outcome::Pass)),
                corpus(Class::Decided(Outcome::Pass)).0,
                corpus(Class::Decided(Outcome::Pass)).1,
                corpus(Class::Decided(Outcome::Review)).0,
                corpus(Class::Decided(Outcome::Rerun)).0,
                corpus(Class::Failed).0,
            );
            for leak in &measured.leaks {
                let _ = write!(
                    text,
                    "\t{}\t{}\t{:.4}\t{:.4}\t{:.4}",
                    leak.passed,
                    leak.leaks,
                    leak.rate(),
                    upper_bound(leak.leaks, leak.passed),
                    leak.weighted_rate()
                );
            }
            text.push('\n');
        }
        text
    }

    /// `triage.md`: the triage in aggregate.
    fn summary(&self) -> String {
        let mut text = String::from(
            "# Draft triage (PROP-0004 phase 0)\n\nResearch-grade stage-1 detectors over DNA's \
             evidence, scored against the truth; not a DNA contract and not production triage. \
             Retired samples are left out.\n\n",
        );
        if let Some(policy) = self.policy() {
            let _ = writeln!(
                text,
                "Reported policy `{}`: {}.",
                policy.name,
                policy.parameters()
            );
        }
        let Some(measured) = self.reported() else {
            return text;
        };
        text.push_str(
            "\n## Samples\n\nA sample without an admitted read is `rerun` with every trace; \
             `failed` is a sample DNA could not analyse or call.\n\n\
             | Outcome | Samples | Share | Corpus-weighted | 95 % CI |\n\
             |---|---:|---:|---:|---:|\n",
        );
        for class in Class::ALL {
            let (rate, ci95) = measured
                .estimate(class)
                .map_or((0.0, 0.0), |estimate| (estimate.rate, estimate.ci95));
            let _ = writeln!(
                text,
                "| {} | {} | {:.4} | {rate:.4} | ± {ci95:.4} |",
                class.label(),
                measured.count(class),
                measured.share(class),
            );
        }
        strata_section(&mut text, measured);
        self.leak_section(&mut text, measured);
        self.region_section(&mut text, measured);
        self.class_section(&mut text);
        self.detector_section(&mut text);
        self.trace_section(&mut text);
        self.sweep_section(&mut text);
        text
    }

    fn leak_section(&self, text: &mut String, measured: &Measurement) {
        text.push_str(
            "\n## Leaks\n\nA leak is a passed sample that is not concordant. The upper bound is \
             the one-sided 95 % Clopper–Pearson bound on the leak rate among passed samples; \
             the corpus rate weights passed samples by `weight_corpus`.\n\n\
             | Measure | Passed | Leaks | Rate | Upper 95 % | Corpus rate |\n\
             |---|---:|---:|---:|---:|---:|\n",
        );
        for measure in Measure::ALL {
            let leak = measured.leaks[measure.index()];
            let _ = writeln!(
                text,
                "| {} | {} | {} | {:.4} | {:.4} | {:.4} |",
                measure.label(),
                leak.passed,
                leak.leaks,
                leak.rate(),
                upper_bound(leak.leaks, leak.passed),
                leak.weighted_rate()
            );
        }
        let mut composition = BTreeMap::<(Measure, &str), usize>::new();
        for (index, row, _, class) in self.scored() {
            if class != Class::Decided(Outcome::Pass) || row.comparisons.is_err() {
                continue;
            }
            let cases = self
                .by_row
                .get(&index)
                .map_or(&[][..], Vec::as_slice)
                .iter()
                .map(|case| &self.cases[*case])
                .collect::<Vec<_>>();
            for measure in Measure::ALL {
                if row.status(measure) == Status::Concordant {
                    continue;
                }
                let at = cases
                    .iter()
                    .filter(|case| case.measure == measure)
                    .collect::<Vec<_>>();
                let class = if at.is_empty() {
                    "no_position"
                } else if at.iter().all(|case| case.detail.class == "iupac") {
                    "iupac"
                } else if at.iter().all(|case| {
                    case.detail.class == "iupac" || Window::of(&case.position) != Window::Other
                }) {
                    "repeat_window"
                } else {
                    "other"
                };
                *composition.entry((measure, class)).or_default() += 1;
            }
        }
        text.push_str(
            "\nLeaked samples by what differs: `iupac` only truth IUPAC alleles, `repeat_window` \
             only IUPAC alleles or positions in a repeat window (`hv1_polyc`, `hv2_polyc`, \
             `hv3_polyc`, `ac_repeat`), `other` anything else, `no_position` not concordant \
             without a discordant position (`N/A`).\n\n\
             | Measure | iupac | repeat_window | other | no_position |\n|---|---:|---:|---:|---:|\n",
        );
        for measure in Measure::ALL {
            let count = |class: &str| composition.get(&(measure, class)).copied().unwrap_or(0);
            let _ = writeln!(
                text,
                "| {} | {} | {} | {} | {} |",
                measure.label(),
                count("iupac"),
                count("repeat_window"),
                count("other"),
                count("no_position")
            );
        }
    }

    fn region_section(&self, text: &mut String, measured: &Measurement) {
        let mut leaks = BTreeMap::<(Measure, Option<usize>), usize>::new();
        for case in self.cases.iter().filter(|case| case.leaks()) {
            let region = case.region.map(|(index, _)| index);
            *leaks.entry((case.measure, region)).or_default() += 1;
            *leaks.entry((case.measure, None)).or_default() += 1;
        }
        text.push_str(
            "\n## Regions\n\nOutcomes of the regions of decided samples with region evidence, \
             as counts and as shares weighted by the samples' `weight_corpus`, and discordant \
             positions inside passed regions (region leaks) per measure.\n\n\
             | Region | Pass | Review | Rerun | Pass corpus | Review corpus | Rerun corpus | Leaks primary | Leaks strict | Leaks IUPAC-compatible |\n\
             |---|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n",
        );
        let mut total = Tally::default();
        let lines = self
            .labels
            .iter()
            .zip(&measured.regions)
            .enumerate()
            .map(|(index, (label, tally))| (label.as_str(), tally, Some(index)))
            .collect::<Vec<_>>();
        for (_, tally, _) in &lines {
            total.add(tally);
        }
        for (label, tally, index) in lines.into_iter().chain([("all", &total, None)]) {
            let leak = |measure: Measure| leaks.get(&(measure, index)).copied().unwrap_or(0);
            let _ = writeln!(
                text,
                "| {label} | {} | {} | {} | {:.4} | {:.4} | {:.4} | {} | {} | {} |",
                tally.count(Outcome::Pass),
                tally.count(Outcome::Review),
                tally.count(Outcome::Rerun),
                tally.weighted(Outcome::Pass),
                tally.weighted(Outcome::Review),
                tally.weighted(Outcome::Rerun),
                leak(Measure::Primary),
                leak(Measure::Strict),
                leak(Measure::IupacCompatible),
            );
        }
    }

    fn class_section(&self, text: &mut String) {
        let mut classes =
            BTreeMap::<(Measure, &str, &str, &str, &str), (usize, BTreeSet<usize>)>::new();
        for case in self.cases.iter().filter(|case| case.leaks()) {
            let entry = classes
                .entry((
                    case.measure,
                    self.region_label(case.region.map(|(index, _)| index)),
                    Window::of(&case.position).label(),
                    case.detail.kind,
                    case.detail.class,
                ))
                .or_default();
            entry.0 += 1;
            entry.1.insert(case.row);
        }
        text.push_str(
            "\n## Leak classes\n\nDiscordant positions inside passed regions, by measure, region, \
             window, kind, and allele class (`triage_leaks.tsv` lists them).\n\n\
             | Measure | Region | Window | Kind | Class | Positions | Samples |\n\
             |---|---|---|---|---|---:|---:|\n",
        );
        let mut sorted = classes.into_iter().collect::<Vec<_>>();
        sorted.sort_by(|(left_key, left), (right_key, right)| {
            (left_key.0, std::cmp::Reverse(left.0), left_key).cmp(&(
                right_key.0,
                std::cmp::Reverse(right.0),
                right_key,
            ))
        });
        for ((measure, region, window, kind, class), (positions, samples)) in sorted {
            let _ = writeln!(
                text,
                "| {} | {region} | {window} | {kind} | {class} | {positions} | {} |",
                measure.label(),
                samples.len()
            );
        }
    }

    fn detector_section(&self, text: &mut String) {
        let discordant = |measure: Measure| {
            self.scored()
                .filter(|(_, row, decision, _)| decision.is_some() && row.comparisons.is_ok())
                .filter(|(_, row, _, _)| row.status(measure) != Status::Concordant)
                .count()
        };
        let rows = self.detector_rows();
        let _ = write!(
            text,
            "\n## Detectors\n\nPer region detector, over analysed and compared samples: samples \
             it flags, not-concordant samples among them (catches), concordant samples it flags \
             (burden), and primary discordant positions it marks itself (hits). Not concordant: \
             {} primary, {} strict, {} IUPAC-compatible.\n\n\
             Applied detectors, with the catches no other applied detector makes (unique) and \
             the concordant samples only it holds back (sole burden: clearing it would pass \
             them):\n\n\
             | Detector | Samples | Catches primary | Unique primary | Catches strict | Unique strict | Catches IUPAC | Unique IUPAC | Burden | Sole burden | Hits |\n\
             |---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n",
            discordant(Measure::Primary),
            discordant(Measure::Strict),
            discordant(Measure::IupacCompatible),
        );
        for row in rows.iter().filter(|row| row.applied) {
            row.line(text);
        }
        text.push_str(
            "\nMeasured but not applied, with what applying each would add: the not-concordant \
             samples it would newly catch and the concordant samples it would newly hold back \
             (lost A):\n\n\
             | Detector | Samples | Catches primary | Added primary | Catches strict | Added strict | Catches IUPAC | Added IUPAC | Burden | Added burden | Hits |\n\
             |---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n",
        );
        for row in rows.iter().filter(|row| !row.applied) {
            row.line(text);
        }
    }

    fn trace_section(&self, text: &mut String) {
        let mut findings = BTreeMap::<TraceFinding, (usize, BTreeSet<usize>, usize)>::new();
        let (mut poor, mut needed, mut unlabelled) = (0, 0, 0);
        for (index, row, decision, _) in self.scored() {
            let (Some(decision), Some(evidence)) = (decision, evidence(row)) else {
                continue;
            };
            for (trace_index, (trace, found)) in
                evidence.traces.iter().zip(&decision.traces).enumerate()
            {
                if found.is_empty() {
                    continue;
                }
                poor += 1;
                let rerun = decision.rerun.contains(&trace_index);
                needed += usize::from(rerun);
                unlabelled += usize::from(trace.primer.is_none());
                for finding in found {
                    let entry = findings.entry(*finding).or_default();
                    entry.0 += 1;
                    entry.1.insert(index);
                    entry.2 += usize::from(rerun);
                }
            }
        }
        let _ = write!(
            text,
            "\n## Traces\n\n{poor} poor traces; {needed} of them are asked for by a rerun, the \
             others stay diagnostics; {unlabelled} carry no primer in their name.\n\n\
             | Finding | Traces | Samples | In a rerun |\n|---|---:|---:|---:|\n"
        );
        for finding in TraceFinding::ALL {
            let (traces, samples, rerun) = findings
                .get(&finding)
                .map_or((0, 0, 0), |(traces, samples, rerun)| {
                    (*traces, samples.len(), *rerun)
                });
            let _ = writeln!(
                text,
                "| {} | {traces} | {samples} | {rerun} |",
                finding.label()
            );
        }
        text.push_str(
            "\nPrimer footprints learned from the admitted traces of the active samples (proxy \
             for which region a poor trace would have read):\n\n\
             | Primer | Strand | Regions | Admitted traces |\n|---|---|---|---:|\n",
        );
        for (primer, footprint) in self.footprints.iter() {
            let regions = self
                .labels
                .iter()
                .zip(&footprint.regions)
                .filter(|(_, reads)| **reads)
                .map(|(label, _)| label.as_str())
                .collect::<Vec<_>>()
                .join(", ");
            let _ = writeln!(
                text,
                "| {primer} | {} | {regions} | {} |",
                footprint.strand.label(),
                footprint.traces
            );
        }
    }

    fn sweep_section(&self, text: &mut String) {
        text.push_str(
            "\n## Sweep (risk–coverage)\n\nEach policy changes one thing against the draft. \
             Leaks are passed samples not concordant, with the one-sided 95 % upper bound.\n\n\
             | Policy | Pass | Pass share | Pass corpus | Review corpus | Rerun corpus | Leaks primary | Upper | Leaks strict | Upper | Leaks IUPAC | Upper |\n\
             |---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n",
        );
        for (policy, measured) in self.policies.iter().zip(&self.measurements) {
            let corpus = |class: Class| {
                measured
                    .estimate(class)
                    .map_or(0.0, |estimate| estimate.rate)
            };
            let _ = write!(
                text,
                "| {} | {} | {:.4} | {:.4} | {:.4} | {:.4} |",
                policy.name,
                measured.count(Class::Decided(Outcome::Pass)),
                measured.share(Class::Decided(Outcome::Pass)),
                corpus(Class::Decided(Outcome::Pass)),
                corpus(Class::Decided(Outcome::Review)),
                corpus(Class::Decided(Outcome::Rerun)),
            );
            for leak in &measured.leaks {
                let _ = write!(
                    text,
                    " {}/{} | {:.4} |",
                    leak.leaks,
                    leak.passed,
                    upper_bound(leak.leaks, leak.passed)
                );
            }
            text.push('\n');
        }
    }
}

/// `traces` value of `rerun_for` for a sample without an admitted read.
const NO_ADMITTED_READ: &str = "no_admitted_read";

fn strata_section(text: &mut String, measured: &Measurement) {
    let Some(pass) = measured.estimate(Class::Decided(Outcome::Pass)) else {
        return;
    };
    text.push_str(
        "\n### Pass by stratum\n\n| Stratum | Samples | Pass | Review | Rerun |\n\
         |---|---:|---:|---:|---:|\n",
    );
    for stratum in &pass.strata {
        let rate = |class: Class| {
            measured
                .estimate(class)
                .and_then(|estimate| {
                    estimate
                        .strata
                        .iter()
                        .find(|other| other.name == stratum.name)
                })
                .map_or(0.0, |other| other.rate)
        };
        let _ = writeln!(
            text,
            "| {} | {} | {:.3} | {:.3} | {:.3} |",
            stratum.name,
            stratum.samples,
            stratum.rate,
            rate(Class::Decided(Outcome::Review)),
            rate(Class::Decided(Outcome::Rerun)),
        );
    }
}

/// Every discordant position of the scored rows that have region evidence
/// and a comparison, under every measure, with its region's outcome under
/// the reported policy and the sample's class.
fn cases(rows: &[Row<'_>], decisions: &[Option<SampleDecision>], classes: &[Class]) -> Vec<Case> {
    let mut cases = Vec::new();
    for (index, ((row, decision), class)) in rows.iter().zip(decisions).zip(classes).enumerate() {
        if row.sample.retired {
            continue;
        }
        let (Some(decision), Some(evidence)) = (decision, evidence(row)) else {
            continue;
        };
        for measure in Measure::ALL {
            let Ok(comparison) = row.comparison(measure) else {
                continue;
            };
            for position in compare::sorted_positions(&comparison.positions) {
                let anchor = compare::anchor(position);
                let region = anchor.and_then(|anchor| evidence.region_of(anchor));
                let flags = region
                    .and_then(|region| {
                        let evidence = &evidence.regions[region];
                        let anchor = usize::try_from(anchor?).ok()?;
                        evidence
                            .positions
                            .get(anchor.checked_sub(evidence.start)?)
                            .copied()
                    })
                    .map(|at| {
                        Detector::all()
                            .into_iter()
                            .filter(|detector| detector.marks(at))
                            .collect::<Vec<_>>()
                    })
                    .unwrap_or_default();
                cases.push(Case {
                    row: index,
                    sample: *class,
                    measure,
                    region: region
                        .and_then(|region| Some((region, decision.regions.get(region)?.outcome))),
                    position: position.to_owned(),
                    detail: Detail::of(comparison, position),
                    flags,
                });
            }
        }
    }
    cases
}

fn finding_labels(findings: &[TraceFinding]) -> String {
    findings
        .iter()
        .map(|finding| finding.label())
        .collect::<Vec<_>>()
        .join("+")
}

/// A/B/C, leaks, and region outcomes of one policy.
struct Measurement {
    samples: usize,
    classes: BTreeMap<Class, usize>,
    estimates: Vec<(Class, Estimate)>,
    /// Per measure, in [`Measure::ALL`] order.
    leaks: [Leak; 3],
    /// Per region, the outcomes of decided samples with region evidence.
    regions: Vec<Tally>,
}

impl Measurement {
    /// Measures one policy's decisions of `rows` over the scored (not
    /// retired) rows.
    fn of(rows: &[Row<'_>], decisions: &[Option<SampleDecision>], regions: usize) -> Self {
        let mut classes = BTreeMap::<Class, usize>::new();
        let mut scored = Vec::new();
        let mut leaks = [Leak::default(); 3];
        let mut tallies = vec![Tally::default(); regions];
        for (row, decision) in rows.iter().zip(decisions) {
            if row.sample.retired {
                continue;
            }
            let decision = decision.as_ref();
            let class = class(row, decision);
            *classes.entry(class).or_default() += 1;
            let weight = Weighting::Corpus.of(row.sample);
            scored.push((row.sample.stratum.as_str(), weight, class));
            if let (Class::Decided(_), Some(decision)) = (class, decision) {
                for (tally, region) in tallies.iter_mut().zip(&decision.regions) {
                    tally.record(region.outcome, weight);
                }
            }
            if class == Class::Decided(Outcome::Pass) && row.comparisons.is_ok() {
                for measure in Measure::ALL {
                    let leak = &mut leaks[measure.index()];
                    leak.passed += 1;
                    leak.weight += weight;
                    if row.status(measure) != Status::Concordant {
                        leak.leaks += 1;
                        leak.leaked_weight += weight;
                    }
                }
            }
        }
        let estimates = Class::ALL
            .into_iter()
            .filter_map(|wanted| {
                let observations = scored
                    .iter()
                    .map(|(stratum, weight, class)| Observation {
                        stratum,
                        weight: *weight,
                        event: *class == wanted,
                    })
                    .collect::<Vec<_>>();
                estimate::estimate(&observations).map(|estimate| (wanted, estimate))
            })
            .collect();
        Self {
            samples: scored.len(),
            classes,
            estimates,
            leaks,
            regions: tallies,
        }
    }

    fn count(&self, class: Class) -> usize {
        self.classes.get(&class).copied().unwrap_or(0)
    }

    fn share(&self, class: Class) -> f64 {
        if self.samples == 0 {
            0.0
        } else {
            self.count(class) as f64 / self.samples as f64
        }
    }

    fn estimate(&self, class: Class) -> Option<&Estimate> {
        self.estimates
            .iter()
            .find(|(of, _)| *of == class)
            .map(|(_, estimate)| estimate)
    }
}

/// Region outcomes: counts and summed corpus weights.
#[derive(Debug, Clone, Default)]
struct Tally {
    counts: BTreeMap<Outcome, usize>,
    weights: BTreeMap<Outcome, f64>,
}

impl Tally {
    fn record(&mut self, outcome: Outcome, weight: f64) {
        *self.counts.entry(outcome).or_default() += 1;
        *self.weights.entry(outcome).or_default() += weight;
    }

    fn add(&mut self, other: &Self) {
        for (outcome, count) in &other.counts {
            *self.counts.entry(*outcome).or_default() += count;
        }
        for (outcome, weight) in &other.weights {
            *self.weights.entry(*outcome).or_default() += weight;
        }
    }

    fn count(&self, outcome: Outcome) -> usize {
        self.counts.get(&outcome).copied().unwrap_or(0)
    }

    /// The outcome's share of the regions, weighted by `weight_corpus`.
    fn weighted(&self, outcome: Outcome) -> f64 {
        let total = self.weights.values().sum::<f64>();
        if total > 0.0 {
            self.weights.get(&outcome).copied().unwrap_or(0.0) / total
        } else {
            0.0
        }
    }
}

/// Passed samples and those not concordant under one measure.
#[derive(Debug, Clone, Copy, Default)]
struct Leak {
    passed: usize,
    leaks: usize,
    weight: f64,
    leaked_weight: f64,
}

impl Leak {
    fn rate(self) -> f64 {
        if self.passed == 0 {
            0.0
        } else {
            self.leaks as f64 / self.passed as f64
        }
    }

    fn weighted_rate(self) -> f64 {
        if self.weight > 0.0 {
            self.leaked_weight / self.weight
        } else {
            0.0
        }
    }
}

/// One row of the detector table.
#[derive(Debug, Clone, Copy)]
struct DetectorRow {
    detector: Detector,
    applied: bool,
    samples: usize,
    regions: usize,
    positions: usize,
    catches: [usize; 3],
    burden: usize,
    /// Not-concordant samples, per measure, that no other applied detector
    /// flags: unique catches of an applied detector, added catches of one the
    /// policy does not apply.
    alone: [usize; 3],
    /// Concordant samples no other applied detector flags: the sole burden of
    /// an applied detector, the added burden (lost A) of one not applied.
    alone_burden: usize,
    hits: usize,
}

impl DetectorRow {
    const fn new(detector: Detector, applied: bool) -> Self {
        Self {
            detector,
            applied,
            samples: 0,
            regions: 0,
            positions: 0,
            catches: [0; 3],
            burden: 0,
            alone: [0; 3],
            alone_burden: 0,
            hits: 0,
        }
    }

    fn line(&self, text: &mut String) {
        let _ = writeln!(
            text,
            "| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |",
            self.detector.label(),
            self.samples,
            self.catches[0],
            self.alone[0],
            self.catches[1],
            self.alone[1],
            self.catches[2],
            self.alone[2],
            self.burden,
            self.alone_burden,
            self.hits
        );
    }
}

/// One discordant position of one sample under one measure.
struct Case {
    /// The row's index.
    row: usize,
    /// The sample's class under the reported policy.
    sample: Class,
    measure: Measure,
    /// The region's index and outcome; `None` outside every region.
    region: Option<(usize, Outcome)>,
    position: String,
    detail: Detail,
    /// The region detectors that mark the position itself.
    flags: Vec<Detector>,
}

impl Case {
    /// Whether the position lies inside a passed region.
    fn leaks(&self) -> bool {
        self.region
            .is_some_and(|(_, outcome)| outcome == Outcome::Pass)
    }
}

/// What differs at one discordant position.
struct Detail {
    /// `dna_unique`, `truth_unique`, `both`, `uncovered` (a truth-only
    /// allele outside DNA's callable intervals), or `none`.
    kind: &'static str,
    /// `iupac` (a truth IUPAC allele), `insertion`, `deletion`, or
    /// `substitution`.
    class: &'static str,
    dna: Vec<Allele>,
    truth: Vec<Allele>,
}

impl Detail {
    fn of(comparison: &Comparison, position: &str) -> Self {
        let at = |alleles: &[Allele]| {
            alleles
                .iter()
                .filter(|allele| allele.position == position)
                .cloned()
                .collect::<Vec<_>>()
        };
        let dna = at(&comparison.dna_unique);
        let uncovered = at(&comparison.uncovered);
        let truth = at(&comparison.truth_unique);
        let kind = match (dna.is_empty(), truth.is_empty()) {
            (false, false) => "both",
            (false, true) => "dna_unique",
            (true, false) if !uncovered.is_empty() => "uncovered",
            (true, false) => "truth_unique",
            (true, true) => "none",
        };
        let all = dna.iter().chain(&truth).collect::<Vec<_>>();
        let class = if truth.iter().any(Allele::is_iupac) {
            "iupac"
        } else if all.iter().any(|allele| allele.reference == "-") {
            "insertion"
        } else if all.iter().any(|allele| allele.alternate == "-") {
            "deletion"
        } else {
            "substitution"
        };
        Self {
            kind,
            class,
            dna,
            truth,
        }
    }
}

/// The one-sided `1 − ALPHA` Clopper–Pearson upper bound on a rate with
/// `events` of `trials`: the `p` at which observing at most `events` has
/// probability `ALPHA`; `1` when every trial is an event or there is none.
pub(crate) fn upper_bound(events: usize, trials: usize) -> f64 {
    if trials == 0 || events >= trials {
        return 1.0;
    }
    let (mut low, mut high) = (events as f64 / trials as f64, 1.0_f64);
    for _ in 0..200 {
        let middle = f64::midpoint(low, high);
        if binomial_cdf(events, trials, middle) > ALPHA {
            low = middle;
        } else {
            high = middle;
        }
    }
    high
}

/// `P(X ≤ events)` for `X ~ Binomial(trials, p)`, `0 < p < 1`.
fn binomial_cdf(events: usize, trials: usize, p: f64) -> f64 {
    if p <= 0.0 {
        return 1.0;
    }
    if p >= 1.0 {
        return 0.0;
    }
    let (log_p, log_q) = (p.ln(), (-p).ln_1p());
    let n = trials as f64;
    let mut log_choose = 0.0_f64;
    let mut sum = 0.0;
    for i in 0..=events {
        let k = i as f64;
        if i > 0 {
            log_choose += (n - k + 1.0).ln() - k.ln();
        }
        sum += (log_choose + k * log_p + (n - k) * log_q).exp();
    }
    sum.min(1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(left: f64, right: f64) -> bool {
        (left - right).abs() < 1e-9
    }

    #[test]
    fn clopper_pearson_upper_bound() {
        // Zero events: 1 − α^(1/n), the "3/n" rule.
        assert!(close(upper_bound(0, 100), 1.0 - 0.05_f64.powf(0.01)));
        assert!(close(upper_bound(0, 1000), 1.0 - 0.05_f64.powf(0.001)));
        // One event in ten: P(X ≤ 1) = 0.05 at the bound.
        let bound = upper_bound(1, 10);
        assert!(close(binomial_cdf(1, 10, bound), 0.05));
        assert!((bound - 0.3942).abs() < 1e-4);
        assert!(close(upper_bound(3, 3), 1.0));
        assert!(close(upper_bound(0, 0), 1.0));
    }

    #[test]
    fn details_name_the_kind_and_class_of_a_difference() {
        let comparison = Comparison {
            concordance: compare::Concordance::Discordant,
            dna_alleles: 2,
            truth_alleles: 2,
            dna_unique: vec![Allele::new("16189", "T", "C"), Allele::new("310", "-", "C")],
            truth_unique: vec![
                Allele::new("16189", "T", "Y"),
                Allele::new("263", "A", "G"),
                Allele::new("16519", "T", "C"),
            ],
            reasons: BTreeSet::new(),
            positions: BTreeSet::new(),
            iupac_only: false,
            uncovered: vec![Allele::new("16519", "T", "C")],
        };
        let detail = Detail::of(&comparison, "16189");
        assert_eq!((detail.kind, detail.class), ("both", "iupac"));
        let detail = Detail::of(&comparison, "310");
        assert_eq!((detail.kind, detail.class), ("dna_unique", "insertion"));
        let detail = Detail::of(&comparison, "263");
        assert_eq!(
            (detail.kind, detail.class),
            ("truth_unique", "substitution")
        );
        let detail = Detail::of(&comparison, "16519");
        assert_eq!((detail.kind, detail.class), ("uncovered", "substitution"));
        assert_eq!(detail.truth, [Allele::new("16519", "T", "C")]);
        assert_eq!(Detail::of(&comparison, "73").kind, "none");
    }

    /// Four samples through decision and score: a clean pass, a passed
    /// sample discordant only at a truth IUPAC allele (a leak under the
    /// primary and strict measures), a one-strand review that catches a
    /// discordant sample, and a DNA failure.
    #[test]
    fn scores_outcomes_leaks_and_catches() -> Result<(), Box<dyn std::error::Error>> {
        use crate::compare::Failure;
        use crate::evidence::{Position, Region, Strand, Trace};
        use crate::runner::{Called, Identity, Outcome as Run};
        use crate::truth::TruthSample;

        let text = crate::manifest::tests::manifest(&[
            "S0\tb1\tclean\t2\t2\t0\t/t/a.ab1",
            "S1\tb1\tclean\t2\t2\t0\t/t/b.ab1",
            "S2\tb1\thard\t1\t1\t0\t/t/c.ab1",
            "S3\tb1\thard\t1\t1\t0\t/t/d.ab1",
            "S4\tb1\thard\t1\t1\t0\t/t/e.ab1",
        ]);
        let manifest = crate::manifest::parse(text.as_bytes())?;
        let both = Position {
            forward: 1,
            reverse: 1,
            ..Position::default()
        };
        let region = |one_strand: Option<usize>| {
            let mut positions = vec![both; 21];
            if let Some(offset) = one_strand {
                positions[offset].reverse = 0;
            }
            Region {
                start: 16180,
                end: 16200,
                positions,
            }
        };
        let called = |alleles: Vec<Allele>, one_strand: Option<usize>| Run {
            traces: Some(1),
            result: Ok(Called {
                rejected: 0,
                alleles,
                intervals: vec![(16024, 16365)],
                identity: Identity {
                    configuration_sha256: "c".into(),
                    reference_name: "r".into(),
                    reference_sha256: "s".into(),
                    profile_id: "p".into(),
                    profile_sha256: "q".into(),
                },
            }),
            triage: Some(Evidence {
                regions: vec![region(one_strand)],
                traces: vec![Trace {
                    name: "x_HV1F.ab1".into(),
                    primer: Some("HV1F".into()),
                    state: TraceState::Admitted {
                        strand: Strand::Forward,
                        clean: vec![21],
                    },
                    span: Some(400),
                }],
            }),
            millis: 0,
        };
        let c16189 = || Allele::new("16189", "T", "C");
        let outcomes = vec![
            called(vec![c16189()], None),
            called(vec![c16189()], None),
            called(vec![c16189()], Some(10)),
            Run {
                traces: Some(1),
                result: Err(Failure::Dna("unreadable".into())),
                triage: None,
                millis: 0,
            },
            // Every read rejected: DNA reports a failure, triage a rerun.
            Run {
                traces: Some(2),
                result: Err(Failure::Dna("all 2 reads were rejected".into())),
                triage: Some(Evidence {
                    regions: Vec::new(),
                    traces: vec![
                        Trace {
                            name: "y_HV1F.ab1".into(),
                            primer: Some("HV1F".into()),
                            state: TraceState::Weak,
                            span: Some(40),
                        },
                        Trace {
                            name: "y_HV1R.ab1".into(),
                            primer: Some("HV1R".into()),
                            state: TraceState::TooShort,
                            span: None,
                        },
                    ],
                }),
                millis: 0,
            },
        ];
        let truth = |alleles: Vec<Allele>| TruthSample {
            batch: "b1".into(),
            status: Some("OK".into()),
            alleles,
            intervals: vec![(16024, 16365)],
        };
        let records = BTreeMap::from([
            ("S0".to_owned(), truth(vec![c16189()])),
            ("S1".to_owned(), truth(vec![Allele::new("16189", "T", "Y")])),
            (
                "S2".to_owned(),
                truth(vec![c16189(), Allele::new("16190", "T", "C")]),
            ),
            ("S3".to_owned(), truth(vec![c16189()])),
            ("S4".to_owned(), truth(vec![c16189()])),
        ]);
        let rows = crate::report::rows(&manifest.samples, &outcomes, &records);
        let triage = Triage::new(&rows, sweep());
        let files = triage.files().into_iter().collect::<BTreeMap<_, _>>();

        let outcomes = files["triage.tsv"]
            .lines()
            .skip(1)
            .map(|line| line.split('\t').nth(4).unwrap_or_default())
            .collect::<Vec<_>>();
        assert_eq!(outcomes, ["pass", "pass", "review", "failed", "rerun"]);
        assert!(
            files["triage.tsv"]
                .lines()
                .last()
                .is_some_and(|line| line.contains("\ty_HV1F.ab1,y_HV1R.ab1\t"))
        );
        assert!(
            files["triage_traces.tsv"]
                .contains("S4\ty_HV1R.ab1\tHV1R\ttoo_short\t\ttrace.too_short\tno_admitted_read")
        );

        let measured = triage.reported().ok_or("a measurement")?;
        assert_eq!(measured.count(Class::Decided(Outcome::Pass)), 2);
        assert_eq!(measured.count(Class::Decided(Outcome::Rerun)), 1);
        assert_eq!(measured.count(Class::Failed), 1);
        // Region shares count the analysed samples' regions only.
        assert_eq!(measured.regions[0].count(Outcome::Pass), 2);
        assert!(close(
            measured.regions[0].weighted(Outcome::Pass),
            4.0 / 5.0
        ));
        let [primary, strict, compatible] = measured.leaks;
        assert_eq!((primary.passed, primary.leaks), (2, 1));
        assert_eq!((strict.passed, strict.leaks), (2, 1));
        assert_eq!((compatible.passed, compatible.leaks), (2, 0));
        // Corpus-weighted A over S0..S4 (weights 2, 2, 1, 1, 1), by stratum.
        let pass = measured
            .estimate(Class::Decided(Outcome::Pass))
            .ok_or("no estimate")?;
        assert!(close(pass.rate, 4.0 / 7.0));

        let leaks = files["triage_leaks.tsv"]
            .lines()
            .skip(1)
            .collect::<Vec<_>>();
        assert_eq!(
            leaks,
            [
                "S1\tclean\tpass\tprimary\t16180-16200\t16189\thv1_polyc\tboth\tiupac\t16189C\t16189Y\t",
                "S1\tclean\tpass\tstrict\t16180-16200\t16189\thv1_polyc\tboth\tiupac\t16189C\t16189Y\t",
            ]
        );

        let detectors = triage.detector_rows();
        let one_strand = detectors
            .iter()
            .find(|row| row.detector == Detector::OneStrand)
            .ok_or("one-strand row")?;
        assert_eq!(
            (
                one_strand.applied,
                one_strand.samples,
                one_strand.catches,
                one_strand.alone,
                one_strand.burden,
                one_strand.hits
            ),
            (true, 1, [1, 1, 1], [1, 1, 1], 0, 1)
        );
        // Not applied: its counters say what applying it would add, and the
        // TSV leaves its unique columns empty.
        let masked_read = detectors
            .iter()
            .find(|row| row.detector == Detector::Unreliable(Cause::MaskedRead))
            .ok_or("masked-read row")?;
        assert!(!masked_read.applied);
        assert!(files["triage_detectors.tsv"].contains(
            "region.unreliable:masked_read\t0\t0\t0\t0\t0\t0\t0\t0\t-\t-\t-\t-\t0\t0\t0\t0\t0"
        ));
        assert!(files["triage.md"].contains("| primary | 2 | 1 | 0.5000 |"));
        assert!(
            files["triage_sweep.tsv"]
                .lines()
                .nth(1)
                .is_some_and(|line| line.starts_with("draft\t"))
        );
        Ok(())
    }

    #[test]
    fn the_sweep_starts_with_the_draft_and_names_each_policy_once() {
        let policies = sweep();
        assert_eq!(policies.first(), Some(&Policy::draft()));
        let names = policies
            .iter()
            .map(|policy| policy.name.as_str())
            .collect::<BTreeSet<_>>();
        assert_eq!(names.len(), policies.len());
    }
}
