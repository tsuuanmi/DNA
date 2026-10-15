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

/// The triage of every panel sample under every policy of the sweep.
pub(crate) struct Triage<'a> {
    rows: &'a [Row<'a>],
    footprints: Footprints,
    policies: Vec<Policy>,
    /// Per policy, per row: the decision, or `None` when DNA failed.
    decisions: Vec<Vec<Option<SampleDecision>>>,
}

/// The triage evidence of a row, when DNA analysed the sample.
fn evidence<'a>(row: &'a Row<'a>) -> Option<&'a Evidence> {
    row.outcome
        .result
        .as_ref()
        .ok()
        .and_then(|called| called.triage.as_ref())
}

impl<'a> Triage<'a> {
    /// Decides every row under every policy; the first policy is the one the
    /// per-sample outputs report.
    pub(crate) fn new(rows: &'a [Row<'a>], policies: Vec<Policy>) -> Self {
        let footprints = Footprints::learn(rows.iter().filter_map(evidence));
        let decisions = policies
            .iter()
            .map(|policy| {
                rows.iter()
                    .map(|row| evidence(row).map(|evidence| decide(evidence, policy, &footprints)))
                    .collect()
            })
            .collect();
        Self {
            rows,
            footprints,
            policies,
            decisions,
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

    /// The scored rows with their decision under policy `index`.
    fn scored(&self, index: usize) -> impl Iterator<Item = (&'a Row<'a>, Option<&SampleDecision>)> {
        self.rows
            .iter()
            .zip(&self.decisions[index])
            .filter(|(row, _)| !row.sample.retired)
            .map(|(row, decision)| (row, decision.as_ref()))
    }

    /// Region labels, from the first sample with evidence.
    fn region_labels(&self) -> Vec<String> {
        self.rows
            .iter()
            .find_map(evidence)
            .map(|evidence| {
                evidence
                    .regions
                    .iter()
                    .map(crate::evidence::Region::label)
                    .collect()
            })
            .unwrap_or_default()
    }

    /// A/B/C and leaks of policy `index`.
    fn measure(&self, index: usize) -> Measurement {
        let mut classes = BTreeMap::<Class, usize>::new();
        let mut observations = Vec::new();
        let mut leaks = [Leak::default(); 3];
        for (row, decision) in self.scored(index) {
            let class = decision.map_or(Class::Failed, |decision| Class::Decided(decision.outcome));
            *classes.entry(class).or_default() += 1;
            observations.push((row, class));
            if class == Class::Decided(Outcome::Pass) && row.comparisons.is_ok() {
                let weight = Weighting::Corpus.of(row.sample);
                for measure in Measure::ALL {
                    let leak = &mut leaks[measure.index()];
                    let leaked = row.status(measure) != Status::Concordant;
                    leak.passed += 1;
                    leak.weight += weight;
                    if leaked {
                        leak.leaks += 1;
                        leak.leaked_weight += weight;
                    }
                }
            }
        }
        let estimates = Class::ALL
            .into_iter()
            .filter_map(|class| {
                let rows = observations
                    .iter()
                    .map(|(row, of)| Observation {
                        stratum: &row.sample.stratum,
                        weight: Weighting::Corpus.of(row.sample),
                        event: *of == class,
                    })
                    .collect::<Vec<_>>();
                estimate::estimate(&rows).map(|estimate| (class, estimate))
            })
            .collect();
        Measurement {
            samples: observations.len(),
            classes,
            estimates,
            leaks,
        }
    }

    /// `triage.tsv`: one row per sample under the reported policy.
    fn samples(&self) -> String {
        let labels = self.region_labels();
        let mut text = String::from(
            "sample_id\tstratum\tweight_corpus\tretired\toutcome\tregions\tfindings\t\
             poor_traces\trerun_traces\tstatus\tstrict_status\tiupac_compatible_status\n",
        );
        for (row, decision) in self.rows.iter().zip(&self.decisions[0]) {
            let (outcome, regions, findings, poor, rerun) = decision.as_ref().map_or_else(
                || {
                    (
                        "failed".to_owned(),
                        String::new(),
                        String::new(),
                        String::new(),
                        String::new(),
                    )
                },
                |decision| {
                    let traces = evidence(row).map_or(&[][..], |evidence| &evidence.traces[..]);
                    let regions = labels
                        .iter()
                        .zip(&decision.regions)
                        .map(|(label, region)| format!("{label}={}", region.outcome.label()))
                        .collect::<Vec<_>>()
                        .join(",");
                    let findings = labels
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
                        .regions
                        .iter()
                        .flat_map(|region| &region.rerun)
                        .copied()
                        .collect::<BTreeSet<_>>()
                        .into_iter()
                        .filter_map(|index| traces.get(index).map(|trace| trace.name.clone()))
                        .collect::<Vec<_>>()
                        .join(",");
                    (
                        decision.outcome.label().to_owned(),
                        regions,
                        findings,
                        poor,
                        rerun,
                    )
                },
            );
            let _ = writeln!(
                text,
                "{}\t{}\t{}\t{}\t{outcome}\t{regions}\t{findings}\t{poor}\t{rerun}\t{}\t{}\t{}",
                row.sample.id,
                row.sample.stratum,
                row.sample.weight_corpus.text,
                u8::from(row.sample.retired),
                row.status(Measure::Primary).label(),
                row.status(Measure::Strict).label(),
                row.status(Measure::IupacCompatible).label(),
            );
        }
        text
    }

    /// `triage_regions.tsv`: one row per sample and region.
    fn regions(&self) -> String {
        let labels = self.region_labels();
        let detectors = Detector::all();
        let mut text = String::from("sample_id\tstratum\tregion\toutcome");
        for detector in &detectors {
            let _ = write!(text, "\t{}", detector.label());
        }
        text.push_str(
            "\trerun_traces\tprimary_positions\tstrict_positions\tiupac_compatible_positions\n",
        );
        for (row, decision) in self.rows.iter().zip(&self.decisions[0]) {
            let (Some(decision), Some(evidence)) = (decision, evidence(row)) else {
                continue;
            };
            for (index, (label, region)) in labels.iter().zip(&decision.regions).enumerate() {
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
        let labels = self.region_labels();
        let mut text =
            String::from("sample_id\ttrace\tprimer\tstate\tcallable_span\tfindings\trerun_for\n");
        for (row, decision) in self.rows.iter().zip(&self.decisions[0]) {
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
                let rerun_for = labels
                    .iter()
                    .zip(&decision.regions)
                    .filter(|(_, region)| region.rerun.contains(&index))
                    .map(|(label, _)| label.as_str())
                    .collect::<Vec<_>>()
                    .join(",");
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

    /// Every discordant position of every scored, compared sample under every
    /// measure, with its region's outcome under the reported policy.
    fn cases(&self) -> Vec<Case<'a>> {
        let mut cases = Vec::new();
        for (row, decision) in self.scored(0) {
            let (Some(decision), Some(evidence)) = (decision, evidence(row)) else {
                continue;
            };
            for measure in Measure::ALL {
                let Ok(comparison) = row.comparison(measure) else {
                    continue;
                };
                for position in compare::sorted_positions(&comparison.positions) {
                    let region =
                        compare::anchor(position).and_then(|anchor| evidence.region_of(anchor));
                    let flags = region
                        .and_then(|region| {
                            let evidence = &evidence.regions[region];
                            let anchor = usize::try_from(compare::anchor(position)?).ok()?;
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
                        row,
                        measure,
                        region: region.map(|index| {
                            (
                                evidence.regions[index].label(),
                                decision.regions[index].outcome,
                            )
                        }),
                        position: position.to_owned(),
                        detail: Detail::of(comparison, position),
                        flags,
                    });
                }
            }
        }
        cases
    }

    /// `triage_leaks.tsv`: discordant positions inside passed regions.
    fn leak_cases(&self) -> String {
        let mut text = String::from(
            "sample_id\tstratum\tsample_outcome\tmeasure\tregion\tposition\twindow\tkind\t\
             class\tdna_alleles\ttruth_alleles\tmarked_by\n",
        );
        for case in self.cases().iter().filter(|case| case.leaks()) {
            let outcome = self.decisions[0]
                .iter()
                .zip(self.rows)
                .find(|(_, row)| row.sample.id == case.row.sample.id)
                .and_then(|(decision, _)| decision.as_ref())
                .map_or("failed", |decision| decision.outcome.label());
            let _ = writeln!(
                text,
                "{}\t{}\t{outcome}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                case.row.sample.id,
                case.row.sample.stratum,
                case.measure.label(),
                case.region
                    .as_ref()
                    .map_or("-", |(label, _)| label.as_str()),
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

    /// Per region detector: samples flagged, catches and unique catches per
    /// measure, burden, and positional hits, under the reported policy's
    /// findings.
    fn detector_rows(&self) -> Vec<DetectorRow> {
        let Some(policy) = self.policy() else {
            return Vec::new();
        };
        let compared = self
            .scored(0)
            .filter_map(|(row, decision)| Some((row, decision?)))
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
        let cases = self.cases();
        Detector::all()
            .into_iter()
            .map(|detector| {
                let others = Detector::all()
                    .into_iter()
                    .filter(|other| *other != detector && policy.applies(*other))
                    .collect::<Vec<_>>();
                let mut row = DetectorRow {
                    detector,
                    applied: policy.applies(detector),
                    ..DetectorRow::default()
                };
                for (sample, decision) in &compared {
                    if !flagged(decision, detector) {
                        continue;
                    }
                    row.samples += 1;
                    row.regions += decision
                        .regions
                        .iter()
                        .filter(|region| {
                            region
                                .findings
                                .get(&detector)
                                .is_some_and(|count| *count > 0)
                        })
                        .count();
                    row.positions += decision
                        .regions
                        .iter()
                        .map(|region| region.findings.get(&detector).copied().unwrap_or(0))
                        .sum::<usize>();
                    let alone = !others.iter().any(|other| flagged(decision, *other));
                    for measure in Measure::ALL {
                        if sample.status(measure) != Status::Concordant {
                            row.catches[measure.index()] += 1;
                            row.unique[measure.index()] += usize::from(alone);
                        }
                    }
                    let concordant = sample.status(Measure::Primary) == Status::Concordant;
                    row.burden += usize::from(concordant);
                    row.sole_burden += usize::from(concordant && alone);
                }
                row.hits = cases
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
            .scored(0)
            .filter_map(|(row, _)| evidence(row))
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

    /// `triage_detectors.tsv`.
    fn detectors(&self) -> String {
        let mut text = String::from(
            "detector\tapplied\tsamples\tregions\tpositions\tcatches_primary\tunique_primary\t\
             catches_strict\tunique_strict\tcatches_iupac_compatible\tunique_iupac_compatible\t\
             burden\tsole_burden\tposition_hits_primary\n",
        );
        for row in self.detector_rows() {
            let _ = writeln!(
                text,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                row.detector.label(),
                u8::from(row.applied),
                row.samples,
                row.regions,
                row.positions,
                row.catches[0],
                row.unique[0],
                row.catches[1],
                row.unique[1],
                row.catches[2],
                row.unique[2],
                row.burden,
                row.sole_burden,
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
        for (index, policy) in self.policies.iter().enumerate() {
            let measured = self.measure(index);
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
        let measured = self.measure(0);
        text.push_str(
            "\n## Samples\n\n| Outcome | Samples | Share | Corpus-weighted | 95 % CI |\n\
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
        Self::strata_section(&mut text, &measured);
        self.leak_section(&mut text, &measured);
        self.region_section(&mut text);
        self.class_section(&mut text);
        self.detector_section(&mut text);
        self.trace_section(&mut text);
        self.sweep_section(&mut text);
        text
    }

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
        let cases = self.cases();
        let mut composition = BTreeMap::<(Measure, &str), usize>::new();
        for (row, decision) in self.scored(0) {
            if decision.is_none_or(|decision| decision.outcome != Outcome::Pass)
                || row.comparisons.is_err()
            {
                continue;
            }
            for measure in Measure::ALL {
                if row.status(measure) == Status::Concordant {
                    continue;
                }
                let at = cases
                    .iter()
                    .filter(|case| case.measure == measure && case.row.sample.id == row.sample.id)
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

    fn region_section(&self, text: &mut String) {
        let labels = self.region_labels();
        let mut counts = vec![BTreeMap::<Outcome, usize>::new(); labels.len()];
        for (_, decision) in self.scored(0) {
            for (count, region) in counts
                .iter_mut()
                .zip(decision.map_or(&[][..], |decision| &decision.regions[..]))
            {
                *count.entry(region.outcome).or_default() += 1;
            }
        }
        let cases = self.cases();
        text.push_str(
            "\n## Regions\n\nOutcomes of the analysed samples' regions, and discordant positions \
             inside passed regions (region leaks) per measure.\n\n\
             | Region | Pass | Review | Rerun | Pass share | Leaks primary | Leaks strict | Leaks IUPAC-compatible |\n\
             |---|---:|---:|---:|---:|---:|---:|---:|\n",
        );
        let mut total = BTreeMap::<Outcome, usize>::new();
        for (label, count) in labels.iter().zip(&counts) {
            for (outcome, value) in count {
                *total.entry(*outcome).or_default() += value;
            }
            region_line(text, label, count, &cases, Some(label));
        }
        region_line(text, "all", &total, &cases, None);
    }

    fn class_section(&self, text: &mut String) {
        let mut classes =
            BTreeMap::<(Measure, String, &str, &str, &str), (usize, BTreeSet<&str>)>::new();
        for case in self.cases().iter().filter(|case| case.leaks()) {
            let entry = classes
                .entry((
                    case.measure,
                    case.region
                        .as_ref()
                        .map_or_else(|| "-".to_owned(), |(label, _)| label.clone()),
                    Window::of(&case.position).label(),
                    case.detail.kind,
                    case.detail.class,
                ))
                .or_default();
            entry.0 += 1;
            entry.1.insert(case.row.sample.id.as_str());
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
            self.scored(0)
                .filter(|(row, decision)| decision.is_some() && row.comparisons.is_ok())
                .filter(|(row, _)| row.status(measure) != Status::Concordant)
                .count()
        };
        let _ = write!(
            text,
            "\n## Detectors\n\nPer region detector, over analysed and compared samples: samples \
             it flags, not-concordant samples among them (catches) and those no other applied \
             detector flags (unique), concordant samples it flags (burden) and those no other \
             applied detector flags (sole burden: the samples clearing it would pass), and \
             primary discordant positions it marks itself (hits). Not concordant: {} primary, \
             {} strict, {} IUPAC-compatible.\n\n\
             | Detector | Applied | Samples | Catches primary | Unique primary | Catches strict | Unique strict | Catches IUPAC | Unique IUPAC | Burden | Sole burden | Hits |\n\
             |---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n",
            discordant(Measure::Primary),
            discordant(Measure::Strict),
            discordant(Measure::IupacCompatible),
        );
        for row in self.detector_rows() {
            let _ = writeln!(
                text,
                "| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |",
                row.detector.label(),
                if row.applied { "yes" } else { "no" },
                row.samples,
                row.catches[0],
                row.unique[0],
                row.catches[1],
                row.unique[1],
                row.catches[2],
                row.unique[2],
                row.burden,
                row.sole_burden,
                row.hits
            );
        }
    }

    fn trace_section(&self, text: &mut String) {
        let mut findings = BTreeMap::<TraceFinding, (usize, BTreeSet<&str>, usize)>::new();
        let (mut poor, mut needed, mut unlabelled) = (0, 0, 0);
        for (row, decision) in self.scored(0) {
            let (Some(decision), Some(evidence)) = (decision, evidence(row)) else {
                continue;
            };
            for (index, (trace, found)) in evidence.traces.iter().zip(&decision.traces).enumerate()
            {
                if found.is_empty() {
                    continue;
                }
                poor += 1;
                let rerun = decision
                    .regions
                    .iter()
                    .any(|region| region.rerun.contains(&index));
                needed += usize::from(rerun);
                unlabelled += usize::from(trace.primer.is_none());
                for finding in found {
                    let entry = findings.entry(*finding).or_default();
                    entry.0 += 1;
                    entry.1.insert(row.sample.id.as_str());
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
        let labels = self.region_labels();
        text.push_str(
            "\nPrimer footprints learned from the admitted traces (proxy for which region a \
             poor trace would have read):\n\n| Primer | Strand | Regions | Admitted traces |\n\
             |---|---|---|---:|\n",
        );
        for (primer, footprint) in self.footprints.iter() {
            let regions = labels
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
        for (index, policy) in self.policies.iter().enumerate() {
            let measured = self.measure(index);
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

fn region_line(
    text: &mut String,
    label: &str,
    counts: &BTreeMap<Outcome, usize>,
    cases: &[Case<'_>],
    region: Option<&String>,
) {
    let count = |outcome: Outcome| counts.get(&outcome).copied().unwrap_or(0);
    let total = Outcome::ALL.into_iter().map(count).sum::<usize>();
    let leaks = |measure: Measure| {
        cases
            .iter()
            .filter(|case| case.measure == measure && case.leaks())
            .filter(|case| {
                region.is_none_or(|region| {
                    case.region
                        .as_ref()
                        .is_some_and(|(label, _)| label == region)
                })
            })
            .count()
    };
    let _ = writeln!(
        text,
        "| {label} | {} | {} | {} | {:.4} | {} | {} | {} |",
        count(Outcome::Pass),
        count(Outcome::Review),
        count(Outcome::Rerun),
        if total == 0 {
            0.0
        } else {
            count(Outcome::Pass) as f64 / total as f64
        },
        leaks(Measure::Primary),
        leaks(Measure::Strict),
        leaks(Measure::IupacCompatible),
    );
}

fn finding_labels(findings: &[TraceFinding]) -> String {
    findings
        .iter()
        .map(|finding| finding.label())
        .collect::<Vec<_>>()
        .join("+")
}

/// A/B/C and leaks of one policy.
struct Measurement {
    samples: usize,
    classes: BTreeMap<Class, usize>,
    estimates: Vec<(Class, Estimate)>,
    /// Per measure, in [`Measure::ALL`] order.
    leaks: [Leak; 3],
}

impl Measurement {
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
    unique: [usize; 3],
    burden: usize,
    /// Concordant samples no other applied detector flags: the samples that
    /// would pass if this detector were cleared.
    sole_burden: usize,
    hits: usize,
}

impl Default for DetectorRow {
    fn default() -> Self {
        Self {
            detector: Detector::Gap,
            applied: false,
            samples: 0,
            regions: 0,
            positions: 0,
            catches: [0; 3],
            unique: [0; 3],
            burden: 0,
            sole_burden: 0,
            hits: 0,
        }
    }
}

/// One discordant position of one sample under one measure.
struct Case<'a> {
    row: &'a Row<'a>,
    measure: Measure,
    /// The region's label and outcome; `None` outside every region.
    region: Option<(String, Outcome)>,
    position: String,
    detail: Detail,
    /// The region detectors that mark the position itself.
    flags: Vec<Detector>,
}

impl Case<'_> {
    /// Whether the position lies inside a passed region.
    fn leaks(&self) -> bool {
        self.region
            .as_ref()
            .is_some_and(|(_, outcome)| *outcome == Outcome::Pass)
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
        ]);
        let rows = crate::report::rows(&manifest.samples, &outcomes, &records);
        let triage = Triage::new(&rows, sweep());
        let files = triage.files().into_iter().collect::<BTreeMap<_, _>>();

        let outcomes = files["triage.tsv"]
            .lines()
            .skip(1)
            .map(|line| line.split('\t').nth(4).unwrap_or_default())
            .collect::<Vec<_>>();
        assert_eq!(outcomes, ["pass", "pass", "review", "failed"]);

        let measured = triage.measure(0);
        assert_eq!(measured.count(Class::Decided(Outcome::Pass)), 2);
        assert_eq!(measured.count(Class::Failed), 1);
        let [primary, strict, compatible] = measured.leaks;
        assert_eq!((primary.passed, primary.leaks), (2, 1));
        assert_eq!((strict.passed, strict.leaks), (2, 1));
        assert_eq!((compatible.passed, compatible.leaks), (2, 0));
        // Corpus-weighted A over S0..S3 (weights 2, 2, 1, 1), by stratum.
        let pass = measured
            .estimate(Class::Decided(Outcome::Pass))
            .ok_or("no estimate")?;
        assert!(close(pass.rate, 4.0 / 6.0));

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
                one_strand.samples,
                one_strand.catches,
                one_strand.unique,
                one_strand.burden,
                one_strand.hits
            ),
            (1, [1, 1, 1], [1, 1, 1], 0, 1)
        );
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
