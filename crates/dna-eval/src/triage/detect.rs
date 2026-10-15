//! Draft stage-1 detectors and the outcome of each region, trace, and sample
//! under one policy (PROP-0004, stage R).
//!
//! The detectors are pure functions of a sample's [`Evidence`]:
//!
//! - trace level: `trace.too_short`, `trace.weak`, `trace.unplaced` (the
//!   read's rejection) and `trace.short` (callable span below
//!   [`Policy::short_calls`]); any of them makes the trace *poor*;
//! - region level: `region.unreliable` per [`Cause`], `region.one_strand`
//!   (positions observed cleanly by one strand only), and `region.gap`
//!   (positions observed cleanly by no read).
//!
//! A region is `rerun` when it has a gap, or a one-strand position while the
//! policy applies `region.one_strand`, and a poor trace could have read it
//! ([`Proxy`]); `review` when any finding the policy applies remains; `pass`
//! otherwise. A sample takes its worst region (rerun > review > pass); a
//! sample without an admitted read has no region evidence and is `rerun` with
//! every poor trace.

use std::collections::{BTreeMap, BTreeSet};

use crate::evidence::{Cause, Evidence, Strand, TraceState, primer_strand};

/// Clean positions a read must observe in a region for its primer to count
/// as reading that region.
const FOOTPRINT_POSITIONS: usize = 50;

/// An outcome, ordered from best to worst.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Outcome {
    Pass,
    Review,
    Rerun,
}

impl Outcome {
    pub(crate) const fn label(self) -> &'static str {
        match self {
            Self::Pass => "pass",
            Self::Review => "review",
            Self::Rerun => "rerun",
        }
    }
}

/// A region-level detector.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Detector {
    Unreliable(Cause),
    OneStrand,
    Gap,
}

impl Detector {
    /// Every region detector, in report order.
    pub(crate) fn all() -> Vec<Self> {
        Cause::ALL
            .into_iter()
            .map(Self::Unreliable)
            .chain([Self::OneStrand, Self::Gap])
            .collect()
    }

    pub(crate) fn label(self) -> String {
        match self {
            Self::Unreliable(cause) => format!("region.unreliable:{}", cause.label()),
            Self::OneStrand => "region.one_strand".to_owned(),
            Self::Gap => "region.gap".to_owned(),
        }
    }

    /// Whether the detector marks the position.
    pub(crate) fn marks(self, position: crate::evidence::Position) -> bool {
        match self {
            Self::Unreliable(cause) => position.has(cause),
            Self::OneStrand => position.strands() == 1,
            Self::Gap => position.strands() == 0,
        }
    }
}

/// A trace-level finding.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum TraceFinding {
    TooShort,
    Weak,
    Unplaced,
    /// Rejected for a cause this draft does not know.
    Rejected,
    Short,
}

impl TraceFinding {
    pub(crate) const ALL: [Self; 5] = [
        Self::TooShort,
        Self::Weak,
        Self::Unplaced,
        Self::Rejected,
        Self::Short,
    ];

    pub(crate) const fn label(self) -> &'static str {
        match self {
            Self::TooShort => "trace.too_short",
            Self::Weak => "trace.weak",
            Self::Unplaced => "trace.unplaced",
            Self::Rejected => "trace.rejected",
            Self::Short => "trace.short",
        }
    }
}

/// Which poor traces could have read a region.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Proxy {
    /// A poor trace whose primer reads the region on a strand the region
    /// lacks, by the primer's [`Footprints`].
    Primer,
    /// Any poor trace of the sample.
    Any,
}

impl Proxy {
    pub(crate) const fn label(self) -> &'static str {
        match self {
            Self::Primer => "primer",
            Self::Any => "any",
        }
    }
}

/// The detectors a policy applies and their parameters.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Policy {
    /// Short label of the policy.
    pub(crate) name: String,
    /// `trace.short`: an admitted or rejected trace whose callable span has
    /// fewer calls (`<`) is poor; `0` turns the detector off.
    pub(crate) short_calls: usize,
    /// The causes of `region.unreliable` the policy applies.
    pub(crate) causes: BTreeSet<Cause>,
    /// Whether `region.one_strand` applies.
    pub(crate) one_strand: bool,
    /// Which poor traces could have read a region.
    pub(crate) proxy: Proxy,
}

impl Policy {
    /// The draft default: `trace.short` below 100 calls, every unreliable
    /// cause except `masked_read`, one-strand positions, and the primer proxy.
    pub(crate) fn draft() -> Self {
        Self {
            name: "draft".to_owned(),
            short_calls: 100,
            causes: [
                Cause::N,
                Cause::Contested,
                Cause::Masked,
                Cause::PhaseLoss,
                Cause::ReferenceFrame,
            ]
            .into(),
            one_strand: true,
            proxy: Proxy::Primer,
        }
    }

    /// Whether the policy applies the region detector.
    pub(crate) fn applies(&self, detector: Detector) -> bool {
        match detector {
            Detector::Unreliable(cause) => self.causes.contains(&cause),
            Detector::OneStrand => self.one_strand,
            Detector::Gap => true,
        }
    }

    /// The parameters, as `key=value` pairs separated by spaces.
    pub(crate) fn parameters(&self) -> String {
        format!(
            "short_calls={} causes={} one_strand={} proxy={}",
            self.short_calls,
            if self.causes.is_empty() {
                "-".to_owned()
            } else {
                self.causes
                    .iter()
                    .map(|cause| cause.label())
                    .collect::<Vec<_>>()
                    .join(",")
            },
            self.one_strand,
            self.proxy.label()
        )
    }
}

/// What each labelled primer reads, learned from the admitted traces of the
/// whole run: its strand, and the regions it reads.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct Footprints(BTreeMap<String, Footprint>);

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Footprint {
    pub(crate) strand: Strand,
    /// Per region, whether the primer reads it.
    pub(crate) regions: Vec<bool>,
    /// Admitted traces the footprint is learned from.
    pub(crate) traces: usize,
}

impl Footprints {
    /// A primer's strand is the one its name gives, its regions those that at
    /// least half of its admitted traces observe cleanly at
    /// [`FOOTPRINT_POSITIONS`] positions or more.
    ///
    /// The caller passes the active samples. The footprints are not left out
    /// per sample: a poor trace, the one a footprint is applied to, never
    /// contributes to one (only admitted traces do), and a sample's own
    /// admitted traces are a small share of a primer's (thousands per common
    /// primer on a panel), so leaving them out would cost a recomputation per
    /// sample for no measurable change.
    pub(crate) fn learn<'a>(samples: impl IntoIterator<Item = &'a Evidence>) -> Self {
        let mut tallies = BTreeMap::<String, (usize, Vec<usize>)>::new();
        for evidence in samples {
            for trace in &evidence.traces {
                let (Some(primer), TraceState::Admitted { clean, .. }) =
                    (&trace.primer, &trace.state)
                else {
                    continue;
                };
                let (traces, regions) = tallies
                    .entry(primer.clone())
                    .or_insert_with(|| (0, vec![0; clean.len()]));
                *traces += 1;
                for (count, positions) in regions.iter_mut().zip(clean) {
                    *count += usize::from(*positions >= FOOTPRINT_POSITIONS);
                }
            }
        }
        Self(
            tallies
                .into_iter()
                .filter_map(|(primer, (traces, regions))| {
                    Some((
                        primer.clone(),
                        Footprint {
                            strand: primer_strand(&primer)?,
                            regions: regions.iter().map(|count| 2 * count >= traces).collect(),
                            traces,
                        },
                    ))
                })
                .collect(),
        )
    }

    pub(crate) fn iter(&self) -> impl Iterator<Item = (&String, &Footprint)> {
        self.0.iter()
    }

    fn get(&self, primer: Option<&String>) -> Option<&Footprint> {
        self.0.get(primer?)
    }
}

/// The outcome of one region.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct RegionDecision {
    pub(crate) outcome: Outcome,
    /// Positions each detector marks, every detector included, applied by
    /// the policy or not.
    pub(crate) findings: BTreeMap<Detector, usize>,
    /// The poor traces (indices into the evidence's traces) the rerun asks
    /// for; empty unless the region is `rerun`.
    pub(crate) rerun: Vec<usize>,
}

/// The outcome of one sample.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct SampleDecision {
    pub(crate) outcome: Outcome,
    /// Whether any trace was admitted; without one the sample has no region
    /// evidence and is `rerun` with every poor trace.
    pub(crate) admitted: bool,
    pub(crate) regions: Vec<RegionDecision>,
    /// Per trace, its findings.
    pub(crate) traces: Vec<Vec<TraceFinding>>,
    /// The poor traces (indices into the evidence's traces) the sample's rerun
    /// asks for, in order.
    pub(crate) rerun: Vec<usize>,
}

/// The findings of one trace under `policy`.
fn trace_findings(state: &TraceState, span: Option<usize>, policy: &Policy) -> Vec<TraceFinding> {
    let mut findings = match state {
        TraceState::Admitted { .. } => Vec::new(),
        TraceState::TooShort => vec![TraceFinding::TooShort],
        TraceState::Weak => vec![TraceFinding::Weak],
        TraceState::Unplaced(_) => vec![TraceFinding::Unplaced],
        TraceState::Rejected => vec![TraceFinding::Rejected],
    };
    if span.is_some_and(|span| span < policy.short_calls) {
        findings.push(TraceFinding::Short);
    }
    findings
}

/// Decides one sample under `policy`.
pub(crate) fn decide(
    evidence: &Evidence,
    policy: &Policy,
    footprints: &Footprints,
) -> SampleDecision {
    let traces = evidence
        .traces
        .iter()
        .map(|trace| trace_findings(&trace.state, trace.span, policy))
        .collect::<Vec<_>>();
    let poor = traces
        .iter()
        .enumerate()
        .filter(|(_, findings)| !findings.is_empty())
        .map(|(index, _)| index)
        .collect::<Vec<_>>();
    let regions = evidence
        .regions
        .iter()
        .enumerate()
        .map(|(index, region)| {
            let mut findings = Detector::all()
                .into_iter()
                .map(|detector| (detector, 0))
                .collect::<BTreeMap<_, _>>();
            let mut missing = BTreeSet::new();
            for position in &region.positions {
                for (detector, count) in &mut findings {
                    *count += usize::from(detector.marks(*position));
                }
                match (position.forward > 0, position.reverse > 0) {
                    (false, false) => {
                        missing.insert(Strand::Forward);
                        missing.insert(Strand::Reverse);
                    }
                    (true, false) if policy.one_strand => {
                        missing.insert(Strand::Reverse);
                    }
                    (false, true) if policy.one_strand => {
                        missing.insert(Strand::Forward);
                    }
                    _ => {}
                }
            }
            let rerun = if missing.is_empty() {
                Vec::new()
            } else {
                poor.iter()
                    .copied()
                    .filter(|trace| match policy.proxy {
                        Proxy::Any => true,
                        Proxy::Primer => footprints
                            .get(evidence.traces[*trace].primer.as_ref())
                            .is_some_and(|footprint| {
                                footprint.regions.get(index).copied().unwrap_or(false)
                                    && missing.contains(&footprint.strand)
                            }),
                    })
                    .collect::<Vec<_>>()
            };
            let outcome = if !rerun.is_empty() {
                Outcome::Rerun
            } else if findings
                .iter()
                .any(|(detector, count)| *count > 0 && policy.applies(*detector))
            {
                Outcome::Review
            } else {
                Outcome::Pass
            };
            RegionDecision {
                outcome,
                findings,
                rerun,
            }
        })
        .collect::<Vec<_>>();
    let admitted = evidence
        .traces
        .iter()
        .any(|trace| matches!(trace.state, TraceState::Admitted { .. }));
    let (outcome, rerun) = if admitted {
        let rerun = regions
            .iter()
            .flat_map(|region| region.rerun.iter().copied())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        let worst = regions
            .iter()
            .map(|region| region.outcome)
            .max()
            .unwrap_or(Outcome::Review);
        (worst, rerun)
    } else if poor.is_empty() {
        (Outcome::Review, Vec::new())
    } else {
        (Outcome::Rerun, poor)
    };
    SampleDecision {
        outcome,
        admitted,
        regions,
        traces,
        rerun,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::evidence::{Position, Region, Trace};

    fn position(forward: u8, reverse: u8) -> Position {
        Position {
            forward,
            reverse,
            ..Position::default()
        }
    }

    fn region(positions: Vec<Position>) -> Region {
        Region {
            start: 1,
            end: positions.len(),
            positions,
        }
    }

    fn trace(name: &str, state: TraceState, span: Option<usize>) -> Trace {
        Trace {
            name: name.into(),
            primer: crate::evidence::primer(name),
            state,
            span,
        }
    }

    fn admitted(strand: Strand, clean: usize) -> TraceState {
        TraceState::Admitted {
            strand,
            clean: vec![clean],
        }
    }

    /// A sample with `regions` and `traces`, plus an admitted trace without a
    /// primer after them.
    fn sample(regions: Vec<Region>, mut traces: Vec<Trace>) -> Evidence {
        traces.push(trace("ok.ab1", admitted(Strand::Forward, 0), Some(500)));
        Evidence { regions, traces }
    }

    #[test]
    fn trace_short_applies_below_the_threshold_only() {
        let policy = Policy::draft();
        let state = admitted(Strand::Forward, 100);
        assert_eq!(
            trace_findings(&state, Some(99), &policy),
            [TraceFinding::Short]
        );
        assert!(trace_findings(&state, Some(100), &policy).is_empty());
        assert!(trace_findings(&state, Some(101), &policy).is_empty());
        assert!(trace_findings(&state, None, &policy).is_empty());
        let off = Policy {
            short_calls: 0,
            ..Policy::draft()
        };
        assert!(trace_findings(&state, Some(0), &off).is_empty());
        assert_eq!(
            trace_findings(&TraceState::Weak, Some(10), &policy),
            [TraceFinding::Weak, TraceFinding::Short]
        );
        assert_eq!(
            trace_findings(&TraceState::Unplaced("orientation_tie"), Some(500), &policy),
            [TraceFinding::Unplaced]
        );
        assert_eq!(
            trace_findings(&TraceState::TooShort, None, &policy),
            [TraceFinding::TooShort]
        );
    }

    #[test]
    fn a_region_read_cleanly_on_both_strands_passes() {
        let evidence = sample(vec![region(vec![position(1, 1); 5])], Vec::new());
        let decision = decide(&evidence, &Policy::draft(), &Footprints::default());
        assert_eq!(decision.outcome, Outcome::Pass);
        assert!(
            decision.regions[0]
                .findings
                .values()
                .all(|count| *count == 0)
        );
    }

    #[test]
    fn one_strand_reviews_or_passes_by_policy() {
        let evidence = sample(
            vec![region(vec![position(1, 1), position(2, 0)])],
            Vec::new(),
        );
        let decision = decide(&evidence, &Policy::draft(), &Footprints::default());
        assert_eq!(decision.outcome, Outcome::Review);
        assert_eq!(decision.regions[0].findings[&Detector::OneStrand], 1);
        let lenient = Policy {
            one_strand: false,
            ..Policy::draft()
        };
        assert_eq!(
            decide(&evidence, &lenient, &Footprints::default()).outcome,
            Outcome::Pass
        );
    }

    #[test]
    fn an_unreliable_cause_reviews_only_when_applied() {
        let mut masked = position(1, 1);
        masked.masked_reads = 1;
        let evidence = sample(vec![region(vec![masked])], Vec::new());
        let draft = Policy::draft();
        let decision = decide(&evidence, &draft, &Footprints::default());
        assert_eq!(decision.outcome, Outcome::Pass);
        assert_eq!(
            decision.regions[0].findings[&Detector::Unreliable(Cause::MaskedRead)],
            1
        );
        let mut causes = draft.causes.clone();
        causes.insert(Cause::MaskedRead);
        let strict = Policy { causes, ..draft };
        assert_eq!(
            decide(&evidence, &strict, &Footprints::default()).outcome,
            Outcome::Review
        );
    }

    #[test]
    fn a_gap_reruns_when_a_poor_trace_could_read_it() {
        // Learned footprints: HV1R reads region 0 on the reverse strand.
        let learned = Evidence {
            regions: Vec::new(),
            traces: vec![
                trace("a_HV1R.ab1", admitted(Strand::Reverse, 50), Some(500)),
                trace("b_HV1R.ab1", admitted(Strand::Reverse, 49), Some(500)),
                trace("c_HV2F.ab1", admitted(Strand::Forward, 0), Some(500)),
            ],
        };
        let footprints = Footprints::learn([&learned]);
        let hv1r = footprints.iter().find(|(primer, _)| *primer == "HV1R");
        assert_eq!(
            hv1r.map(|(_, footprint)| (footprint.strand, footprint.regions.clone())),
            Some((Strand::Reverse, vec![true]))
        );

        // Forward-only region, with a weak HV1R and a weak HV2F trace.
        let evidence = sample(
            vec![region(vec![position(1, 0), position(1, 1)])],
            vec![
                trace("s_HV1R.ab1", TraceState::Weak, Some(50)),
                trace("s_HV2F.ab1", TraceState::Weak, Some(50)),
            ],
        );
        let decision = decide(&evidence, &Policy::draft(), &footprints);
        assert_eq!(decision.outcome, Outcome::Rerun);
        assert_eq!(decision.regions[0].rerun, [0]);
        let any = Policy {
            proxy: Proxy::Any,
            ..Policy::draft()
        };
        assert_eq!(
            decide(&evidence, &any, &footprints).regions[0].rerun,
            [0, 1]
        );
        // Without one-strand, the region passes and the poor traces stay
        // diagnostics.
        let lenient = Policy {
            one_strand: false,
            ..Policy::draft()
        };
        let decision = decide(&evidence, &lenient, &footprints);
        assert_eq!(decision.outcome, Outcome::Pass);
        assert_eq!(
            decision.traces[0],
            [TraceFinding::Weak, TraceFinding::Short]
        );

        // A gap needs any strand: still rerun, even without one-strand.
        let gap = Evidence {
            regions: vec![region(vec![position(0, 0)])],
            traces: evidence.traces.clone(),
        };
        assert_eq!(decide(&gap, &lenient, &footprints).outcome, Outcome::Rerun);
        // A forward poor trace cannot supply a missing reverse strand.
        let forward_only = sample(
            vec![region(vec![position(1, 0)])],
            vec![trace("s_HV2F.ab1", TraceState::Weak, None)],
        );
        let mut footprints_hv2f = Footprints::learn([&Evidence {
            regions: Vec::new(),
            traces: vec![trace("c_HV2F.ab1", admitted(Strand::Forward, 80), None)],
        }]);
        assert_eq!(
            decide(&forward_only, &Policy::draft(), &footprints_hv2f).outcome,
            Outcome::Review
        );
        footprints_hv2f = Footprints::default();
        assert_eq!(
            decide(&forward_only, &Policy::draft(), &footprints_hv2f).outcome,
            Outcome::Review
        );
    }

    #[test]
    fn the_sample_takes_its_worst_region() {
        let evidence = sample(
            vec![
                region(vec![position(1, 1)]),
                region(vec![position(1, 0)]),
                region(vec![position(0, 0)]),
            ],
            vec![trace("s_HV1R.ab1", TraceState::TooShort, None)],
        );
        let decision = decide(
            &evidence,
            &Policy {
                proxy: Proxy::Any,
                ..Policy::draft()
            },
            &Footprints::default(),
        );
        assert_eq!(
            decision
                .regions
                .iter()
                .map(|region| region.outcome)
                .collect::<Vec<_>>(),
            [Outcome::Pass, Outcome::Rerun, Outcome::Rerun]
        );
        assert_eq!(decision.outcome, Outcome::Rerun);
        assert_eq!(decision.rerun, [0]);
        let review = decide(
            &sample(evidence.regions.clone(), Vec::new()),
            &Policy::draft(),
            &Footprints::default(),
        );
        assert_eq!(review.outcome, Outcome::Review);
    }

    #[test]
    fn a_sample_without_an_admitted_read_reruns_every_poor_trace() {
        let evidence = Evidence {
            regions: Vec::new(),
            traces: vec![
                trace("s.ab1", TraceState::Weak, Some(40)),
                trace("s_HV1R.ab1", TraceState::TooShort, None),
            ],
        };
        let decision = decide(&evidence, &Policy::draft(), &Footprints::default());
        assert!(!decision.admitted);
        assert_eq!(
            (decision.outcome, decision.rerun),
            (Outcome::Rerun, vec![0, 1])
        );
    }
}
