//! What each `region.minor_signal` rule of the sweep buys and costs: per
//! reported policy with a rule, against the same policy without it (decided
//! by [`decide`](crate::triage::detect::decide) like any other policy), the
//! samples it moves from pass to review, the leaks it removes, its burden,
//! the corpus-weighted A it costs, and its recall of truth IUPAC sites.

use std::fmt::Write as _;

use super::{Class, Detail, Leak, Triage, class, evidence, upper_bound};
use crate::compare::{self, Measure};
use crate::estimate::{self, Observation};
use crate::report::{Row, Status, Weighting};
use crate::triage::detect::{Detector, Outcome, Policy, SampleDecision};
use crate::triage::minor::Rule;

/// One minor-signal rule's effect under one policy, over the scored rows.
#[derive(Debug, Clone, Copy, Default)]
pub(super) struct Effect {
    /// Corpus-weighted A without the rule, and with it.
    pass_without: f64,
    pass_with: f64,
    /// Passed samples and leaks without the rule and with it, per measure in
    /// [`Measure::ALL`] order.
    leaks_without: [Leak; Measure::COUNT],
    leaks_with: [Leak; Measure::COUNT],
    /// Compared samples the rule alone sends from pass to review, and the
    /// concordant ones among them (burden).
    moved: usize,
    burden: usize,
    /// Primary truth IUPAC positions inside a region that passes without
    /// the rule, and those the rule marks itself (recall).
    sites: usize,
    caught: usize,
    /// Primary truth IUPAC positions inside any region, and those marked.
    sites_all: usize,
    caught_all: usize,
}

/// Whether the minor-signal measurement reports policy `index`: a policy
/// with a rule that is the reported (first) policy or has a rule other than
/// the reported policy's.
pub(super) fn reports(index: usize, policy: &Policy, reported: Option<Rule>) -> bool {
    policy.minor.is_some() && (index == 0 || policy.minor != reported)
}

/// Measures `policy`'s minor-signal rule from its decisions of `rows` and
/// the decisions of the same policy without the rule; `None` when the
/// policy applies none.
pub(super) fn effect(
    rows: &[Row<'_>],
    decisions: &[Option<SampleDecision>],
    without: &[Option<SampleDecision>],
    policy: &Policy,
) -> Option<Effect> {
    let rule = Detector::MinorSignal(policy.minor?);
    let mut effect = Effect::default();
    let mut shares_without = Vec::new();
    let mut shares_with = Vec::new();
    for ((row, decision), cleared) in rows.iter().zip(decisions).zip(without) {
        if row.sample.retired {
            continue;
        }
        let (decision, cleared) = (decision.as_ref(), cleared.as_ref());
        let passed = class(row, decision) == Class::Decided(Outcome::Pass);
        let free = class(row, cleared) == Class::Decided(Outcome::Pass);
        let weight = Weighting::Corpus.of(row.sample);
        let stratum = row.sample.stratum.as_str();
        shares_without.push(Observation {
            stratum,
            weight,
            event: free,
        });
        shares_with.push(Observation {
            stratum,
            weight,
            event: passed,
        });
        let (Some(cleared), Some(evidence), Ok(primary)) =
            (cleared, evidence(row), row.comparison(Measure::Primary))
        else {
            continue;
        };
        for measure in Measure::ALL {
            let leaked = row.status(measure) != Status::Concordant;
            for (leaks, counted) in [
                (&mut effect.leaks_without, free),
                (&mut effect.leaks_with, passed),
            ] {
                if counted {
                    let leak = &mut leaks[measure.index()];
                    leak.passed += 1;
                    leak.weight += weight;
                    if leaked {
                        leak.leaks += 1;
                        leak.leaked_weight += weight;
                    }
                }
            }
        }
        if free && !passed {
            effect.moved += 1;
            effect.burden += usize::from(row.status(Measure::Primary) == Status::Concordant);
        }
        for position in compare::sorted_positions(&primary.positions) {
            if Detail::of(primary, position).class != "iupac" {
                continue;
            }
            let Some((region, offset)) =
                compare::anchor(position).and_then(|anchor| evidence.locate(anchor))
            else {
                continue;
            };
            let (Some(at), Some(outcome)) = (
                evidence.regions[region].positions.get(offset).copied(),
                cleared.regions.get(region),
            ) else {
                continue;
            };
            let marked = rule.marks(at, evidence.minors(region, offset));
            effect.sites_all += 1;
            effect.caught_all += usize::from(marked);
            if outcome.outcome == Outcome::Pass {
                effect.sites += 1;
                effect.caught += usize::from(marked);
            }
        }
    }
    let rate = |observations: &[Observation<'_>]| {
        estimate::estimate(observations).map_or(0.0, |estimate| estimate.rate)
    };
    effect.pass_without = rate(&shares_without);
    effect.pass_with = rate(&shares_with);
    Some(effect)
}

impl Triage<'_> {
    /// The policies the minor-signal measurement reports, with the rule and
    /// its effect.
    fn minor_rows(&self) -> Vec<(&Policy, Rule, Effect)> {
        self.policies
            .iter()
            .zip(&self.minor)
            .filter_map(|(policy, effect)| Some((policy, policy.minor?, (*effect)?)))
            .collect()
    }

    /// `triage_minor.tsv`: the reported policy and each policy with another
    /// minor-signal rule.
    pub(super) fn minor_table(&self) -> String {
        let mut text = String::from(
            "policy\trule\tpass_corpus_without\tpass_corpus\tpass_corpus_cost\tmoved\tburden\t\
             iupac_sites\tiupac_caught\tiupac_sites_all\tiupac_caught_all",
        );
        for measure in Measure::ALL {
            let label = measure.label();
            let _ = write!(
                text,
                "\tpassed_without_{label}\tleaks_without_{label}\tpassed_{label}\tleaks_{label}\t\
                 leak_upper95_{label}\tleak_rate_corpus_{label}"
            );
        }
        text.push('\n');
        for (policy, rule, effect) in self.minor_rows() {
            let _ = write!(
                text,
                "{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}",
                policy.name,
                rule.label(),
                effect.pass_without,
                effect.pass_with,
                effect.pass_without - effect.pass_with,
                effect.moved,
                effect.burden,
                effect.sites,
                effect.caught,
                effect.sites_all,
                effect.caught_all,
            );
            for (without, with) in effect.leaks_without.iter().zip(&effect.leaks_with) {
                let _ = write!(
                    text,
                    "\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}",
                    without.passed,
                    without.leaks,
                    with.passed,
                    with.leaks,
                    upper_bound(with.leaks, with.passed),
                    with.weighted_rate()
                );
            }
            text.push('\n');
        }
        text
    }

    /// The minor-signal section of `triage.md`.
    pub(super) fn minor_section(&self, text: &mut String) {
        text.push_str(
            "\n## Minor signal\n\nThe reported policy and each policy with another \
             `region.minor_signal` rule, against the same policy without the rule: corpus-weighted \
             A without and with the rule and its cost, compared samples the rule alone moves from \
             pass to review (moved) and the concordant ones among them (burden), primary leaks \
             without and with the rule (one-sided 95 % upper bound), IUPAC-compatible leaks, and \
             truth IUPAC positions inside regions that pass without the rule that the rule marks \
             itself (recall).\n\n\
             | Policy | Rule | A without | A | Cost | Moved | Burden | Leaks without | Leaks | Upper | IUPAC-compatible leaks | IUPAC sites caught |\n\
             |---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n",
        );
        let primary = Measure::Primary.index();
        let compatible = Measure::IupacCompatible.index();
        for (policy, rule, effect) in self.minor_rows() {
            let (without, with) = (effect.leaks_without[primary], effect.leaks_with[primary]);
            let _ = writeln!(
                text,
                "| {} | `{}` | {:.4} | {:.4} | {:.4} | {} | {} | {}/{} | {}/{} | {:.4} | {}/{} | {}/{} |",
                policy.name,
                rule.label(),
                effect.pass_without,
                effect.pass_with,
                effect.pass_without - effect.pass_with,
                effect.moved,
                effect.burden,
                without.leaks,
                without.passed,
                with.leaks,
                with.passed,
                upper_bound(with.leaks, with.passed),
                effect.leaks_with[compatible].leaks,
                effect.leaks_with[compatible].passed,
                effect.caught,
                effect.sites,
            );
        }
    }
}
