//! `dna-eval power`: the passes and attainable bounds an approval run can
//! expect, projected from a development run onto a panel's strata.
//!
//! It reads no truth of the panel: only the panel manifest's stratum sizes
//! and the development run's `triage.tsv`, or for a candidate against a base
//! policy its `triage_policies.tsv`. The panel's rate is the stratified rate
//! of the development rates, weighted by the panel's stratum sizes. Its
//! lower bound is the Korn–Graubard one-sided 95 % lower bound, the method
//! of the gates, applied once to the whole panel. Development rates come
//! from the panel the rules were tuned on, so even that bound may be
//! optimistic.

use std::collections::BTreeMap;
use std::fmt::Write as _;

use super::stats::{binomial_at_most, korn_graubard_lower, largest_count};
use crate::error::{Error, Result};
use crate::manifest::Manifest;
use crate::triage::score::upper_bound;

/// Largest leak count the table reports.
const MAX_LEAKS: usize = 10;

/// Largest leak count the marginal table reports: a candidate passes few new
/// samples.
const MAX_MARGINAL_LEAKS: usize = 3;

/// A marginal bound whose chance of being met at a true leak rate is asked
/// for.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Target {
    /// The true leak rate assumed among newly passed samples.
    pub(crate) true_rate: f64,
    /// The candidate marginal upper bounds.
    pub(crate) bounds: Vec<f64>,
}

/// Per stratum: active development samples and how many passed.
fn dev_passes(triage: &str) -> BTreeMap<String, (usize, usize)> {
    let mut lines = triage.lines();
    let header = lines
        .next()
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .unwrap_or_default();
    let column = |name: &str| header.iter().position(|field| *field == name);
    let (Some(stratum), Some(retired), Some(outcome)) =
        (column("stratum"), column("retired"), column("outcome"))
    else {
        return BTreeMap::new();
    };
    let mut strata = BTreeMap::<String, (usize, usize)>::new();
    for fields in lines.map(|line| line.split('\t').collect::<Vec<_>>()) {
        if fields.get(retired).is_some_and(|value| *value == "1") {
            continue;
        }
        let Some(name) = fields.get(stratum) else {
            continue;
        };
        let entry = strata.entry((*name).to_owned()).or_default();
        entry.0 += 1;
        entry.1 += usize::from(fields.get(outcome).is_some_and(|value| *value == "pass"));
    }
    strata
}

/// Per stratum: active development samples and how many the `candidate`
/// policy passes that the `base` policy does not, from `triage_policies.tsv`.
fn dev_new_passes(
    policies: &str,
    base: &str,
    candidate: &str,
) -> Result<BTreeMap<String, (usize, usize)>> {
    let mut lines = policies.lines();
    let header = lines
        .next()
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .unwrap_or_default();
    let column = |name: &str| -> Result<usize> {
        let mut found = header
            .iter()
            .enumerate()
            .filter(|(_, field)| **field == name)
            .map(|(index, _)| index);
        match (found.next(), found.next()) {
            (Some(index), None) => Ok(index),
            _ => Err(Error::Policy {
                name: name.to_owned(),
            }),
        }
    };
    let (stratum, retired) = (column("stratum")?, column("retired")?);
    let (base, candidate) = (column(base)?, column(candidate)?);
    let mut strata = BTreeMap::<String, (usize, usize)>::new();
    for fields in lines.map(|line| line.split('\t').collect::<Vec<_>>()) {
        if fields.get(retired).is_some_and(|value| *value == "1") {
            continue;
        }
        let Some(name) = fields.get(stratum) else {
            continue;
        };
        let passes = |index: usize| fields.get(index).is_some_and(|value| *value == "pass");
        let entry = strata.entry((*name).to_owned()).or_default();
        entry.0 += 1;
        entry.1 += usize::from(passes(candidate) && !passes(base));
    }
    Ok(strata)
}

/// Active panel samples per stratum.
fn panel_sizes(panel: &Manifest) -> BTreeMap<&str, usize> {
    let mut sizes = BTreeMap::<&str, usize>::new();
    for sample in panel.samples.iter().filter(|sample| !sample.retired) {
        *sizes.entry(sample.stratum.as_str()).or_default() += 1;
    }
    sizes
}

/// The development rates projected onto a panel.
#[derive(Debug, Clone, PartialEq)]
struct Projection {
    /// Active panel samples.
    size: usize,
    /// The stratified rate and its one-sided 95 % Korn–Graubard lower bound.
    rate: f64,
    lower: f64,
    /// The panel samples counted at the rate and at the lower bound.
    expected: usize,
    cautious: usize,
}

/// Projects development counts per stratum onto the panel's strata, writing
/// one Markdown table row per panel stratum to `text`. A panel stratum
/// without development samples counts at rate zero and is marked.
fn project_rates(
    dev: &BTreeMap<String, (usize, usize)>,
    sizes: &BTreeMap<&str, usize>,
    text: &mut String,
) -> Projection {
    let size = sizes.values().sum::<usize>();
    let (mut rate, mut variance, mut samples) = (0.0, 0.0, 0_usize);
    for (stratum, panel) in sizes {
        let (count, events) = dev.get(*stratum).copied().unwrap_or((0, 0));
        let share = *panel as f64 / size as f64;
        let stratum_rate = if count == 0 {
            0.0
        } else {
            events as f64 / count as f64
        };
        if count > 0 {
            rate += share * stratum_rate;
            variance += share * share * stratum_rate * (1.0 - stratum_rate) / count as f64;
            samples += count;
        }
        let note = if count == 0 {
            " (no development samples)"
        } else {
            ""
        };
        let _ = writeln!(
            text,
            "| {stratum}{note} | {panel} | {count} | {events} | {stratum_rate:.4} | {} |",
            largest_count(*panel, stratum_rate)
        );
    }
    let lower = korn_graubard_lower(rate, variance, samples);
    Projection {
        size,
        rate,
        lower,
        expected: largest_count(size, rate),
        cautious: largest_count(size, lower),
    }
}

/// The projection as Markdown.
pub(crate) fn project(triage: &str, panel: &Manifest) -> String {
    let mut text = String::from(
        "# Approval power\n\n\
         | Stratum | Panel samples | Dev samples | Dev passes | Dev pass rate | Expected passes |\n\
         |---|---:|---:|---:|---:|---:|\n",
    );
    let projection = project_rates(&dev_passes(triage), &panel_sizes(panel), &mut text);
    summary(&mut text, &projection, "passes");
    leak_table(
        &mut text,
        projection.expected,
        projection.cautious,
        MAX_LEAKS,
    );
    text
}

/// The projection of the samples a `candidate` policy passes that a `base`
/// policy does not: the samples a marginal endpoint against `base` counts,
/// with the marginal bound attainable at each leak count and, for a
/// `target`, the chance that each of its bounds is met.
pub(crate) fn project_candidate(
    policies: &str,
    panel: &Manifest,
    base: &str,
    candidate: &str,
    target: Option<&Target>,
) -> Result<String> {
    let dev = dev_new_passes(policies, base, candidate)?;
    let mut text = format!(
        "# Approval power: `{candidate}` against `{base}`\n\n\
         Newly passed samples: passed by `{candidate}` and not by `{base}`.\n\n\
         | Stratum | Panel samples | Dev samples | Dev newly passed | Dev rate | Expected newly passed |\n\
         |---|---:|---:|---:|---:|---:|\n"
    );
    let projection = project_rates(&dev, &panel_sizes(panel), &mut text);
    summary(&mut text, &projection, "newly passed samples");
    let _ = write!(
        text,
        "The pre-registration records the lower-bound projection, {} samples: a marginal \
         `upper_bound` below {:.4}, the bound at 0 leaks among them, cannot be met.\n\n",
        projection.cautious,
        upper_bound(0, projection.cautious)
    );
    leak_table(
        &mut text,
        projection.expected,
        projection.cautious,
        MAX_MARGINAL_LEAKS,
    );
    if let Some(target) = target {
        meet_table(&mut text, projection.cautious, target);
    }
    Ok(text)
}

/// The panel totals of a projection.
fn summary(text: &mut String, projection: &Projection, counted: &str) {
    let _ = write!(
        text,
        "\nPanel: {} active samples. Stratified rate {:.4}, Korn–Graubard one-sided 95 % \
         lower bound {:.4}. Expected {counted}: {} (lower-bound projection {}). The \
         development rates come from the panel the rules were tuned on, so even the lower \
         bound may be optimistic.\n\n",
        projection.size,
        projection.rate,
        projection.lower,
        projection.expected,
        projection.cautious
    );
}

/// The count-based bound attainable at 0 to `maximum` leaks among the
/// expected and the lower-bound passes, as a Markdown table.
fn leak_table(text: &mut String, expected: usize, cautious: usize, maximum: usize) {
    text.push_str(
        "| Leaks | Count UB at expected passes | Count UB at lower-bound passes |\n\
         |---:|---:|---:|\n",
    );
    for leaks in 0..=maximum {
        let _ = writeln!(
            text,
            "| {leaks} | {:.4} | {:.4} |",
            upper_bound(leaks, expected),
            upper_bound(leaks, cautious)
        );
    }
}

/// The most leaks among `samples` whose count-based bound meets `bound`.
fn allowed_leaks(samples: usize, bound: f64) -> Option<usize> {
    (0..=samples)
        .take_while(|leaks| upper_bound(*leaks, samples) <= bound)
        .last()
}

/// Per bound of `target`, the leaks it allows among `samples` and the
/// chance of at most that many at the true rate, as a Markdown table.
fn meet_table(text: &mut String, samples: usize, target: &Target) {
    let _ = write!(
        text,
        "\nChance that a marginal bound is met among the {samples} lower-bound samples at a \
         true leak rate of {:.4}:\n\n\
         | Bound | Leaks allowed | P(met) |\n\
         |---:|---:|---:|\n",
        target.true_rate
    );
    for bound in &target.bounds {
        match allowed_leaks(samples, *bound) {
            Some(leaks) => {
                let _ = writeln!(
                    text,
                    "| {bound:.4} | {leaks} | {:.4} |",
                    binomial_at_most(leaks, samples, target.true_rate)
                );
            }
            None => {
                let _ = writeln!(text, "| {bound:.4} | none | 0.0000 |");
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts_active_passes_per_stratum() {
        let triage = "sample_id\tstratum\tweight_corpus\tretired\toutcome\n\
                      a\tclean\t1\t0\tpass\n\
                      b\tclean\t1\t0\treview\n\
                      c\tclean\t1\t1\tpass\n\
                      d\tiupac\t1\t0\tpass\n";
        let strata = dev_passes(triage);
        assert_eq!(strata["clean"], (2, 1));
        assert_eq!(strata["iupac"], (1, 1));
    }

    #[test]
    fn counts_only_samples_the_candidate_newly_passes() -> Result<()> {
        let policies = "sample_id\tstratum\tretired\tdraft\twide\n\
                        a\tclean\t0\tpass\tpass\n\
                        b\tclean\t0\treview\tpass\n\
                        c\tclean\t0\treview\treview\n\
                        d\tclean\t1\treview\tpass\n\
                        e\tiupac\t0\tpass\treview\n";
        let strata = dev_new_passes(policies, "draft", "wide")?;
        assert_eq!(strata["clean"], (3, 1));
        assert_eq!(strata["iupac"], (1, 0));
        assert!(matches!(
            dev_new_passes(policies, "draft", "missing"),
            Err(Error::Policy { name }) if name == "missing"
        ));
        let repeated = "sample_id\tstratum\tretired\tdraft\tdraft\n";
        assert!(dev_new_passes(repeated, "draft", "draft").is_err());
        Ok(())
    }

    #[test]
    fn projects_one_joint_lower_bound_over_the_panel() {
        // Two strata at 10 % of 100 development samples each, onto 1,000
        // panel samples each; a third panel stratum has no development
        // samples and counts at zero.
        let dev = BTreeMap::from([("a".to_owned(), (100, 10)), ("b".to_owned(), (100, 10))]);
        let sizes = BTreeMap::from([("a", 1000), ("b", 1000), ("c", 0)]);
        let mut text = String::new();
        let projection = project_rates(&dev, &sizes, &mut text);
        assert!((projection.rate - 0.1).abs() < 1e-12);
        assert_eq!(projection.expected, 200);
        // The joint bound is the Clopper–Pearson bound of 20 in 200, above
        // the sum of the two strata's own bounds.
        let joint = super::super::stats::beta_lower(20.0, 200.0);
        assert!((projection.lower - joint).abs() < 1e-9);
        let summed = 2 * largest_count(1000, super::super::stats::beta_lower(10.0, 100.0));
        assert!(projection.cautious > summed);
        assert!(text.contains("| c (no development samples) | 0 | 0 | 0 | 0.0000 | 0 |"));
    }

    #[test]
    fn a_bound_is_met_with_the_leaks_it_allows() {
        // 300 samples: 0 leaks give a bound of about 0.99 %, 1 leak 1.57 %.
        assert_eq!(allowed_leaks(300, 0.01), Some(0));
        assert_eq!(allowed_leaks(300, 0.016), Some(1));
        assert_eq!(allowed_leaks(300, 0.005), None);
        let mut text = String::new();
        let target = Target {
            true_rate: 0.002,
            bounds: vec![0.005, 0.016],
        };
        meet_table(&mut text, 300, &target);
        let met = binomial_at_most(1, 300, 0.002);
        assert!(text.contains("| 0.0050 | none | 0.0000 |"));
        assert!(text.contains(&format!("| 0.0160 | 1 | {met:.4} |")));
    }
}
