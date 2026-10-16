//! `dna-eval power`: the passes and attainable bounds an approval run can
//! expect, projected from a development run onto a panel's strata.
//!
//! It reads no truth of the panel: only the panel manifest's stratum sizes
//! and the development run's `triage.tsv`, or for a candidate against a base
//! policy its `triage_policies.tsv`. Development pass rates are optimistic,
//! because the rules were tuned there, so the projection also uses each
//! stratum's one-sided 95 % Clopper–Pearson lower bound.

use std::collections::BTreeMap;
use std::fmt::Write as _;

use super::stats::largest_count;
use crate::error::{Error, Result};
use crate::manifest::Manifest;
use crate::triage::score::upper_bound;

/// Largest leak count the table reports.
const MAX_LEAKS: usize = 10;

/// Largest leak count the marginal table reports: a candidate passes few new
/// samples.
const MAX_MARGINAL_LEAKS: usize = 3;

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

/// Per panel stratum, its size, the development rate and its lower bound,
/// and the passes projected at each, as a Markdown table body; with the
/// totals at the rate and at the lower bound.
fn stratum_table(
    dev: &BTreeMap<String, (usize, usize)>,
    sizes: &BTreeMap<&str, usize>,
    text: &mut String,
) -> (usize, usize) {
    let (mut expected, mut cautious) = (0_usize, 0_usize);
    for (stratum, size) in sizes {
        let (samples, passed) = dev.get(*stratum).copied().unwrap_or((0, 0));
        let rate = if samples == 0 {
            0.0
        } else {
            passed as f64 / samples as f64
        };
        let lower = lower_bound(passed, samples);
        let (point, low) = (largest_count(*size, rate), largest_count(*size, lower));
        expected += point;
        cautious += low;
        let _ = writeln!(
            text,
            "| {stratum} | {size} | {rate:.4} | {lower:.4} | {point} | {low} |"
        );
    }
    (expected, cautious)
}

/// The one-sided 95 % lower bound on a rate with `events` of `trials`.
fn lower_bound(events: usize, trials: usize) -> f64 {
    if trials == 0 || events == 0 {
        return 0.0;
    }
    1.0 - upper_bound(trials - events, trials)
}

/// The projection as Markdown.
pub(crate) fn project(triage: &str, panel: &Manifest) -> String {
    let mut text = String::from(
        "# Approval power\n\n\
         | Stratum | Panel samples | Dev pass rate | Dev lower bound | Expected passes | Lower-bound passes |\n\
         |---|---:|---:|---:|---:|---:|\n",
    );
    let (expected, cautious) = stratum_table(&dev_passes(triage), &panel_sizes(panel), &mut text);
    let _ = write!(
        text,
        "\nExpected passes: {expected} (lower-bound projection {cautious}).\n\n"
    );
    leak_table(&mut text, expected, cautious, MAX_LEAKS);
    text
}

/// The projection of the samples a `candidate` policy passes that a `base`
/// policy does not: the samples a marginal endpoint against `base` counts,
/// with the marginal bound attainable at each leak count.
pub(crate) fn project_candidate(
    policies: &str,
    panel: &Manifest,
    base: &str,
    candidate: &str,
) -> Result<String> {
    let dev = dev_new_passes(policies, base, candidate)?;
    let mut text = format!(
        "# Approval power: `{candidate}` against `{base}`\n\n\
         Newly passed samples: passed by `{candidate}` and not by `{base}`.\n\n\
         | Stratum | Panel samples | Dev newly passed rate | Dev lower bound | Expected newly passed | Lower-bound newly passed |\n\
         |---|---:|---:|---:|---:|---:|\n"
    );
    let (expected, cautious) = stratum_table(&dev, &panel_sizes(panel), &mut text);
    let _ = write!(
        text,
        "\nExpected newly passed: {expected} (lower-bound projection {cautious}). The \
         pre-registration records the lower-bound projection: a marginal `upper_bound` \
         below {:.4}, the bound at 0 leaks among {cautious} samples, cannot be met.\n\n",
        upper_bound(0, cautious)
    );
    leak_table(&mut text, expected, cautious, MAX_MARGINAL_LEAKS);
    Ok(text)
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
        assert!(lower_bound(1, 2) < 0.5 && lower_bound(0, 5) == 0.0);
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
}
