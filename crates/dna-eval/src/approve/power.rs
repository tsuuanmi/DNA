//! `dna-eval power`: the passes and attainable bounds an approval run can
//! expect, projected from a development run onto a panel's strata.
//!
//! It reads no truth of the panel: only the panel manifest's stratum sizes
//! and the development run's `triage.tsv`. Development pass rates are
//! optimistic, because the rules were tuned there, so the projection also
//! uses each stratum's one-sided 95 % Clopper–Pearson lower bound.

use std::collections::BTreeMap;
use std::fmt::Write as _;

use super::stats::largest_count;
use crate::manifest::Manifest;
use crate::triage::score::upper_bound;

/// Largest leak count the table reports.
const MAX_LEAKS: usize = 10;

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

/// The one-sided 95 % lower bound on a rate with `events` of `trials`.
fn lower_bound(events: usize, trials: usize) -> f64 {
    if trials == 0 || events == 0 {
        return 0.0;
    }
    1.0 - upper_bound(trials - events, trials)
}

/// The projection as Markdown.
pub(crate) fn project(triage: &str, panel: &Manifest) -> String {
    let dev = dev_passes(triage);
    let mut sizes = BTreeMap::<&str, usize>::new();
    for sample in panel.samples.iter().filter(|sample| !sample.retired) {
        *sizes.entry(sample.stratum.as_str()).or_default() += 1;
    }
    let mut text = String::from(
        "# Approval power\n\n\
         | Stratum | Panel samples | Dev pass rate | Dev lower bound | Expected passes | Lower-bound passes |\n\
         |---|---:|---:|---:|---:|---:|\n",
    );
    let (mut expected, mut cautious) = (0_usize, 0_usize);
    for (stratum, size) in &sizes {
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
    let _ = write!(
        text,
        "\nExpected passes: {expected} (lower-bound projection {cautious}).\n\n\
         | Leaks | Count UB at expected passes | Count UB at lower-bound passes |\n\
         |---:|---:|---:|\n"
    );
    for leaks in 0..=MAX_LEAKS {
        let _ = writeln!(
            text,
            "| {leaks} | {:.4} | {:.4} |",
            upper_bound(leaks, expected),
            upper_bound(leaks, cautious)
        );
    }
    text
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
}
