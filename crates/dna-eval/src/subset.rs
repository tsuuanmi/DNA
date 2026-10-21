//! A quick subset of a panel: a deterministic, stratified sample of its
//! staying rows, written as a manifest of its own, for a first-pass
//! `dna-eval diff` before the full panel. Its verdict is a pre-filter, never
//! the gate: a candidate that passes it still has to pass the full panel.
//!
//! Each stratum keeps its share of `rows`, at least `minimum` rows (all of
//! them when it has fewer), so the rare strata are always present. Within a
//! stratum the rows with the smallest SHA-256 of `salt`, a tab, and the
//! sample identifier are kept. The kept rows are written verbatim and in the
//! parent's order; the header keeps the parent's fields, records the
//! subset's parameters and parent, and recomputes the row count, the body
//! hash, and the truth digest over the kept rows.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

use sha2::{Digest, Sha256};

use crate::manifest::{Manifest, Sample, encode};
use crate::truth;

/// The parent fields a subset's header carries over unchanged.
const CARRIED: [&str; 8] = [
    "panel",
    "version",
    "salt",
    "created",
    "truth",
    "truth_database_sha256",
    "truth_created",
    "truth_samples",
];

/// A subset manifest and how many rows each stratum kept.
pub(crate) struct Subset {
    pub(crate) text: String,
    pub(crate) rows: usize,
    pub(crate) strata: BTreeMap<String, (usize, usize)>,
}

/// The rows each stratum keeps: its proportional share of `rows`, at least
/// `minimum`, at most all of its `available` rows.
fn allocation(available: usize, total: usize, rows: usize, minimum: usize) -> usize {
    let share = (available * rows).div_ceil(total.max(1));
    share.max(minimum).min(available)
}

/// The quick subset of `parent`'s staying rows.
pub(crate) fn subset(parent: &Manifest, rows: usize, minimum: usize, salt: &str) -> Subset {
    let staying = parent
        .samples
        .iter()
        .filter(|sample| !sample.retired)
        .collect::<Vec<_>>();
    let mut by_stratum = BTreeMap::<&str, Vec<&Sample>>::new();
    for sample in &staying {
        by_stratum
            .entry(sample.stratum.as_str())
            .or_default()
            .push(sample);
    }
    let rank = |sample: &Sample| {
        format!(
            "{:x}",
            Sha256::digest(format!("{salt}\t{}", sample.id).as_bytes())
        )
    };
    let mut kept = BTreeSet::<&str>::new();
    let mut strata = BTreeMap::new();
    for (stratum, samples) in &mut by_stratum {
        let take = allocation(samples.len(), staying.len(), rows, minimum);
        samples.sort_by_cached_key(|sample| rank(sample));
        kept.extend(samples.iter().take(take).map(|sample| sample.id.as_str()));
        strata.insert((*stratum).to_owned(), (take, samples.len()));
    }
    let chosen = parent
        .samples
        .iter()
        .filter(|sample| kept.contains(sample.id.as_str()))
        .collect::<Vec<_>>();
    let mut body = String::new();
    for sample in &chosen {
        body.push_str(&sample.fields.join("\t"));
        body.push('\n');
    }
    let body_sha256 = format!("{:x}", Sha256::digest(body.as_bytes()));
    let truth_digest = truth::panel_digest(chosen.iter().map(|sample| {
        (
            sample.id.as_str(),
            sample.truth_digest.as_deref().unwrap_or("-"),
        )
    }));
    let mut text = String::new();
    let carried = CARRIED
        .iter()
        .filter_map(|key| {
            parent
                .header
                .get(*key)
                .map(|value| format!("{key}={}", encode(value)))
        })
        .collect::<Vec<_>>()
        .join(" ");
    let _ = writeln!(text, "# {carried} rows={} retired=0", chosen.len());
    let _ = writeln!(
        text,
        "# subset=quick subset_rows={rows} subset_minimum={minimum} subset_salt={} \
         parent_body_sha256={} truth_digest={truth_digest}",
        encode(salt),
        parent.body_sha256,
    );
    let _ = writeln!(
        text,
        "# body_sha256={body_sha256} (sha256 of the data rows below the column header)"
    );
    text.push_str(&parent.columns.join("\t"));
    text.push('\n');
    text.push_str(&body);
    Subset {
        text,
        rows: chosen.len(),
        strata,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::manifest::parse;

    fn manifest(rows: &[(&str, &str, &str)]) -> Result<Manifest, Box<dyn std::error::Error>> {
        let columns = "sample_id\tbatch\tstratum\tweight\tweight_corpus\tanchor_candidate\ttrace_paths\tretired\n";
        let mut body = String::new();
        for (id, stratum, retired) in rows {
            let _ = writeln!(body, "{id}\tb\t{stratum}\t1\t1\t0\t/t/{id}.ab1\t{retired}");
        }
        let text = format!(
            "# panel=dev version=v9 salt=s rows={}\n# body_sha256={:x}\n{columns}{body}",
            rows.len(),
            Sha256::digest(body.as_bytes())
        );
        Ok(parse(text.as_bytes())?)
    }

    #[test]
    fn keeps_every_stratum_and_parses_as_a_manifest() -> Result<(), Box<dyn std::error::Error>> {
        let mut rows = (0..90)
            .map(|index| (format!("C{index}"), "clean", "0"))
            .collect::<Vec<_>>();
        rows.extend((0..2).map(|index| (format!("R{index}"), "rare", "0")));
        rows.push(("X0".to_owned(), "rare", "1"));
        let borrowed = rows
            .iter()
            .map(|(id, stratum, retired)| (id.as_str(), *stratum, *retired))
            .collect::<Vec<_>>();
        let parent = manifest(&borrowed)?;
        let quick = subset(&parent, 10, 5, "q");
        // Clean keeps its share, rare all its staying rows, never the retired one.
        assert_eq!(quick.strata.get("clean"), Some(&(10, 90)));
        assert_eq!(quick.strata.get("rare"), Some(&(2, 2)));
        assert_eq!(quick.rows, 12);
        let reread = parse(quick.text.as_bytes())?;
        assert_eq!(reread.samples.len(), 12);
        assert!(reread.samples.iter().all(|sample| sample.id != "X0"));
        assert_eq!(reread.header.get("panel").map(String::as_str), Some("dev"));
        assert_eq!(
            reread.header.get("subset_minimum").map(String::as_str),
            Some("5")
        );
        // The same parameters keep the same rows; another salt others.
        assert_eq!(subset(&parent, 10, 5, "q").text, quick.text);
        assert_ne!(subset(&parent, 10, 5, "r").text, quick.text);
        Ok(())
    }
}
