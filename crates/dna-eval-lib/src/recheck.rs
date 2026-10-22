//! The re-check list: difference blocks whose truth is under re-check.
//!
//! Raw-trace adjudication finds some unreviewed blocks to be the truth's
//! error, not DNA's, and some samples await a blind re-check of the truth.
//! `run --recheck` names them so that the pass measures count such a block in
//! its own class ([`crate::workload::Class::Recheck`]) instead of as a DNA
//! error. The list is a TSV file with the header
//! `sample_id\tfirst\tlast\treason\tsource`: per row a sample, the 1-based
//! inclusive reference positions of the block, why its truth is under
//! re-check, and where the adjudication is recorded.

use std::collections::{BTreeMap, BTreeSet};

/// The header of a re-check list.
pub(crate) const HEADER: &str = "sample_id\tfirst\tlast\treason\tsource";

/// A parsed re-check list: per sample, its listed spans.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct Recheck {
    spans: BTreeMap<String, Vec<(i64, i64)>>,
    rows: usize,
}

impl Recheck {
    /// Parses a re-check list: [`HEADER`], then one row per block with every
    /// cell non-empty and `1 <= first <= last`. Blank lines, other columns,
    /// and a repeated (`sample_id`, `first`, `last`) are errors.
    pub(crate) fn parse(text: &str) -> Result<Self, String> {
        let mut lines = text.lines().enumerate();
        match lines.next() {
            Some((_, header)) if header == HEADER => {}
            Some((_, header)) => {
                return Err(format!(
                    "line 1: expected the header `{}`, got `{header}`",
                    HEADER.replace('\t', "<TAB>")
                ));
            }
            None => return Err("no header".to_owned()),
        }
        let mut seen = BTreeSet::new();
        let mut spans = BTreeMap::<String, Vec<(i64, i64)>>::new();
        let mut rows = 0;
        for (index, line) in lines {
            let number = index + 1;
            let cells = line.split('\t').collect::<Vec<_>>();
            let [sample, first, last, reason, source] = cells.as_slice() else {
                return Err(format!(
                    "line {number}: expected 5 cells, got {}",
                    cells.len()
                ));
            };
            if [sample, reason, source]
                .iter()
                .any(|cell| cell.trim().is_empty())
            {
                return Err(format!(
                    "line {number}: sample_id, reason, and source must not be empty"
                ));
            }
            let position = |cell: &str| cell.parse::<i64>().ok().filter(|value| *value >= 1);
            let (Some(first), Some(last)) = (position(first), position(last)) else {
                return Err(format!(
                    "line {number}: first and last must be positions from 1, got `{first}` and \
                     `{last}`"
                ));
            };
            if first > last {
                return Err(format!("line {number}: first {first} is after last {last}"));
            }
            if !seen.insert(((*sample).to_owned(), first, last)) {
                return Err(format!(
                    "line {number}: {sample} {first}-{last} is listed twice"
                ));
            }
            spans
                .entry((*sample).to_owned())
                .or_default()
                .push((first, last));
            rows += 1;
        }
        Ok(Self { spans, rows })
    }

    /// The rows of the list.
    pub(crate) const fn rows(&self) -> usize {
        self.rows
    }

    /// The spans listed for `sample`, in file order.
    pub(crate) fn spans(&self, sample: &str) -> &[(i64, i64)] {
        self.spans.get(sample).map_or(&[], Vec::as_slice)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn list(rows: &[&str]) -> String {
        let mut text = format!("{HEADER}\n");
        for row in rows {
            text.push_str(row);
            text.push('\n');
        }
        text
    }

    #[test]
    fn a_list_gives_each_sample_its_spans() -> Result<(), String> {
        let recheck = Recheck::parse(&list(&[
            "S1\t573\t573\ttruth wrong: both strands C8\tresults/a.tsv",
            "S2\t303\t315\tblind re-check pending\tresults/b.tsv",
            "S1\t16057\t16057\ttruth N\tresults/a.tsv",
        ]))?;
        assert_eq!(recheck.rows(), 3);
        assert_eq!(recheck.spans("S1"), [(573, 573), (16057, 16057)]);
        assert_eq!(recheck.spans("S2"), [(303, 315)]);
        assert!(recheck.spans("S3").is_empty());
        assert_eq!(Recheck::parse(&list(&[])).map(|list| list.rows()), Ok(0));
        Ok(())
    }

    #[test]
    fn a_malformed_list_is_rejected() {
        for (text, expected) in [
            (String::new(), "no header"),
            (
                "sample_id\tfirst\tlast\treason\tsource\tnote\n".to_owned(),
                "line 1: expected the header",
            ),
            (
                "sample_id\tfirst\tlast\treason\n".to_owned(),
                "line 1: expected the header",
            ),
            (
                list(&["S1\t573\t573\ttruth wrong"]),
                "line 2: expected 5 cells",
            ),
            (list(&["S1\t573\t573\ta\tb\tc"]), "line 2: expected 5 cells"),
            (list(&[""]), "line 2: expected 5 cells"),
            (list(&["S1\t57x\t573\ta\tb"]), "line 2: first and last"),
            (list(&["S1\t0\t573\ta\tb"]), "line 2: first and last"),
            (list(&["S1\t574\t573\ta\tb"]), "line 2: first 574 is after"),
            (list(&["\t573\t573\ta\tb"]), "line 2: sample_id, reason"),
            (list(&["S1\t573\t573\t \tb"]), "line 2: sample_id, reason"),
            (
                list(&["S1\t573\t573\ta\tb", "S1\t573\t573\tc\td"]),
                "line 3: S1 573-573 is listed twice",
            ),
        ] {
            let message = Recheck::parse(&text).err();
            assert!(
                message
                    .as_deref()
                    .is_some_and(|message| message.starts_with(expected)),
                "{text:?}: {message:?}"
            );
        }
    }
}
