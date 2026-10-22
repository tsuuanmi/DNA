//! The version ledger of the regression set: the local index of recorded
//! entries (`index.tsv`, beside a copy of each entry's check under
//! `results/`), and the rows `record` appends to the two committed,
//! append-only baseline files, which name no sample.

use std::fmt::Write as _;
use std::path::Path;

use super::Header;

/// The local index of a ledger directory.
pub(crate) const INDEX: &str = "index.tsv";
/// The header of `index.tsv`.
pub(crate) const INDEX_HEADER: &str =
    "entry\tkind\tbase\trevision\tgolden_sha256\tverdict\trecorded\n";

/// The columns of the committed ledger baseline, one row per entry.
pub(crate) const LEDGER_COLUMNS: &str = "entry\tdate\tkind\tbase\trevision\tdna_version\tplugins\t\
     contracts\tconfiguration_sha256\tprofile_sha256\treference_sha256\tmanifest_body_sha256\t\
     truth_digest\trecheck_sha256\tgolden_sha256\tkinds_sha256\tsamples\tfailed\tclear\t\
     with_check\tchecks\tnotes\talleles\tblocks_non_het\tescapes\tblocks_het\t\
     blocks_length_note\tblocks_truth_n\tblocks_recheck\tchanged_samples\tchanged_lines\t\
     explained_lines\tposthoc_rows\tintent_blob\tverdict";
/// The columns of the committed strata baseline, one row per entry and
/// stratum.
pub(crate) const STRATA_COLUMNS: &str = "entry\tstratum\tsamples\tgolden_sha256\tchanged_samples";

/// The comment line that opens a baseline file of `series`: `kind` is
/// `regress-ledger` or `regress-strata`.
pub(crate) fn baseline_comment(kind: &str, series: &str) -> String {
    format!("# {kind} series={series}")
}

/// One row of `index.tsv`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct IndexRow {
    pub(crate) entry: String,
    pub(crate) kind: String,
    pub(crate) base: String,
    pub(crate) revision: String,
    pub(crate) golden_sha256: String,
    pub(crate) verdict: String,
    pub(crate) recorded: String,
}

impl IndexRow {
    /// The row's line, with its newline.
    pub(crate) fn line(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            self.entry,
            self.kind,
            self.base,
            self.revision,
            self.golden_sha256,
            self.verdict,
            self.recorded
        )
    }
}

/// The rows of `index.tsv`, in recording order: none when the file is
/// absent.
pub(crate) fn read_index(text: Option<&str>) -> Result<Vec<IndexRow>, String> {
    let Some(text) = text else {
        return Ok(Vec::new());
    };
    let body = text
        .strip_prefix(INDEX_HEADER)
        .ok_or("index.tsv does not start with its header")?;
    body.lines()
        .enumerate()
        .map(|(number, line)| {
            let cells = line.split('\t').collect::<Vec<_>>();
            let [
                entry,
                kind,
                base,
                revision,
                golden_sha256,
                verdict,
                recorded,
            ] = cells[..]
            else {
                return Err(format!(
                    "index.tsv line {} has {} cells",
                    number + 2,
                    cells.len()
                ));
            };
            if !valid_entry(entry) {
                return Err(format!(
                    "index.tsv line {} names entry {entry:?}",
                    number + 2
                ));
            }
            Ok(IndexRow {
                entry: entry.to_owned(),
                kind: kind.to_owned(),
                base: base.to_owned(),
                revision: revision.to_owned(),
                golden_sha256: golden_sha256.to_owned(),
                verdict: verdict.to_owned(),
                recorded: recorded.to_owned(),
            })
        })
        .collect()
}

/// Whether `entry` is an entry id: `r` and at least four digits.
pub(crate) fn valid_entry(entry: &str) -> bool {
    entry
        .strip_prefix('r')
        .is_some_and(|digits| digits.len() >= 4 && digits.bytes().all(|byte| byte.is_ascii_digit()))
}

/// The entries of a committed baseline file of `series` with `columns`, in
/// file order (an entry repeats in the strata file, once per stratum).
pub(crate) fn baseline_entries(
    text: &str,
    (kind, series): (&str, &str),
    columns: &str,
) -> Result<Vec<String>, String> {
    let mut lines = text
        .strip_suffix('\n')
        .ok_or("the file does not end in a newline")?
        .split('\n');
    let comment = baseline_comment(kind, series);
    if lines.next() != Some(comment.as_str()) {
        return Err(format!("the first line is not {comment:?}"));
    }
    if lines.next() != Some(columns) {
        return Err("the column header is not the contract's".to_owned());
    }
    let width = columns.split('\t').count();
    lines
        .map(|line| {
            let cells = line.split('\t').collect::<Vec<_>>();
            if cells.len() != width {
                return Err(format!("a row has {} cells, not {width}", cells.len()));
            }
            Ok(cells[0].to_owned())
        })
        .collect()
}

/// The ledger baseline row of `entry`, recorded on `date`, from the run
/// header `header` of its check.
pub(crate) fn ledger_row(entry: &str, (date, kind): (&str, &str), header: &Header) -> String {
    let counts = &header.counts;
    let changes = &header.changes;
    let key = &header.version_key;
    let kinds = header
        .kinds
        .iter()
        .map(|kind| format!("{}={}", kind.kind, kind.sha256))
        .collect::<Vec<_>>()
        .join(";");
    let cells = [
        entry.to_owned(),
        date.to_owned(),
        kind.to_owned(),
        header.head.clone().unwrap_or_else(|| "-".to_owned()),
        header.source.revision.clone().unwrap_or_default(),
        key.dna_version.clone(),
        key.plugins.clone(),
        key.contracts.clone(),
        key.configuration_sha256.clone(),
        header.profile_sha256.clone(),
        header.reference_sha256.clone(),
        header.manifest.body_sha256.clone(),
        header.truth.digest.clone(),
        header
            .recheck
            .as_ref()
            .map_or_else(|| "-".to_owned(), |recheck| recheck.sha256.clone()),
        header.golden_sha256.clone(),
        kinds,
        counts.samples.to_string(),
        counts.failed.to_string(),
        counts.clear.to_string(),
        counts.with_check.to_string(),
        counts.checks.to_string(),
        counts.notes.to_string(),
        counts.alleles.to_string(),
        counts.blocks_non_het.to_string(),
        counts.escapes.to_string(),
        counts.blocks_het.to_string(),
        counts.blocks_length_note.to_string(),
        counts.blocks_truth_n.to_string(),
        counts.blocks_recheck.to_string(),
        changes.changed_samples.to_string(),
        changes.changed_lines.to_string(),
        changes.explained_lines.to_string(),
        changes.posthoc_rows.to_string(),
        header
            .intent
            .as_ref()
            .map_or_else(|| "-".to_owned(), |intent| intent.blob.clone()),
        header.verdict.clone(),
    ];
    let mut line = cells.join("\t");
    line.push('\n');
    line
}

/// The strata baseline rows of `entry`, one per stratum in name order.
pub(crate) fn strata_rows(entry: &str, header: &Header) -> String {
    let mut text = String::new();
    for stratum in &header.strata {
        let _ = writeln!(
            text,
            "{entry}\t{}\t{}\t{}\t{}",
            stratum.stratum, stratum.samples, stratum.golden_sha256, stratum.changed_samples
        );
    }
    text
}

/// The ledger directory's index text, when it has one.
pub(crate) fn index_text(ledger: &Path) -> Result<Option<String>, String> {
    let path = ledger.join(INDEX);
    match std::fs::read_to_string(&path) {
        Ok(text) => Ok(Some(text)),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(None),
        Err(error) => Err(format!("{}: {error}", path.display())),
    }
}

/// The UTC calendar date `YYYY-MM-DD` of `seconds` since the Unix epoch.
pub(crate) fn date(seconds: u64) -> String {
    // Days to the civil date, after H. Hinnant's `civil_from_days`.
    let days = i64::try_from(seconds / 86_400).unwrap_or(i64::MAX) + 719_468;
    let era = days.div_euclid(146_097);
    let day_of_era = days.rem_euclid(146_097);
    let year_of_era =
        (day_of_era - day_of_era / 1_460 + day_of_era / 36_524 - day_of_era / 146_096) / 365;
    let day_of_year = day_of_era - (365 * year_of_era + year_of_era / 4 - year_of_era / 100);
    let shifted_month = (5 * day_of_year + 2) / 153;
    let day = day_of_year - (153 * shifted_month + 2) / 5 + 1;
    let month = if shifted_month < 10 {
        shifted_month + 3
    } else {
        shifted_month - 9
    };
    let year = year_of_era + era * 400 + i64::from(month <= 2);
    format!("{year:04}-{month:02}-{day:02}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dates_are_utc_calendar_days() {
        assert_eq!(date(0), "1970-01-01");
        assert_eq!(date(951_782_400), "2000-02-29");
        assert_eq!(date(1_792_627_200 - 1), "2026-10-21");
        assert_eq!(date(1_792_627_200), "2026-10-22");
    }

    #[test]
    fn entries_are_r_and_four_digits_or_more() {
        for (entry, valid) in [
            ("r0001", true),
            ("r12345", true),
            ("r001", false),
            ("0001", false),
            ("r00a1", false),
        ] {
            assert_eq!(valid_entry(entry), valid, "{entry}");
        }
    }

    #[test]
    fn the_committed_baselines_follow_the_format() -> Result<(), String> {
        baseline_entries(
            include_str!("../../baselines/regress-v1-ledger.tsv"),
            ("regress-ledger", "regress-v1"),
            LEDGER_COLUMNS,
        )?;
        baseline_entries(
            include_str!("../../baselines/regress-v1-strata.tsv"),
            ("regress-strata", "regress-v1"),
            STRATA_COLUMNS,
        )?;
        Ok(())
    }

    #[test]
    fn reads_the_index_and_the_baselines_strictly() -> Result<(), String> {
        assert!(read_index(None)?.is_empty());
        let row = IndexRow {
            entry: "r0001".into(),
            kind: "initial".into(),
            base: "-".into(),
            revision: "a".repeat(40),
            golden_sha256: "b".repeat(64),
            verdict: "REGRESS INITIAL".into(),
            recorded: "2026-10-22".into(),
        };
        let text = format!("{INDEX_HEADER}{}", row.line());
        assert_eq!(read_index(Some(&text))?, vec![row]);
        assert!(read_index(Some("entry\n")).is_err());
        assert!(read_index(Some(&format!("{INDEX_HEADER}x\t-\n"))).is_err());

        let head = format!("# regress-strata series=regress-v1\n{STRATA_COLUMNS}\n");
        let series = ("regress-strata", "regress-v1");
        assert!(baseline_entries(&head, series, STRATA_COLUMNS)?.is_empty());
        assert_eq!(
            baseline_entries(
                &format!("{head}r0001\ta\t2\tx\t0\nr0001\tb\t1\ty\t0\n"),
                series,
                STRATA_COLUMNS
            )?,
            vec!["r0001", "r0001"]
        );
        assert!(baseline_entries(&head, ("regress-strata", "regress-v2"), STRATA_COLUMNS).is_err());
        assert!(baseline_entries(&format!("{head}r0001\ta\n"), series, STRATA_COLUMNS).is_err());
        assert!(baseline_entries(head.trim_end(), series, STRATA_COLUMNS).is_err());
        Ok(())
    }
}
