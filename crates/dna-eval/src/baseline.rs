//! Baseline runs: the per-sample statuses of an earlier run, for anchors and
//! the diff.
//!
//! A baseline is the `scores.tsv` of an earlier run directory. The
//! `sample_id` and `status` columns are required; `strict_status` and
//! `iupac_compatible_status` are read when present. Statuses fall into
//! classes, so a run that scored on the fair scale (`exact`, `error`) stays
//! comparable with one that reports concordance (`concordant`, `discordant`).

use std::collections::BTreeMap;

use crate::compare::Measure;
use crate::error::BaselineError;

/// Each sample's status under every measure the baseline reports.
pub(crate) type Statuses = BTreeMap<String, BTreeMap<Measure, String>>;

/// Name of the scores file in a run directory.
pub(crate) const SCORES: &str = "scores.tsv";

/// The class of a status, across scoring conventions.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Class {
    /// `exact` or `concordant`.
    Good,
    /// `error` or `discordant`.
    Bad,
    /// `not_applicable`.
    NotApplicable,
    /// `failed`.
    Failed,
    /// Any other status.
    Other,
}

impl Class {
    pub(crate) fn of(status: &str) -> Self {
        match status {
            "exact" | "concordant" => Self::Good,
            "error" | "discordant" => Self::Bad,
            "not_applicable" => Self::NotApplicable,
            "failed" => Self::Failed,
            _ => Self::Other,
        }
    }
}

/// Parses one scores file into statuses by sample and measure.
pub(crate) fn parse(bytes: &[u8]) -> Result<Statuses, BaselineError> {
    let text = std::str::from_utf8(bytes).map_err(|_| BaselineError::NotUtf8)?;
    let mut lines = text.lines();
    let columns = lines
        .next()
        .unwrap_or_default()
        .split('\t')
        .collect::<Vec<_>>();
    let index = |column: &'static str| {
        columns
            .iter()
            .position(|name| *name == column)
            .ok_or(BaselineError::MissingColumn { column })
    };
    let sample = index("sample_id")?;
    index(Measure::Primary.column())?;
    let measures = Measure::ALL
        .into_iter()
        .filter_map(|measure| {
            columns
                .iter()
                .position(|name| *name == measure.column())
                .map(|index| (measure, index))
        })
        .collect::<Vec<_>>();
    let mut statuses = BTreeMap::new();
    for (offset, line) in lines.enumerate() {
        let line_number = offset + 2;
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() != columns.len() {
            return Err(BaselineError::Fields {
                line: line_number,
                found: fields.len(),
                expected: columns.len(),
            });
        }
        let row = measures
            .iter()
            .map(|(measure, index)| (*measure, fields[*index].to_owned()))
            .collect();
        if statuses.insert(fields[sample].to_owned(), row).is_some() {
            return Err(BaselineError::DuplicateSample {
                line: line_number,
                sample: fields[sample].to_owned(),
            });
        }
    }
    Ok(statuses)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_statuses_by_column_name() -> Result<(), BaselineError> {
        let statuses = parse(
            b"sample_id\tstratum\tweight\tstatus\tfp\tfn\terror_regions\n\
              S1\tclean\t1\texact\t0\t0\t\n\
              S2\tclean\t1\tfailed\t\t\t\n",
        )?;
        assert_eq!(statuses["S1"][&Measure::Primary], "exact");
        assert!(!statuses["S1"].contains_key(&Measure::Strict));
        assert_eq!(Class::of(&statuses["S2"][&Measure::Primary]), Class::Failed);
        let statuses = parse(b"sample_id\tstatus\tstrict_status\nS1\tconcordant\tdiscordant\n")?;
        assert_eq!(statuses["S1"][&Measure::Strict], "discordant");
        assert_eq!(Class::of("concordant"), Class::of("exact"));
        assert_eq!(Class::of("discordant"), Class::Bad);
        assert_eq!(Class::of("unknown"), Class::Other);
        Ok(())
    }

    #[test]
    fn rejects_files_it_cannot_compare() {
        assert_eq!(
            parse(b"sample_id\tfp\n"),
            Err(BaselineError::MissingColumn { column: "status" })
        );
        assert_eq!(
            parse(b"sample_id\tstatus\tfp\nS1\texact\n"),
            Err(BaselineError::Fields {
                line: 2,
                found: 2,
                expected: 3,
            })
        );
        assert_eq!(
            parse(b"sample_id\tstatus\nS1\texact\nS1\terror\n"),
            Err(BaselineError::DuplicateSample {
                line: 3,
                sample: "S1".into(),
            })
        );
        assert_eq!(parse(&[0xff]), Err(BaselineError::NotUtf8));
    }
}
