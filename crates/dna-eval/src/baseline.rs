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

/// The measures a baseline reports, and each sample's status under them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Statuses {
    /// Measures with a status column, in [`Measure::ALL`] order.
    pub(crate) measures: Vec<Measure>,
    /// Each sample's status under every reported measure.
    pub(crate) samples: BTreeMap<String, BTreeMap<Measure, String>>,
}

impl Statuses {
    /// The statuses of one sample.
    pub(crate) fn get(&self, sample: &str) -> Option<&BTreeMap<Measure, String>> {
        self.samples.get(sample)
    }

    /// Whether the baseline has a status column for `measure`.
    pub(crate) fn reports(&self, measure: Measure) -> bool {
        self.measures.contains(&measure)
    }
}

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
    let mut samples = BTreeMap::new();
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
        if samples.insert(fields[sample].to_owned(), row).is_some() {
            return Err(BaselineError::DuplicateSample {
                line: line_number,
                sample: fields[sample].to_owned(),
            });
        }
    }
    Ok(Statuses {
        measures: measures.iter().map(|(measure, _)| *measure).collect(),
        samples,
    })
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
        assert_eq!(statuses.samples["S1"][&Measure::Primary], "exact");
        assert!(!statuses.samples["S1"].contains_key(&Measure::Strict));
        assert_eq!(statuses.measures, [Measure::Primary]);
        assert!(!statuses.reports(Measure::Strict));
        assert_eq!(
            Class::of(&statuses.samples["S2"][&Measure::Primary]),
            Class::Failed
        );
        let statuses = parse(b"sample_id\tstatus\tstrict_status\nS1\tconcordant\tdiscordant\n")?;
        assert_eq!(statuses.samples["S1"][&Measure::Strict], "discordant");
        assert_eq!(statuses.measures, [Measure::Primary, Measure::Strict]);
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
