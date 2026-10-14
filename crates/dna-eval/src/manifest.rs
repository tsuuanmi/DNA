//! Panel manifests: the samples of one test panel, their stratum, weights,
//! and traces.
//!
//! A manifest starts with `#` comment lines of `key=value` fields (panel,
//! version, salt, truth and body hashes), then one tab-separated column
//! header, then one row per sample. `body_sha256` is the SHA-256 of every byte
//! after the column header line and must match. From version 2 a `retired`
//! column marks samples kept for continuity but left out of estimates; a
//! manifest without it retires nothing.

use std::collections::{BTreeMap, BTreeSet};
use std::path::PathBuf;

use sha2::{Digest, Sha256};

use crate::error::ManifestError;

/// Largest accepted number of sample rows.
pub(crate) const MAX_ROWS: usize = 1_000_000;
/// Largest accepted number of traces of one sample.
pub(crate) const MAX_TRACES: usize = 64;

/// Columns the harness reads; other columns are allowed and ignored.
const SAMPLE_ID: &str = "sample_id";
const BATCH: &str = "batch";
const STRATUM: &str = "stratum";
const WEIGHT: &str = "weight";
const WEIGHT_CORPUS: &str = "weight_corpus";
const ANCHOR: &str = "anchor_candidate";
const TRACES: &str = "trace_paths";
/// Optional column: `1` for a retired sample.
pub(crate) const RETIRED: &str = "retired";

/// One parsed and verified panel manifest.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Manifest {
    /// `key=value` fields of the comment header.
    pub(crate) header: BTreeMap<String, String>,
    /// The verified SHA-256 of the body.
    pub(crate) body_sha256: String,
    /// Column names in order.
    pub(crate) columns: Vec<String>,
    /// Samples in manifest order.
    pub(crate) samples: Vec<Sample>,
}

/// One sample row of a panel manifest.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Sample {
    pub(crate) id: String,
    /// The batch the panel assigns the sample to.
    pub(crate) batch: String,
    pub(crate) stratum: String,
    /// Inverse sampling share within the panel's pool, as written.
    pub(crate) weight: Weight,
    /// Inverse sampling share scaled to the whole corpus, as written.
    pub(crate) weight_corpus: Weight,
    /// Whether the panel proposes the sample as an anchor.
    pub(crate) anchor: bool,
    /// Trace paths in manifest order.
    pub(crate) traces: Vec<PathBuf>,
    /// Whether the sample is retired: kept, but left out of estimates.
    pub(crate) retired: bool,
    /// The row's fields as written, in column order.
    pub(crate) fields: Vec<String>,
}

/// A weight with the text it was written as, so outputs echo it exactly.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Weight {
    pub(crate) value: f64,
    pub(crate) text: String,
}

/// Parses and verifies one manifest.
pub(crate) fn parse(bytes: &[u8]) -> Result<Manifest, ManifestError> {
    let text = std::str::from_utf8(bytes).map_err(|_| ManifestError::NotUtf8)?;
    let mut header = BTreeMap::new();
    let mut offset = 0;
    let mut line_number = 0;
    let columns = loop {
        let Some(line) = next_line(text, offset) else {
            return Err(ManifestError::MissingColumns);
        };
        line_number += 1;
        offset += line.len();
        let line = line.trim_end_matches(['\n', '\r']);
        if let Some(comment) = line.strip_prefix('#') {
            for (key, value) in comment
                .split_whitespace()
                .filter_map(|field| field.split_once('='))
            {
                header.insert(key.to_owned(), value.to_owned());
            }
        } else {
            break line.split('\t').collect::<Vec<_>>();
        }
    };
    let body = &text[offset..];
    let declared = header
        .get("body_sha256")
        .ok_or(ManifestError::MissingBodyHash)?;
    let actual = format!("{:x}", Sha256::digest(body.as_bytes()));
    if *declared != actual {
        return Err(ManifestError::BodyHash {
            declared: declared.clone(),
            actual,
        });
    }
    let index = |column: &'static str| {
        columns
            .iter()
            .position(|name| *name == column)
            .ok_or(ManifestError::MissingColumn { column })
    };
    let positions = [
        index(SAMPLE_ID)?,
        index(BATCH)?,
        index(STRATUM)?,
        index(WEIGHT)?,
        index(WEIGHT_CORPUS)?,
        index(ANCHOR)?,
        index(TRACES)?,
    ];
    let retired = columns.iter().position(|name| *name == RETIRED);

    let mut samples = Vec::new();
    let mut seen = BTreeSet::new();
    for line in body.lines() {
        line_number += 1;
        if samples.len() == MAX_ROWS {
            return Err(ManifestError::TooManyRows { maximum: MAX_ROWS });
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() != columns.len() {
            return Err(ManifestError::Fields {
                line: line_number,
                found: fields.len(),
                expected: columns.len(),
            });
        }
        let sample = row(&fields, positions, retired, line_number)?;
        if !seen.insert(sample.id.clone()) {
            return Err(ManifestError::DuplicateSample {
                line: line_number,
                sample: sample.id,
            });
        }
        samples.push(sample);
    }
    if let Some(declared) = header.get("rows")
        && declared.parse::<usize>().ok() != Some(samples.len())
    {
        return Err(ManifestError::RowCount {
            declared: declared.clone(),
            actual: samples.len(),
        });
    }
    Ok(Manifest {
        header,
        body_sha256: actual,
        columns: columns.iter().map(|column| (*column).to_owned()).collect(),
        samples,
    })
}

/// The line that starts at `offset`, with its line ending.
fn next_line(text: &str, offset: usize) -> Option<&str> {
    let rest = text.get(offset..).filter(|rest| !rest.is_empty())?;
    Some(rest.find('\n').map_or(rest, |end| &rest[..=end]))
}

fn row(
    fields: &[&str],
    positions: [usize; 7],
    retired: Option<usize>,
    line: usize,
) -> Result<Sample, ManifestError> {
    let [id, batch, stratum, weight, weight_corpus, anchor, traces] =
        positions.map(|index| fields[index]);
    if id.is_empty() {
        return Err(ManifestError::EmptySample { line });
    }
    let flag = |value: &str, column: &'static str| match value {
        "0" => Ok(false),
        "1" => Ok(true),
        _ => Err(ManifestError::Value {
            line,
            column,
            value: value.to_owned(),
        }),
    };
    let anchor = flag(anchor, ANCHOR)?;
    let retired = retired
        .map(|index| flag(fields[index], RETIRED))
        .transpose()?
        .unwrap_or(false);
    let traces = traces
        .split(';')
        .filter(|path| !path.is_empty())
        .map(PathBuf::from)
        .collect::<Vec<_>>();
    if traces.is_empty() || traces.len() > MAX_TRACES {
        return Err(ManifestError::Traces {
            line,
            found: traces.len(),
            maximum: MAX_TRACES,
        });
    }
    Ok(Sample {
        id: id.to_owned(),
        batch: batch.to_owned(),
        stratum: stratum.to_owned(),
        weight: parse_weight(weight, WEIGHT, line)?,
        weight_corpus: parse_weight(weight_corpus, WEIGHT_CORPUS, line)?,
        anchor,
        traces,
        retired,
        fields: fields.iter().map(|field| (*field).to_owned()).collect(),
    })
}

fn parse_weight(text: &str, column: &'static str, line: usize) -> Result<Weight, ManifestError> {
    match text.parse::<f64>() {
        Ok(value) if value.is_finite() && value >= 0.0 => Ok(Weight {
            value,
            text: text.to_owned(),
        }),
        _ => Err(ManifestError::Value {
            line,
            column,
            value: text.to_owned(),
        }),
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    const COLUMNS: &str =
        "sample_id\tbatch\tstratum\tweight\tweight_corpus\tanchor_candidate\ttrace_paths\n";

    /// A manifest with a correct body hash for the given rows.
    pub(crate) fn manifest(rows: &[&str]) -> String {
        let body = rows.iter().flat_map(|row| [*row, "\n"]).collect::<String>();
        format!(
            "# panel=dev version=v1 salt=s rows={}\n# truth_sha256=abc\n# body_sha256={:x} (of the rows)\n{COLUMNS}{body}",
            rows.len(),
            Sha256::digest(body.as_bytes())
        )
    }

    #[test]
    fn parses_a_verified_manifest() -> Result<(), ManifestError> {
        let text = manifest(&[
            "S2\tb1\tclean\t2.5\t3\t1\t/t/a.ab1;/t/b.ab1",
            "S1\tb1\tlegacy_addon:clean\t0.0\t0.0\t0\t/t/c.ab1;",
        ]);
        let parsed = parse(text.as_bytes())?;
        assert_eq!(parsed.header["panel"], "dev");
        assert_eq!(parsed.header["truth_sha256"], "abc");
        assert_eq!(parsed.samples.len(), 2);
        let first = &parsed.samples[0];
        assert_eq!(first.id, "S2");
        assert_eq!(first.weight.value, 2.5);
        assert_eq!(first.weight_corpus.text, "3");
        assert!(first.anchor);
        assert_eq!(
            first.traces,
            [PathBuf::from("/t/a.ab1"), PathBuf::from("/t/b.ab1")]
        );
        assert_eq!(parsed.samples[1].traces.len(), 1);
        assert!(!parsed.samples[1].anchor);
        assert!(!parsed.samples[0].retired);
        assert_eq!(parsed.samples[0].fields[0], "S2");
        assert_eq!(parsed.columns.len(), 7);
        Ok(())
    }

    #[test]
    fn reads_the_optional_retired_column() -> Result<(), ManifestError> {
        let columns = "sample_id\tbatch\tstratum\tweight\tweight_corpus\tanchor_candidate\ttrace_paths\tretired\n";
        let body = "S1\tb\tclean\t0\t0\t0\t/t/a.ab1\t1\nS2\tb\tclean\t2\t2\t0\t/t/b.ab1\t0\n";
        let text = format!(
            "# body_sha256={:x}\n{columns}{body}",
            Sha256::digest(body.as_bytes())
        );
        let parsed = parse(text.as_bytes())?;
        assert!(parsed.samples[0].retired);
        assert!(!parsed.samples[1].retired);
        let body = "S1\tb\tclean\t0\t0\t0\t/t/a.ab1\tyes\n";
        let text = format!(
            "# body_sha256={:x}\n{columns}{body}",
            Sha256::digest(body.as_bytes())
        );
        assert_eq!(
            parse(text.as_bytes()),
            Err(ManifestError::Value {
                line: 3,
                column: RETIRED,
                value: "yes".into(),
            })
        );
        Ok(())
    }

    #[test]
    fn rejects_a_body_that_does_not_match_its_hash() {
        let text = manifest(&["S1\tb\tclean\t1\t1\t0\t/t/a.ab1"]).replace("/t/a", "/t/x");
        assert!(matches!(
            parse(text.as_bytes()),
            Err(ManifestError::BodyHash { .. })
        ));
    }

    #[test]
    fn rejects_malformed_manifests() {
        let row = "S1\tb\tclean\t1\t1\t0\t/t/a.ab1";
        assert_eq!(
            parse(b"# body_sha256=x\n"),
            Err(ManifestError::MissingColumns)
        );
        assert_eq!(
            parse(format!("{COLUMNS}{row}\n").as_bytes()),
            Err(ManifestError::MissingBodyHash)
        );
        let cases = [
            (
                manifest(&[row, row]),
                ManifestError::DuplicateSample {
                    line: 6,
                    sample: "S1".into(),
                },
            ),
            (
                manifest(&["S1\tb\tclean\t1\t1\t0"]),
                ManifestError::Fields {
                    line: 5,
                    found: 6,
                    expected: 7,
                },
            ),
            (
                manifest(&["S1\tb\tclean\t-1\t1\t0\t/t/a.ab1"]),
                ManifestError::Value {
                    line: 5,
                    column: "weight",
                    value: "-1".into(),
                },
            ),
            (
                manifest(&["S1\tb\tclean\t1\tNaN\t0\t/t/a.ab1"]),
                ManifestError::Value {
                    line: 5,
                    column: "weight_corpus",
                    value: "NaN".into(),
                },
            ),
            (
                manifest(&["S1\tb\tclean\t1\t1\tyes\t/t/a.ab1"]),
                ManifestError::Value {
                    line: 5,
                    column: "anchor_candidate",
                    value: "yes".into(),
                },
            ),
            (
                manifest(&["S1\tb\tclean\t1\t1\t0\t;"]),
                ManifestError::Traces {
                    line: 5,
                    found: 0,
                    maximum: MAX_TRACES,
                },
            ),
            (
                manifest(&["\tb\tclean\t1\t1\t0\t/t/a.ab1"]),
                ManifestError::EmptySample { line: 5 },
            ),
            (
                manifest(&[row]).replace("rows=1", "rows=2"),
                ManifestError::RowCount {
                    declared: "2".into(),
                    actual: 1,
                },
            ),
            (
                manifest(&[row]).replace("anchor_candidate", "anchor"),
                ManifestError::MissingColumn {
                    column: "anchor_candidate",
                },
            ),
        ];
        for (text, expected) in cases {
            assert_eq!(parse(text.as_bytes()), Err(expected));
        }
        assert_eq!(parse(&[0xff, b'\n']), Err(ManifestError::NotUtf8));
    }
}
