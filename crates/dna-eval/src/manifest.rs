//! Panel manifests: the samples of one test panel, their stratum, weights,
//! and traces.
//!
//! A manifest starts with `#` comment lines of `key=value` fields (panel,
//! version, salt, truth and body hashes; in a value `%`, `=`, and white space
//! are written as `%XX` escapes of their UTF-8 bytes), then one tab-separated
//! column
//! header, then one row per sample. `body_sha256` is the SHA-256 of every byte
//! after the column header line and must match. From version 2 a `retired`
//! column marks samples kept for continuity but left out of estimates; a
//! manifest without it retires nothing. From version 3 a `truth_digest`
//! column records each sample's truth digest when the panel was built, and
//! the header's `truth_digest` the panel's.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::path::PathBuf;
use std::str::FromStr;

use serde::{Deserialize, Serialize};
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
/// Optional column: the sample's truth digest when the panel was built.
pub(crate) const TRUTH_DIGEST: &str = "truth_digest";
/// Optional column: the batch's run date, `YYYY-MM-DD`, or empty when
/// unknown.
pub(crate) const BATCH_DATE: &str = "batch_date";

/// A calendar month, written `YYYY-MM`: the unit of the current data
/// (`--current-since`) and of a batch's run date. Parsed once, here.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(try_from = "String", into = "String")]
pub(crate) struct Month {
    year: u16,
    month: u8,
}

impl Month {
    /// The month of a date written `YYYY-MM-DD`, with a day of 01 to 31.
    pub(crate) fn of_date(date: &str) -> Option<Self> {
        let (month, day) = date.split_at_checked(7)?;
        let day = day.strip_prefix('-')?;
        let valid_day = day.len() == 2
            && day.bytes().all(|byte| byte.is_ascii_digit())
            && (1..=31).contains(&day.parse::<u8>().ok()?);
        valid_day.then(|| month.parse().ok()).flatten()
    }
}

impl FromStr for Month {
    type Err = String;

    fn from_str(text: &str) -> Result<Self, Self::Err> {
        let invalid = || format!("{text:?} is not a month written YYYY-MM");
        let (year, month) = text.split_once('-').ok_or_else(invalid)?;
        let digits = |part: &str, length: usize| {
            part.len() == length && part.bytes().all(|byte| byte.is_ascii_digit())
        };
        if !digits(year, 4) || !digits(month, 2) {
            return Err(invalid());
        }
        let month = month.parse::<u8>().map_err(|_| invalid())?;
        if !(1..=12).contains(&month) {
            return Err(invalid());
        }
        Ok(Self {
            year: year.parse().map_err(|_| invalid())?,
            month,
        })
    }
}

impl TryFrom<String> for Month {
    type Error = String;

    fn try_from(text: String) -> Result<Self, Self::Error> {
        text.parse()
    }
}

impl From<Month> for String {
    fn from(month: Month) -> Self {
        month.to_string()
    }
}

impl std::fmt::Display for Month {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(formatter, "{:04}-{:02}", self.year, self.month)
    }
}

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
    /// The sample's truth digest when the panel was built, if recorded.
    pub(crate) truth_digest: Option<String>,
    /// The batch's run date as written, when the manifest has the column.
    pub(crate) batch_date: Option<String>,
    /// The month of the batch's run date, when it has one.
    pub(crate) batch_month: Option<Month>,
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
                let value = decode(value).ok_or_else(|| ManifestError::HeaderValue {
                    line: line_number,
                    value: value.to_owned(),
                })?;
                header.insert(key.to_owned(), value);
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
    let optional = [
        columns.iter().position(|name| *name == RETIRED),
        columns.iter().position(|name| *name == TRUTH_DIGEST),
        columns.iter().position(|name| *name == BATCH_DATE),
    ];

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
        let sample = row(&fields, positions, optional, line_number)?;
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

/// A header value as `parse` reads it back: `%`, `=`, and white space become
/// `%XX` escapes of their UTF-8 bytes.
pub(crate) fn encode(value: &str) -> String {
    let mut encoded = String::with_capacity(value.len());
    for symbol in value.chars() {
        if symbol == '%' || symbol == '=' || symbol.is_whitespace() {
            let mut bytes = [0; 4];
            for byte in symbol.encode_utf8(&mut bytes).bytes() {
                let _ = write!(encoded, "%{byte:02X}");
            }
        } else {
            encoded.push(symbol);
        }
    }
    encoded
}

/// The value of an encoded header value, or `None` for a malformed escape.
fn decode(value: &str) -> Option<String> {
    let mut bytes = Vec::with_capacity(value.len());
    let mut rest = value.as_bytes();
    while let Some((&byte, tail)) = rest.split_first() {
        if byte == b'%' {
            let hex = std::str::from_utf8(tail.get(..2)?).ok()?;
            bytes.push(u8::from_str_radix(hex, 16).ok()?);
            rest = tail.get(2..)?;
        } else {
            bytes.push(byte);
            rest = tail;
        }
    }
    String::from_utf8(bytes).ok()
}

/// The line that starts at `offset`, with its line ending.
fn next_line(text: &str, offset: usize) -> Option<&str> {
    let rest = text.get(offset..).filter(|rest| !rest.is_empty())?;
    Some(rest.find('\n').map_or(rest, |end| &rest[..=end]))
}

fn row(
    fields: &[&str],
    positions: [usize; 7],
    [retired, truth_digest, batch_date]: [Option<usize>; 3],
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
    let truth_digest = truth_digest
        .map(|index| {
            let value = fields[index];
            if value.len() == 64
                && value
                    .bytes()
                    .all(|byte| matches!(byte, b'0'..=b'9' | b'a'..=b'f'))
            {
                Ok(value.to_owned())
            } else {
                Err(ManifestError::Value {
                    line,
                    column: TRUTH_DIGEST,
                    value: value.to_owned(),
                })
            }
        })
        .transpose()?;
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
        truth_digest,
        batch_month: batch_date
            .map(|index| fields[index])
            .filter(|date| !date.is_empty())
            .map(|date| {
                Month::of_date(date).ok_or_else(|| ManifestError::Value {
                    line,
                    column: BATCH_DATE,
                    value: date.to_owned(),
                })
            })
            .transpose()?,
        batch_date: batch_date.map(|index| fields[index].to_owned()),
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
    fn reads_the_optional_batch_date_column() -> Result<(), ManifestError> {
        let columns = "sample_id\tbatch\tbatch_date\tstratum\tweight\tweight_corpus\tanchor_candidate\ttrace_paths\n";
        let body = "S1\tb\t2025-05-25\tclean\t0\t0\t0\t/t/a.ab1\n";
        let text = format!(
            "# body_sha256={:x}\n{columns}{body}",
            Sha256::digest(body.as_bytes())
        );
        let parsed = parse(text.as_bytes())?;
        assert_eq!(parsed.samples[0].batch_date.as_deref(), Some("2025-05-25"));
        assert_eq!(parsed.samples[0].batch_month, "2025-05".parse().ok());
        // An empty date is unknown; a malformed one stops the parse.
        for (date, month) in [("", None), ("2026-01-31", "2026-01".parse().ok())] {
            let body = format!("S1\tb\t{date}\tclean\t0\t0\t0\t/t/a.ab1\n");
            let text = format!(
                "# body_sha256={:x}\n{columns}{body}",
                Sha256::digest(body.as_bytes())
            );
            assert_eq!(parse(text.as_bytes())?.samples[0].batch_month, month);
        }
        for date in [
            "2026-1-01",
            "2026-13-01",
            "2026-01",
            "2026-01-00",
            "2026/01/01",
            "20260101",
        ] {
            let body = format!("S1\tb\t{date}\tclean\t0\t0\t0\t/t/a.ab1\n");
            let text = format!(
                "# body_sha256={:x}\n{columns}{body}",
                Sha256::digest(body.as_bytes())
            );
            assert!(
                matches!(
                    parse(text.as_bytes()),
                    Err(ManifestError::Value {
                        column: BATCH_DATE,
                        ..
                    })
                ),
                "{date}"
            );
        }
        Ok(())
    }

    #[test]
    fn months_are_year_dash_month() {
        let month = |text: &str| text.parse::<Month>().ok();
        assert!(month("2026-01").is_some() && month("2025-12").is_some());
        for text in [
            "2026-1",
            "2026-13",
            "2026/01",
            "2026-00",
            "26-01",
            "2026-01-01",
            "",
        ] {
            assert!(month(text).is_none(), "{text}");
        }
        assert!(month("2025-12") < month("2026-01"));
        assert_eq!(
            month("2026-01").map(|month| month.to_string()).as_deref(),
            Some("2026-01")
        );
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
        assert_eq!(parsed.samples[0].truth_digest, None);
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

    #[test]
    fn header_values_round_trip_through_escapes() -> Result<(), ManifestError> {
        let value = "/data/my stats=v1/100%\tdone é";
        let encoded = encode(value);
        assert!(!encoded.contains([' ', '=', '\t']));
        let body = "S1\tb\tclean\t1\t1\t0\t/t/a.ab1\n";
        let text = format!(
            "# truth={encoded} body_sha256={:x}\n{COLUMNS}{body}",
            Sha256::digest(body.as_bytes())
        );
        assert_eq!(parse(text.as_bytes())?.header["truth"], value);
        let text = format!(
            "# truth=100%G body_sha256={:x}\n{COLUMNS}{body}",
            Sha256::digest(body.as_bytes())
        );
        assert_eq!(
            parse(text.as_bytes()),
            Err(ManifestError::HeaderValue {
                line: 1,
                value: "100%G".into(),
            })
        );
        Ok(())
    }

    #[test]
    fn reads_the_optional_truth_digest_column() -> Result<(), ManifestError> {
        let columns = "sample_id\tbatch\tstratum\tweight\tweight_corpus\tanchor_candidate\ttrace_paths\ttruth_digest\n";
        let digest = "a".repeat(64);
        let body = format!("S1\tb\tclean\t1\t1\t0\t/t/a.ab1\t{digest}\n");
        let text = format!(
            "# body_sha256={:x}\n{columns}{body}",
            Sha256::digest(body.as_bytes())
        );
        assert_eq!(
            parse(text.as_bytes())?.samples[0].truth_digest,
            Some(digest)
        );
        let body = "S1\tb\tclean\t1\t1\t0\t/t/a.ab1\tABC\n";
        let text = format!(
            "# body_sha256={:x}\n{columns}{body}",
            Sha256::digest(body.as_bytes())
        );
        assert_eq!(
            parse(text.as_bytes()),
            Err(ManifestError::Value {
                line: 3,
                column: TRUTH_DIGEST,
                value: "ABC".into(),
            })
        );
        Ok(())
    }
}
