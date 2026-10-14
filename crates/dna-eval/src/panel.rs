//! A new panel version from an earlier one and the truth database.
//!
//! Membership is frozen: the new manifest keeps every row of the parent, in
//! order. A sample the database has no record of is retired (`retired = 1`,
//! both weights `0`), and the other samples of its stratum share the
//! stratum's weight, so each stratum still sums to its population. Strata
//! without a newly retired sample keep their weights as written. The header
//! names the database and its SHA-256 as the truth.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::Path;

use sha2::{Digest, Sha256};

use crate::manifest::{Manifest, RETIRED, Sample};
use crate::truth::TruthDatabase;

/// The text of the new manifest and what changed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Derived {
    pub(crate) text: String,
    pub(crate) body_sha256: String,
    pub(crate) rows: usize,
    /// Samples retired in this version.
    pub(crate) retired: usize,
    /// Samples whose database batch differs from the panel's.
    pub(crate) batch_differs: usize,
}

/// The truth a new version is built against.
pub(crate) struct Source<'a> {
    pub(crate) path: &'a Path,
    pub(crate) sha256: &'a str,
    pub(crate) database: &'a TruthDatabase,
}

/// Derives version `version` of `parent`, created on `created`.
pub(crate) fn derive(
    parent: &Manifest,
    truth: &Source<'_>,
    version: &str,
    created: &str,
) -> Derived {
    let missing = |sample: &Sample| !truth.database.samples.contains_key(&sample.id);
    let retired = |sample: &Sample| sample.retired || missing(sample);
    let newly = parent
        .samples
        .iter()
        .filter(|sample| retired(sample) && !sample.retired)
        .count();

    let mut columns = parent.columns.clone();
    let retired_column = columns
        .iter()
        .position(|column| column == RETIRED)
        .unwrap_or_else(|| {
            columns.push(RETIRED.to_owned());
            columns.len() - 1
        });
    let weight_columns =
        ["weight", "weight_corpus"].map(|name| columns.iter().position(|column| column == name));

    // Per stratum: the weight total, the samples that stay, and whether a
    // sample is newly retired.
    let mut strata = BTreeMap::<&str, ([f64; 2], usize, bool)>::new();
    for sample in &parent.samples {
        let entry = strata.entry(sample.stratum.as_str()).or_default();
        entry.0[0] += sample.weight.value;
        entry.0[1] += sample.weight_corpus.value;
        if retired(sample) {
            entry.2 |= !sample.retired;
        } else {
            entry.1 += 1;
        }
    }

    let mut body = String::new();
    for sample in &parent.samples {
        let mut fields = sample.fields.clone();
        fields.resize(columns.len(), String::new());
        let gone = retired(sample);
        fields[retired_column] = u8::from(gone).to_string();
        if let Some((totals, staying, changed)) = strata.get(sample.stratum.as_str())
            && *changed
        {
            for (index, total) in weight_columns.iter().zip(totals) {
                if let Some(index) = index {
                    fields[*index] = if gone || *staying == 0 {
                        "0".to_owned()
                    } else {
                        format!("{:.4}", total / *staying as f64)
                    };
                }
            }
        }
        body.push_str(&fields.join("\t"));
        body.push('\n');
    }
    let body_sha256 = format!("{:x}", Sha256::digest(body.as_bytes()));

    let field = |key: &str| parent.header.get(key).map_or("-", String::as_str);
    let info = |key: &str| truth.database.info.get(key).map_or("-", String::as_str);
    let total_retired = parent
        .samples
        .iter()
        .filter(|sample| retired(sample))
        .count();
    let mut text = String::new();
    let _ = writeln!(
        text,
        "# panel={} version={version} salt={} created={created} rows={} retired={total_retired} \
         parent_version={} parent_body_sha256={}",
        field("panel"),
        field("salt"),
        parent.samples.len(),
        field("version"),
        parent.body_sha256,
    );
    let _ = writeln!(
        text,
        "# truth={} truth_sha256={} truth_created={} truth_samples={}",
        truth.path.display(),
        truth.sha256,
        info("created_at"),
        info("samples"),
    );
    let _ = writeln!(
        text,
        "# body_sha256={body_sha256} (sha256 of the data rows below the column header)"
    );
    text.push_str(&columns.join("\t"));
    text.push('\n');
    text.push_str(&body);

    Derived {
        text,
        body_sha256,
        rows: parent.samples.len(),
        retired: newly,
        batch_differs: parent
            .samples
            .iter()
            .filter(|sample| {
                truth
                    .database
                    .samples
                    .get(&sample.id)
                    .is_some_and(|record| record.batch != sample.batch)
            })
            .count(),
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use super::*;
    use crate::manifest;

    #[test]
    fn retires_samples_without_truth_and_keeps_strata_whole()
    -> Result<(), Box<dyn std::error::Error>> {
        let parent = manifest::parse(
            manifest::tests::manifest(&[
                "S1\tb1\tclean\t3\t6\t1\t/t/a.ab1",
                "S2\tb1\tclean\t3\t6\t0\t/t/b.ab1",
                "S3\tb1\tclean\t3\t6\t0\t/t/c.ab1",
                "S4\tb2\thard\t2.5\t2.5\t0\t/t/d.ab1",
            ])
            .as_bytes(),
        )?;
        // S2 is not in the database; S4's batch differs.
        let directory = tempfile::tempdir()?;
        let mut sql = String::from(crate::truth::tests::FIXTURE);
        sql.push_str(
            "DELETE FROM samples; INSERT INTO samples VALUES ('S1', 'b1', 0, 0, 0, 0, 0),
             ('S3', 'b1', 0, 0, 0, 0, 0), ('S4', 'b9', 0, 0, 0, 0, 0);",
        );
        let path = crate::truth::tests::database(directory.path(), &sql)?;
        let database = crate::truth::load(&path, &BTreeSet::from(["S1", "S2", "S3", "S4"]))?;
        let source = Source {
            path: Path::new("/data/statistics.sqlite"),
            sha256: "abc",
            database: &database,
        };
        let derived = derive(&parent, &source, "v2", "2026-10-15");
        assert_eq!(
            (derived.rows, derived.retired, derived.batch_differs),
            (4, 1, 1)
        );

        let child = manifest::parse(derived.text.as_bytes())?;
        assert_eq!(child.body_sha256, derived.body_sha256);
        assert_eq!(child.header["version"], "v2");
        assert_eq!(child.header["parent_version"], "v1");
        assert_eq!(child.header["parent_body_sha256"], parent.body_sha256);
        assert_eq!(child.header["truth_sha256"], "abc");
        assert_eq!(child.header["truth_samples"], "3");
        assert_eq!(child.header["retired"], "1");
        let ids = child
            .samples
            .iter()
            .map(|sample| {
                (
                    sample.id.as_str(),
                    sample.retired,
                    sample.weight.text.as_str(),
                )
            })
            .collect::<Vec<_>>();
        assert_eq!(
            ids,
            [
                ("S1", false, "4.5000"),
                ("S2", true, "0"),
                ("S3", false, "4.5000"),
                ("S4", false, "2.5"),
            ]
        );
        assert_eq!(child.samples[0].weight_corpus.text, "9.0000");
        assert_eq!(child.samples[0].traces, parent.samples[0].traces);

        // Deriving again keeps the retired sample retired and the weights.
        let again = derive(&child, &source, "v3", "2026-10-16");
        assert_eq!(again.retired, 0);
        let grandchild = manifest::parse(again.text.as_bytes())?;
        assert_eq!(grandchild.samples[0].weight.text, "4.5000");
        assert!(grandchild.samples[1].retired);
        assert_eq!(grandchild.columns, child.columns);
        Ok(())
    }
}
