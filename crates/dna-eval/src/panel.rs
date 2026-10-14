//! A new panel version from an earlier one and the truth database.
//!
//! Membership is frozen: the new manifest keeps every row of the parent, in
//! order. A sample the database has no record of is retired (`retired = 1`,
//! both weights `0`), and the weights of the other samples of its stratum are
//! scaled in proportion, `w' = w · total / Σ staying w`, so each stratum still
//! sums to its population and keeps its differences within; a zero weight
//! stays zero. When every staying weight is zero the retired weight cannot be
//! carried, the weights stay zero, and the stratum is reported. Strata without
//! a newly retired sample keep their weights as written. Each row
//! records its sample's truth digest (`truth_digest`), and the header the
//! panel's, so a run can tell whether the truth of the panel changed; the
//! database path, file SHA-256, and creation time are recorded as provenance
//! only.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::Path;

use sha2::{Digest, Sha256};

use crate::manifest::{Manifest, RETIRED, Sample, TRUTH_DIGEST, encode};
use crate::truth::{self, TruthDatabase};

/// The text of the new manifest and what changed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Derived {
    pub(crate) text: String,
    pub(crate) body_sha256: String,
    /// The panel's truth digest.
    pub(crate) truth_digest: String,
    pub(crate) rows: usize,
    /// Samples retired in this version.
    pub(crate) retired: usize,
    /// Samples whose database batch differs from the panel's.
    pub(crate) batch_differs: usize,
    /// Strata that lost weight because every staying sample weighs zero.
    pub(crate) stranded: Vec<String>,
}

/// One stratum's weight totals, over all samples and over the staying ones,
/// per weight column, and whether a sample in it is newly retired.
#[derive(Default)]
struct Stratum {
    total: [f64; 2],
    staying: [f64; 2],
    changed: bool,
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
    let mut column = |name: &str| {
        columns
            .iter()
            .position(|column| column == name)
            .unwrap_or_else(|| {
                columns.push(name.to_owned());
                columns.len() - 1
            })
    };
    let retired_column = column(RETIRED);
    let digest_column = column(TRUTH_DIGEST);
    let weight_columns =
        ["weight", "weight_corpus"].map(|name| columns.iter().position(|column| column == name));

    let mut strata = BTreeMap::<&str, Stratum>::new();
    for sample in &parent.samples {
        let stratum = strata.entry(sample.stratum.as_str()).or_default();
        let weights = [sample.weight.value, sample.weight_corpus.value];
        for (total, weight) in stratum.total.iter_mut().zip(weights) {
            *total += weight;
        }
        if retired(sample) {
            stratum.changed |= !sample.retired;
        } else {
            for (staying, weight) in stratum.staying.iter_mut().zip(weights) {
                *staying += weight;
            }
        }
    }

    let digests = parent
        .samples
        .iter()
        .map(|sample| truth::digest(truth.database.samples.get(&sample.id)))
        .collect::<Vec<_>>();
    let mut body = String::new();
    for (sample, digest) in parent.samples.iter().zip(&digests) {
        let mut fields = sample.fields.clone();
        fields.resize(columns.len(), String::new());
        let gone = retired(sample);
        fields[retired_column] = u8::from(gone).to_string();
        fields[digest_column].clone_from(digest);
        if let Some(stratum) = strata.get(sample.stratum.as_str())
            && stratum.changed
        {
            let weights = [sample.weight.value, sample.weight_corpus.value];
            for (column, index) in weight_columns.iter().enumerate() {
                if let Some(index) = index {
                    let staying = stratum.staying[column];
                    let weight = if gone || staying <= 0.0 {
                        0.0
                    } else {
                        weights[column] * stratum.total[column] / staying
                    };
                    fields[*index] = if weight == 0.0 {
                        "0".to_owned()
                    } else {
                        format!("{weight:.4}")
                    };
                }
            }
        }
        body.push_str(&fields.join("\t"));
        body.push('\n');
    }
    let body_sha256 = format!("{:x}", Sha256::digest(body.as_bytes()));
    let truth_digest = truth::panel_digest(
        parent
            .samples
            .iter()
            .zip(&digests)
            .map(|(sample, digest)| (sample.id.as_str(), digest.as_str())),
    );

    let field = |key: &str| encode(parent.header.get(key).map_or("-", String::as_str));
    let info = |key: &str| encode(truth.database.info.get(key).map_or("-", String::as_str));
    let (version, created) = (encode(version), encode(created));
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
        "# truth={} truth_digest={truth_digest} truth_database_sha256={} truth_created={} \
         truth_samples={}",
        encode(&truth.path.display().to_string()),
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
        truth_digest,
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
        stranded: strata
            .iter()
            .filter(|(_, stratum)| {
                stratum.changed
                    && (0..2)
                        .any(|column| stratum.total[column] > 0.0 && stratum.staying[column] <= 0.0)
            })
            .map(|(name, _)| (*name).to_owned())
            .collect(),
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
        let mut sql = String::from(truth::tests::FIXTURE);
        sql.push_str(
            "DELETE FROM samples; INSERT INTO samples VALUES ('S1', 'b1', 0, 0, 0, 0, 0),
             ('S3', 'b1', 0, 0, 0, 0, 0), ('S4', 'b9', 0, 0, 0, 0, 0);",
        );
        let path = truth::tests::database(directory.path(), &sql)?;
        let database = truth::load(&path, &BTreeSet::from(["S1", "S2", "S3", "S4"]))?;
        let source = Source {
            path: Path::new("/data/lab stats=v1/statistics.sqlite"),
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
        assert_eq!(child.header["truth_database_sha256"], "abc");
        // The path with a space and '=' reads back.
        assert_eq!(
            child.header["truth"],
            "/data/lab stats=v1/statistics.sqlite"
        );
        assert_eq!(child.header["truth_digest"], derived.truth_digest);
        assert_eq!(
            child.samples[1].truth_digest.as_deref(),
            Some(truth::digest(None).as_str())
        );
        assert_eq!(
            child.samples[0].truth_digest.as_deref(),
            Some(truth::digest(database.samples.get("S1")).as_str())
        );
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

        assert!(derived.stranded.is_empty());

        // Deriving again keeps the retired sample retired and the weights.
        let again = derive(&child, &source, "v3", "2026-10-16");
        assert_eq!(again.retired, 0);
        let grandchild = manifest::parse(again.text.as_bytes())?;
        assert_eq!(grandchild.samples[0].weight.text, "4.5000");
        assert!(grandchild.samples[1].retired);
        assert_eq!(grandchild.columns, child.columns);
        assert_eq!(again.truth_digest, derived.truth_digest);
        Ok(())
    }

    #[test]
    fn scales_staying_weights_in_proportion() -> Result<(), Box<dyn std::error::Error>> {
        // Weights 3, 3, 0 in one stratum and 0 in another; the second 3 and
        // the lone zero-weight sample have no truth.
        let parent = manifest::parse(
            manifest::tests::manifest(&[
                "S1\tb1\tclean\t3\t3\t0\t/t/a.ab1",
                "S2\tb1\tclean\t3\t3\t0\t/t/b.ab1",
                "S3\tb1\tclean\t0\t0\t0\t/t/c.ab1",
                "S4\tb1\thard\t2\t2\t0\t/t/d.ab1",
                "S5\tb1\thard\t0\t0\t0\t/t/e.ab1",
            ])
            .as_bytes(),
        )?;
        let directory = tempfile::tempdir()?;
        let mut sql = String::from(truth::tests::FIXTURE);
        sql.push_str(
            "DELETE FROM samples; INSERT INTO samples VALUES ('S1', 'b1', 0, 0, 0, 0, 0),
             ('S3', 'b1', 0, 0, 0, 0, 0), ('S5', 'b1', 0, 0, 0, 0, 0);",
        );
        let path = truth::tests::database(directory.path(), &sql)?;
        let database = truth::load(&path, &BTreeSet::from(["S1", "S2", "S3", "S4", "S5"]))?;
        let source = Source {
            path: Path::new("/data/statistics.sqlite"),
            sha256: "abc",
            database: &database,
        };
        let derived = derive(&parent, &source, "v2", "2026-10-15");
        let child = manifest::parse(derived.text.as_bytes())?;
        let weights = child
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
        // The clean total 6 moves to S1; S3 stays zero. In hard nothing stays
        // with weight, so the weight is lost and the stratum reported.
        assert_eq!(
            weights,
            [
                ("S1", false, "6.0000"),
                ("S2", true, "0"),
                ("S3", false, "0"),
                ("S4", true, "0"),
                ("S5", false, "0"),
            ]
        );
        assert_eq!(derived.stranded, ["hard"]);
        Ok(())
    }
}
