//! Reviewer truth from the lab's statistics database (`statistics.sqlite`).
//!
//! The database is opened read-only. Its `info` table must declare
//! `schema_version` 1. Per sample, the harness reads:
//!
//! - `samples (sample_id, batch, …)`: the batch of the reviewed result;
//! - `variants (sample_id, kind, pos, base, ref, seq)`: one row per reviewed
//!   allele, with `-` for an inserted or deleted base;
//! - `intervals (sample_id, region, start_pos, end_pos)`: the analysed
//!   intervals, closed and 1-based;
//! - `sample_status (batch, sample_id, status, …)`: the sample's pipeline
//!   status in its batch (`OK` when its reviewed result is complete).
//!
//! Only the panel's samples are read.
//!
//! A sample's truth digest is the SHA-256 of a canonical text of exactly what
//! the harness uses: the batch, the status, the alleles, and the analysed
//! intervals, each sorted. A panel's truth digest is the SHA-256 of its
//! samples' digests in panel order, so a database refresh that changes no
//! panel sample keeps it.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::path::Path;

use rusqlite::{Connection, OpenFlags, OptionalExtension, params};
use sha2::{Digest, Sha256};

use crate::compare::Allele;
use crate::error::TruthError;
use crate::manifest::Manifest;

/// Largest accepted database file.
pub(crate) const MAX_TRUTH_BYTES: u64 = 4 * 1024 * 1024 * 1024;
/// Largest accepted number of alleles of one sample.
const MAX_ALLELES: usize = 10_000;
/// The database schema this module reads.
const SCHEMA_VERSION: &str = "1";

/// One sample's reviewer truth.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct TruthSample {
    pub(crate) batch: String,
    /// Pipeline status in that batch, when recorded.
    pub(crate) status: Option<String>,
    pub(crate) alleles: Vec<Allele>,
    /// Analysed intervals over every region, closed and 1-based.
    pub(crate) intervals: Vec<(i64, i64)>,
}

/// The panel's records and the database's `info` fields.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct TruthDatabase {
    pub(crate) info: BTreeMap<String, String>,
    pub(crate) samples: BTreeMap<String, TruthSample>,
}

/// The truth digest of one sample; a sample without a record has the digest
/// of `missing`.
pub(crate) fn digest(record: Option<&TruthSample>) -> String {
    let mut text = String::new();
    match record {
        None => text.push_str("missing\n"),
        Some(record) => {
            let mut alleles = record.alleles.iter().collect::<Vec<_>>();
            alleles.sort();
            let mut intervals = record.intervals.clone();
            intervals.sort_unstable();
            let _ = write!(
                text,
                "batch\t{}\nstatus\t{}\n",
                record.batch,
                record.status.as_deref().unwrap_or("")
            );
            for allele in alleles {
                let _ = writeln!(
                    text,
                    "allele\t{}\t{}\t{}",
                    allele.position, allele.reference, allele.alternate
                );
            }
            for (start, end) in intervals {
                let _ = writeln!(text, "interval\t{start}\t{end}");
            }
        }
    }
    format!("{:x}", Sha256::digest(text.as_bytes()))
}

/// The truth digest of a panel from its samples' identifiers and digests,
/// in panel order.
pub(crate) fn panel_digest<'a>(samples: impl IntoIterator<Item = (&'a str, &'a str)>) -> String {
    let mut hasher = Sha256::new();
    for (sample, digest) in samples {
        hasher.update(format!("{sample}\t{digest}\n").as_bytes());
    }
    format!("{:x}", hasher.finalize())
}

/// How the panel's truth compares with the digests its manifest records.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct TruthCheck {
    /// The panel digest the manifest header records, if any.
    pub(crate) recorded: Option<String>,
    /// The panel digest of the current database.
    pub(crate) current: String,
    /// Samples whose recorded digest differs from the current one, in panel
    /// order.
    pub(crate) changed: Vec<ChangedSample>,
}

/// One sample whose truth changed since the panel was built.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ChangedSample {
    pub(crate) sample: String,
    pub(crate) stratum: String,
    pub(crate) retired: bool,
    pub(crate) recorded: String,
    pub(crate) current: String,
}

/// Whether the truth of the panel is the one its manifest records.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum TruthState {
    /// The panel digest equals the recorded one.
    Unchanged,
    /// The panel digest differs: a new panel version is due.
    Changed,
    /// The manifest records no truth digest (panel versions 1 and 2).
    NotRecorded,
}

impl TruthState {
    pub(crate) const fn label(self) -> &'static str {
        match self {
            Self::Unchanged => "unchanged",
            Self::Changed => "changed",
            Self::NotRecorded => "not_recorded",
        }
    }
}

impl TruthCheck {
    /// Compares the database's truth of the manifest's samples with the
    /// digests the manifest records.
    pub(crate) fn new(manifest: &Manifest, database: &TruthDatabase) -> Self {
        let digests = manifest
            .samples
            .iter()
            .map(|sample| digest(database.samples.get(&sample.id)))
            .collect::<Vec<_>>();
        let changed = manifest
            .samples
            .iter()
            .zip(&digests)
            .filter_map(|(sample, current)| {
                let recorded = sample.truth_digest.as_ref()?;
                (recorded != current).then(|| ChangedSample {
                    sample: sample.id.clone(),
                    stratum: sample.stratum.clone(),
                    retired: sample.retired,
                    recorded: recorded.clone(),
                    current: current.clone(),
                })
            })
            .collect();
        Self {
            recorded: manifest.header.get("truth_digest").cloned(),
            current: panel_digest(
                manifest
                    .samples
                    .iter()
                    .zip(&digests)
                    .map(|(sample, digest)| (sample.id.as_str(), digest.as_str())),
            ),
            changed,
        }
    }

    pub(crate) fn state(&self) -> TruthState {
        match &self.recorded {
            None => TruthState::NotRecorded,
            Some(recorded) if *recorded == self.current && self.changed.is_empty() => {
                TruthState::Unchanged
            }
            Some(_) => TruthState::Changed,
        }
    }

    /// `truth_changes.tsv`: the samples whose truth changed.
    pub(crate) fn changes(&self) -> String {
        let mut text =
            String::from("sample_id\tstratum\tretired\trecorded_digest\tcurrent_digest\n");
        for change in &self.changed {
            let _ = writeln!(
                text,
                "{}\t{}\t{}\t{}\t{}",
                change.sample,
                change.stratum,
                u8::from(change.retired),
                change.recorded,
                change.current
            );
        }
        text
    }
}

/// Reads the records of the `wanted` samples from the database at `path`.
pub(crate) fn load(path: &Path, wanted: &BTreeSet<&str>) -> Result<TruthDatabase, TruthError> {
    let connection = Connection::open_with_flags(
        path,
        OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_NO_MUTEX,
    )
    .map_err(TruthError::Open)?;
    read(&connection, wanted)
}

fn read(connection: &Connection, wanted: &BTreeSet<&str>) -> Result<TruthDatabase, TruthError> {
    let info = connection
        .prepare("SELECT key, value FROM info")
        .and_then(|mut statement| {
            statement
                .query_map([], |row| {
                    Ok((row.get::<_, String>(0)?, row.get::<_, String>(1)?))
                })?
                .collect::<Result<BTreeMap<_, _>, _>>()
        })
        .map_err(TruthError::Query)?;
    let version = info.get("schema_version").cloned().unwrap_or_default();
    if version != SCHEMA_VERSION {
        return Err(TruthError::Schema {
            found: version,
            expected: SCHEMA_VERSION,
        });
    }

    let mut batch = connection
        .prepare("SELECT batch FROM samples WHERE sample_id = ?1")
        .map_err(TruthError::Query)?;
    let mut alleles = connection
        .prepare("SELECT pos, ref, seq FROM variants WHERE sample_id = ?1 ORDER BY rowid")
        .map_err(TruthError::Query)?;
    let mut intervals = connection
        .prepare(
            "SELECT start_pos, end_pos FROM intervals WHERE sample_id = ?1 \
             ORDER BY start_pos, end_pos",
        )
        .map_err(TruthError::Query)?;
    let mut status = connection
        .prepare("SELECT status FROM sample_status WHERE sample_id = ?1 AND batch = ?2")
        .map_err(TruthError::Query)?;
    let mut samples = BTreeMap::new();
    for sample in wanted {
        let Some(sample_batch) = batch
            .query_row(params![sample], |row| row.get::<_, String>(0))
            .optional()
            .map_err(TruthError::Query)?
        else {
            continue;
        };
        let rows = alleles
            .query_map(params![sample], |row| {
                Ok(Allele::new(
                    &row.get::<_, String>(0)?,
                    &row.get::<_, String>(1)?,
                    &row.get::<_, String>(2)?,
                ))
            })
            .map_err(TruthError::Query)?;
        let mut sample_alleles = Vec::new();
        for allele in rows {
            if sample_alleles.len() == MAX_ALLELES {
                return Err(TruthError::TooManyAlleles {
                    sample: (*sample).to_owned(),
                    maximum: MAX_ALLELES,
                });
            }
            sample_alleles.push(allele.map_err(TruthError::Query)?);
        }
        let sample_intervals = intervals
            .query_map(params![sample], |row| {
                Ok((row.get::<_, i64>(0)?, row.get::<_, i64>(1)?))
            })
            .and_then(Iterator::collect::<Result<Vec<_>, _>>)
            .map_err(TruthError::Query)?;
        let sample_status = status
            .query_row(params![sample, sample_batch], |row| row.get::<_, String>(0))
            .optional()
            .map_err(TruthError::Query)?;
        samples.insert(
            (*sample).to_owned(),
            TruthSample {
                batch: sample_batch,
                status: sample_status,
                alleles: sample_alleles,
                intervals: sample_intervals,
            },
        );
    }
    Ok(TruthDatabase { info, samples })
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    /// A synthetic database with the statistics schema.
    pub(crate) const FIXTURE: &str = "
        CREATE TABLE info (key TEXT PRIMARY KEY, value TEXT NOT NULL);
        CREATE TABLE samples (sample_id TEXT PRIMARY KEY, batch TEXT NOT NULL, no_snps INTEGER,
            no_ins INTEGER, no_dels INTEGER, no_identicals INTEGER, no_unread INTEGER);
        CREATE TABLE variants (sample_id TEXT NOT NULL, kind TEXT NOT NULL, pos TEXT NOT NULL,
            base INTEGER NOT NULL, ref TEXT NOT NULL, seq TEXT NOT NULL);
        CREATE TABLE intervals (sample_id TEXT NOT NULL, region TEXT NOT NULL,
            start_pos INTEGER NOT NULL, end_pos INTEGER NOT NULL);
        CREATE TABLE sample_status (batch TEXT NOT NULL, sample_id TEXT NOT NULL,
            barcode TEXT NOT NULL, status TEXT NOT NULL, in_metadata INTEGER NOT NULL,
            in_seq INTEGER NOT NULL, seq_reject_reason TEXT NOT NULL, in_json INTEGER NOT NULL,
            in_fasta INTEGER NOT NULL, PRIMARY KEY (batch, sample_id));
        INSERT INTO info VALUES ('schema_version', '1'), ('samples', '3');
        INSERT INTO samples VALUES ('S1', 'b1', 2, 1, 1, 0, 0), ('S2', 'b2', 0, 0, 0, 0, 0),
            ('S3', 'b1', 1, 0, 0, 0, 0);
        INSERT INTO variants VALUES ('S1', 'snp', '263', 263, 'A', 'G'),
            ('S1', 'insertion', '309.1', 309, '-', 'C'), ('S1', 'deletion', '249', 249, 'A', '-'),
            ('S1', 'snp', '16189', 16189, 'T', 'Y'), ('S3', 'snp', '73', 73, 'A', 'G');
        INSERT INTO intervals VALUES ('S1', 'HV2', 73, 340), ('S1', 'HV1', 16024, 16365);
        INSERT INTO sample_status VALUES ('b1', 'S1', 'x', 'OK', 1, 1, '', 1, 1),
            ('b9', 'S2', 'x', 'OK', 1, 1, '', 1, 1);
    ";

    /// Builds the fixture database in a file under `directory`.
    pub(crate) fn database(
        directory: &Path,
        sql: &str,
    ) -> Result<std::path::PathBuf, rusqlite::Error> {
        let path = directory.join("statistics.sqlite");
        Connection::open(&path)?.execute_batch(sql)?;
        Ok(path)
    }

    #[test]
    fn reads_the_panel_samples_read_only() -> Result<(), Box<dyn std::error::Error>> {
        let directory = tempfile::tempdir()?;
        let path = database(directory.path(), FIXTURE)?;
        let truth = load(&path, &BTreeSet::from(["S1", "S2", "S9"]))?;
        assert_eq!(truth.info["samples"], "3");
        assert_eq!(truth.samples.len(), 2);
        let first = &truth.samples["S1"];
        assert_eq!(first.batch, "b1");
        assert_eq!(first.status.as_deref(), Some("OK"));
        assert_eq!(
            first.alleles,
            [
                Allele::new("263", "A", "G"),
                Allele::new("309.1", "-", "C"),
                Allele::new("249", "A", "-"),
                Allele::new("16189", "T", "Y"),
            ]
        );
        assert_eq!(first.intervals, [(73, 340), (16024, 16365)]);
        // S2's status row belongs to another batch.
        assert_eq!(truth.samples["S2"].status, None);
        assert!(truth.samples["S2"].alleles.is_empty());
        // The connection cannot write.
        let connection = Connection::open_with_flags(&path, OpenFlags::SQLITE_OPEN_READ_ONLY)?;
        assert!(connection.execute("DELETE FROM variants", []).is_err());
        Ok(())
    }

    #[test]
    fn digests_exactly_the_truth_the_harness_uses() {
        let record = TruthSample {
            batch: "b1".into(),
            status: Some("OK".into()),
            alleles: vec![
                Allele::new("263", "A", "G"),
                Allele::new("16189.0", "T", "C"),
            ],
            intervals: vec![(16024, 16365), (73, 340)],
        };
        let reordered = TruthSample {
            alleles: vec![Allele::new("16189", "T", "C"), Allele::new("263", "A", "G")],
            intervals: vec![(73, 340), (16024, 16365)],
            ..record.clone()
        };
        // Row order and position writing do not matter; every field does.
        assert_eq!(digest(Some(&record)), digest(Some(&reordered)));
        let changes = [
            TruthSample {
                batch: "b2".into(),
                ..record.clone()
            },
            TruthSample {
                status: None,
                ..record.clone()
            },
            TruthSample {
                alleles: vec![Allele::new("263", "A", "G")],
                ..record.clone()
            },
            TruthSample {
                intervals: vec![(73, 341), (16024, 16365)],
                ..record.clone()
            },
        ];
        for changed in &changes {
            assert_ne!(digest(Some(&record)), digest(Some(changed)));
        }
        assert_ne!(digest(None), digest(Some(&record)));
        assert_eq!(digest(None).len(), 64);
        let first = digest(Some(&record));
        let second = digest(None);
        assert_ne!(
            panel_digest([("S1", first.as_str()), ("S2", second.as_str())]),
            panel_digest([("S2", second.as_str()), ("S1", first.as_str())])
        );
    }

    #[test]
    fn finds_the_samples_whose_truth_changed() -> Result<(), Box<dyn std::error::Error>> {
        let directory = tempfile::tempdir()?;
        let wanted = BTreeSet::from(["S1", "S2", "S3"]);
        let path = database(directory.path(), FIXTURE)?;
        let built = load(&path, &wanted)?;
        let parent = crate::manifest::parse(
            crate::manifest::tests::manifest(&[
                "S1\tb1\tclean\t1\t1\t0\t/t/a.ab1",
                "S2\tb2\tclean\t1\t1\t0\t/t/b.ab1",
                "S3\tb1\tclean\t1\t1\t0\t/t/c.ab1",
            ])
            .as_bytes(),
        )?;
        let source = crate::panel::Source {
            path: &path,
            sha256: "abc",
            database: &built,
        };
        let panel = crate::manifest::parse(
            crate::panel::derive(&parent, &source, "v3", "2026-10-15")
                .text
                .as_bytes(),
        )?;
        assert_eq!(
            TruthCheck::new(&parent, &built).state(),
            TruthState::NotRecorded
        );
        let same = TruthCheck::new(&panel, &built);
        assert_eq!(same.state(), TruthState::Unchanged);
        assert!(same.changed.is_empty());

        // A refresh that changes S3's status and one of S1's alleles.
        let refreshed = directory.path().join("refreshed");
        std::fs::create_dir(&refreshed)?;
        let sql = format!(
            "{FIXTURE} UPDATE variants SET seq = 'C' WHERE sample_id = 'S1' AND pos = '263';
             INSERT INTO sample_status VALUES ('b1', 'S3', 'x', 'HOLD', 1, 1, '', 1, 1);"
        );
        let refreshed = load(&database(&refreshed, &sql)?, &wanted)?;
        let check = TruthCheck::new(&panel, &refreshed);
        assert_eq!(check.state(), TruthState::Changed);
        let changed = check
            .changed
            .iter()
            .map(|change| change.sample.as_str())
            .collect::<Vec<_>>();
        assert_eq!(changed, ["S1", "S3"]);
        assert!(
            check
                .changes()
                .starts_with("sample_id\tstratum\tretired\trecorded_digest")
        );
        assert_eq!(check.changes().lines().count(), 3);
        Ok(())
    }

    #[test]
    fn rejects_another_schema_or_a_missing_file() -> Result<(), Box<dyn std::error::Error>> {
        let directory = tempfile::tempdir()?;
        let path = database(
            directory.path(),
            "CREATE TABLE info (key TEXT PRIMARY KEY, value TEXT NOT NULL);
             INSERT INTO info VALUES ('schema_version', '2');",
        )?;
        assert!(matches!(
            load(&path, &BTreeSet::new()),
            Err(TruthError::Schema { .. })
        ));
        assert!(matches!(
            load(&directory.path().join("missing.sqlite"), &BTreeSet::new()),
            Err(TruthError::Open(_))
        ));
        let empty = directory.path().join("empty");
        std::fs::create_dir(&empty)?;
        let path = database(&empty, "CREATE TABLE other (x INTEGER);")?;
        assert!(matches!(
            load(&path, &BTreeSet::new()),
            Err(TruthError::Query(_))
        ));
        Ok(())
    }
}
