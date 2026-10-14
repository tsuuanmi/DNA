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

use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

use rusqlite::{Connection, OpenFlags, OptionalExtension, params};

use crate::compare::Allele;
use crate::error::TruthError;

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
