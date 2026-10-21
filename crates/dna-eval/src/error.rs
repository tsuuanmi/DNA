//! Typed failures of an evaluation run.
//!
//! [`Error`] stops the run: an input that cannot be trusted as a whole (the
//! manifest, the truth database) or an output that cannot be
//! written.
//! A sample that cannot be called or compared is not an error of the run; it
//! becomes a [`crate::compare::Failure`] row.

use std::io;
use std::path::PathBuf;

/// Result of a harness step.
pub(crate) type Result<T> = std::result::Result<T, Error>;

/// A failure that stops the evaluation run.
#[derive(Debug, thiserror::Error)]
pub(crate) enum Error {
    /// An input file could not be read.
    #[error("cannot read {kind} {path}: {source}")]
    Read {
        kind: &'static str,
        path: PathBuf,
        #[source]
        source: io::Error,
    },
    /// An input file is larger than its bound.
    #[error("{kind} {path} has {bytes} bytes, more than the maximum of {maximum}")]
    TooLarge {
        kind: &'static str,
        path: PathBuf,
        bytes: u64,
        maximum: u64,
    },
    /// The panel manifest is invalid.
    #[error("invalid panel manifest: {0}")]
    Manifest(#[from] ManifestError),
    /// The truth file is invalid.
    #[error("invalid truth file: {0}")]
    Truth(#[from] TruthError),
    /// The reference cannot be loaded.
    #[error("reference {path}: {message}")]
    Reference { path: PathBuf, message: String },
    /// An earlier run's written forms are for another truth database or
    /// malformed.
    #[error("forms {path}: {message}")]
    Forms { path: PathBuf, message: String },
    /// The simulated bench cannot run.
    #[error("simulated bench: {0}")]
    SimBench(String),
    /// A batch's run context (`--run-context`) cannot be read.
    #[error("run context {path}: {message}")]
    RunContext { path: PathBuf, message: String },
    /// The forced-window list is invalid.
    #[error("invalid forced windows {path}: {message}")]
    ForcedWindows { path: PathBuf, message: String },
    /// The re-check list is invalid.
    #[error("invalid re-check list {path}: {message}")]
    Recheck { path: PathBuf, message: String },
    /// The bench case list is invalid, or names a panel or sample that no
    /// manifest holds.
    #[error("invalid bench case list {path}: {message}")]
    Cases { path: PathBuf, message: String },
    /// The current-data month is not `YYYY-MM`.
    #[error("invalid current-data month {value:?}: expected YYYY-MM")]
    CurrentSince { value: String },
    /// The review window is not a positive number of positions.
    #[error("invalid review window {window}: it must be at least 1")]
    ReviewWindow { window: i64 },

    /// A review zone width is negative.
    #[error("invalid review zone width {name} {value}: it must not be negative")]
    ZoneWidth { name: &'static str, value: i64 },

    /// The configuration's target profile cannot be read, or is not the one
    /// DNA analysed with.
    #[error("target profile {path}: {message}")]
    Profile { path: PathBuf, message: String },
    /// The baseline scores are invalid.
    #[error("invalid baseline scores: {0}")]
    Baseline(#[from] BaselineError),
    /// A command-line value is not a single token.
    #[error("invalid {name} {value:?}: it must be one token without spaces or '='")]
    Argument { name: &'static str, value: String },
    /// The thread count is outside its bound.
    #[error("thread count {threads} is outside 1..={maximum}")]
    Threads { threads: usize, maximum: usize },
    /// The output directory already holds files.
    #[error("output directory {path} is not empty")]
    OutputNotEmpty { path: PathBuf },
    /// An output could not be written.
    #[error("cannot write {path}: {source}")]
    Write {
        path: PathBuf,
        #[source]
        source: io::Error,
    },
    /// The run header could not be serialized.
    #[error("{0}")]
    Cache(#[from] crate::cache::CacheError),
    #[error("{0}")]
    Diff(#[from] DiffError),
    #[error("cannot serialize the run header: {0}")]
    Serialize(#[source] serde_json::Error),
    /// A harness thread panicked outside the per-sample guard.
    #[error("a harness thread panicked outside a sample")]
    Worker,
    /// A command was given a hold-out panel, which stays reserved.
    #[error(
        "panel {panel:?} is a hold-out panel: it stays reserved, and dna-eval reads development and readiness panels only"
    )]
    HoldOut { panel: String },
}

/// A panel manifest that does not follow the manifest format.
#[derive(Debug, thiserror::Error, PartialEq, Eq)]
pub(crate) enum ManifestError {
    /// The manifest is not UTF-8.
    #[error("the manifest is not UTF-8")]
    NotUtf8,
    /// A header value has a malformed `%` escape.
    #[error("line {line} has a malformed header value {value:?}")]
    HeaderValue { line: usize, value: String },
    /// No column header follows the comment header.
    #[error("the column header is missing")]
    MissingColumns,
    /// The column header lacks a required column.
    #[error("the column header lacks {column}")]
    MissingColumn { column: &'static str },
    /// The comment header does not declare `body_sha256`.
    #[error("the comment header does not declare body_sha256")]
    MissingBodyHash,
    /// The body does not hash to the declared `body_sha256`.
    #[error("the body has sha256 {actual}, the header declares {declared}")]
    BodyHash { declared: String, actual: String },
    /// The body does not hold the declared number of rows.
    #[error("the body has {actual} rows, the header declares {declared}")]
    RowCount { declared: String, actual: usize },
    /// The body holds more rows than the bound.
    #[error("the body has more than {maximum} rows")]
    TooManyRows { maximum: usize },
    /// A row does not have one field per column.
    #[error("line {line} has {found} fields, expected {expected}")]
    Fields {
        line: usize,
        found: usize,
        expected: usize,
    },
    /// A row has an empty sample identifier.
    #[error("line {line} has an empty sample_id")]
    EmptySample { line: usize },
    /// A sample appears on two rows.
    #[error("line {line} repeats sample {sample:?}")]
    DuplicateSample { line: usize, sample: String },
    /// A weight is not a finite non-negative number, or an anchor flag is not
    /// `0` or `1`.
    #[error("line {line} has an invalid {column} {value:?}")]
    Value {
        line: usize,
        column: &'static str,
        value: String,
    },
    /// A row lists no trace or more traces than the bound.
    #[error("line {line} lists {found} traces, expected 1..={maximum}")]
    Traces {
        line: usize,
        found: usize,
        maximum: usize,
    },
}

/// A truth database that cannot be read as the statistics schema.
#[derive(Debug, thiserror::Error)]
pub(crate) enum TruthError {
    /// The database cannot be opened read-only.
    #[error("cannot open the database read-only: {0}")]
    Open(#[source] rusqlite::Error),
    /// A query failed, for example on a missing table or column.
    #[error("query failed: {0}")]
    Query(#[source] rusqlite::Error),
    /// The database declares another schema version.
    #[error("schema version {found:?}, expected {expected:?}")]
    Schema {
        found: String,
        expected: &'static str,
    },
    /// A sample has more alleles than the bound.
    #[error("sample {sample:?} has more than {maximum} alleles")]
    TooManyAlleles { sample: String, maximum: usize },
}

/// Baseline scores that cannot be compared with this run.
#[derive(Debug, thiserror::Error, PartialEq, Eq)]
pub(crate) enum BaselineError {
    /// The scores file is not UTF-8.
    #[error("the scores file is not UTF-8")]
    NotUtf8,
    /// A required column is missing.
    #[error("the column {column} is missing")]
    MissingColumn { column: &'static str },
    /// A row does not have one field per column.
    #[error("line {line} has {found} fields, expected {expected}")]
    Fields {
        line: usize,
        found: usize,
        expected: usize,
    },
    /// A sample appears on two rows.
    #[error("line {line} repeats sample {sample:?}")]
    DuplicateSample { line: usize, sample: String },
}

/// Why two runs cannot be compared.
#[derive(Debug, thiserror::Error)]
pub(crate) enum DiffError {
    /// A run file cannot be read or parsed.
    #[error("{path}: {message}")]
    Run {
        /// The file.
        path: PathBuf,
        /// What is wrong with it.
        message: String,
    },
    /// A run's files count other passes than its own `run.json` records.
    #[error("{0}: unreviewed.tsv and scores.tsv count other passes than run.json records ({1})")]
    Inconsistent(PathBuf, String),
    /// The runs measure different things.
    #[error("the runs are not comparable: {0}; pass --force to compare them anyway")]
    Mismatch(String),
}
