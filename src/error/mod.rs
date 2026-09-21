//! Typed failures exposed by the DNA application boundary.

use std::path::PathBuf;

/// Result type returned by DNA operations.
pub type Result<T> = std::result::Result<T, Error>;

/// A failure in a validated analysis stage.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// A required path is absent, invalid, or has the wrong filesystem type.
    #[error("invalid {kind} path {path}: {reason}")]
    Path {
        kind: &'static str,
        path: PathBuf,
        reason: String,
    },
    /// A file could not be read.
    #[error("failed to read {kind} file {path}: {source}")]
    Read {
        kind: &'static str,
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },
    /// Configuration bytes were not valid TOML.
    #[error("invalid configuration {path}: {source}")]
    ConfigParse {
        path: PathBuf,
        #[source]
        source: toml::de::Error,
    },
    /// A parsed configuration value violates the scientific contract.
    #[error("invalid configuration value: {0}")]
    Config(String),
    /// The ABIF container or one of its required records is invalid.
    #[error("invalid ABIF input: {0}")]
    Abif(String),
    /// The reference FASTA is invalid.
    #[error("invalid reference FASTA: {0}")]
    Fasta(String),
    /// DNA-derived base re-calling failed.
    #[error("base re-calling failed: {0}")]
    Basecalling(String),
    /// Observational signal-quality feature extraction failed.
    #[error("signal processing failed: {0}")]
    DNAProcessing(String),
    /// Quality scoring or end trimming failed.
    #[error("quality control failed: {0}")]
    QualityControl(String),
    /// Pairwise alignment failed or was not uniquely interpretable.
    #[error("alignment failed: {0}")]
    Alignment(String),
    /// Variant extraction or normalization failed.
    #[error("variant calling failed: {0}")]
    Variant(String),
    /// Sample-level read evidence is inconsistent or invalid.
    #[error("sample evidence failed: {0}")]
    Sample(String),
    /// A completed model could not be assembled consistently.
    #[error("failed to assemble analysis report: {0}")]
    Report(String),
    /// A completed model could not be serialized.
    #[error("failed to serialize result JSON: {0}")]
    Serialize(#[from] serde_json::Error),
    /// An operation failed and its terminal error record also could not be persisted.
    #[error("{operation}; additionally failed to persist the operation error log: {logging}")]
    OperationAndLog {
        operation: Box<Error>,
        logging: Box<Error>,
    },
    /// A log directory or append-only log file operation failed.
    #[error("failed to access log path {path}: {source}")]
    Log {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },
    /// An output file operation failed.
    #[error("failed to publish output {path}: {source}")]
    Output {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },
}
