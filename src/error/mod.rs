//! Typed failures exposed by the DNA application boundary.

use std::path::PathBuf;

/// Result type returned by DNA operations.
pub type Result<T> = std::result::Result<T, Error>;

/// A failure in a validated analysis stage.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum Error {
    /// A required path is absent, invalid, or has the wrong filesystem type.
    #[error("invalid {kind} path {path}: {reason}")]
    Path {
        /// Role of the path, for example `"AB1"` or `"output"`.
        kind: &'static str,
        /// The rejected path.
        path: PathBuf,
        /// Why the path was rejected.
        reason: String,
    },
    /// A file could not be read.
    #[error("failed to read {kind} file {path}: {source}")]
    Read {
        /// Role of the file, for example `"AB1"` or `"configuration"`.
        kind: &'static str,
        /// The unreadable file.
        path: PathBuf,
        /// Underlying I/O failure.
        #[source]
        source: std::io::Error,
    },
    /// Configuration bytes were not valid TOML.
    #[error("invalid configuration {path}: {source}")]
    ConfigParse {
        /// The configuration file.
        path: PathBuf,
        /// Underlying TOML parse failure.
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
    /// Variant extraction or configured eligibility failed.
    #[error("variant calling failed: {0}")]
    Variant(String),
    /// Post-calling variant representation normalization failed.
    #[error("variant normalization failed: {0}")]
    VariantNormalization(String),
    /// Target-specific variant nomenclature representation failed.
    #[error("variant nomenclature failed: {0}")]
    VariantNomenclature(String),
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
        /// The original operation failure.
        operation: Box<Error>,
        /// The failure to record it in the operation log.
        logging: Box<Error>,
    },
    /// A log directory or append-only log file operation failed.
    #[error("failed to access log path {path}: {source}")]
    Log {
        /// The log directory or file.
        path: PathBuf,
        /// Underlying I/O failure.
        #[source]
        source: std::io::Error,
    },
    /// An output file operation failed.
    #[error("failed to publish output {path}: {source}")]
    Output {
        /// The output or temporary file.
        path: PathBuf,
        /// Underlying I/O failure.
        #[source]
        source: std::io::Error,
    },
}
