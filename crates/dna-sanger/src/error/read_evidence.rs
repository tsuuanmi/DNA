//! Failures of a recorded read evidence document (`dna.read_evidence/v1`,
//! PROP-0011).

use std::path::PathBuf;

use dna_kernel::error::ForeignError;

/// Why a read evidence document cannot be read back.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum ReadEvidenceError {
    /// The file size is zero or exceeds the accepted maximum.
    #[error("read evidence file size {bytes} is outside 1..={maximum} bytes")]
    FileSize {
        /// Observed file size.
        bytes: u64,
        /// Largest accepted file size.
        maximum: u64,
    },
    /// The document bytes are not JSON of the read evidence's closed shape.
    #[error("{path}: {source}")]
    Parse {
        /// The document.
        path: PathBuf,
        /// Underlying JSON parse failure.
        #[source]
        source: ForeignError,
    },
    /// The document is not a `dna.read_evidence/v1` document.
    #[error("unsupported schema version {found:?}; expected \"dna.read_evidence/v1\"")]
    UnsupportedSchema {
        /// The schema the document names.
        found: String,
    },
    /// The document was written under another configuration than the one
    /// loaded.
    #[error("read evidence configuration {found} is not the loaded configuration {expected}")]
    Configuration {
        /// The document's configuration SHA-256.
        found: String,
        /// The loaded configuration's SHA-256.
        expected: String,
    },
    /// The document was written by other plugins than this build's evidence
    /// workflow.
    #[error("read evidence plugins {found} are not this build's {expected}")]
    Plugins {
        /// The document's plugins, as `id@version`.
        found: String,
        /// This build's evidence plugins, as `id@version`.
        expected: String,
    },
    /// A reason label the Sanger modality does not raise.
    #[error("unknown reason label {label:?}")]
    UnknownReason {
        /// The label as written.
        label: String,
    },
    /// A call is not where call order puts it.
    #[error("call at position {position} has index {index}")]
    CallOrder {
        /// Position of the call in the document.
        position: usize,
        /// The index the call declares.
        index: usize,
    },
    /// A value contradicts the contract or another value of the document.
    #[error("{0}")]
    Inconsistent(&'static str),
}
