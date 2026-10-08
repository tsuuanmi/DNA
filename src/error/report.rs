//! Result assembly failures.

/// Why a completed analysis could not be projected into its result contract.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum ReportError {
    /// A sample read has no UTF-8 file-name stem to display.
    #[error("sample read has no valid UTF-8 filename stem")]
    InvalidReadName,
    /// Two sample reads share one display name.
    #[error("sample read name {name:?} is not unique")]
    DuplicateReadName {
        /// The shared name.
        name: String,
    },
    /// Sample evidence refers to a read index that does not exist.
    #[error("sample evidence references missing read {index}")]
    MissingRead {
        /// Read index.
        index: usize,
    },
    /// A variant call could not be resolved to its call evidence.
    #[error("{0}")]
    CallEvidence(super::CallEvidenceError),
    /// A completed model violates its assembly invariants.
    #[error("{0}")]
    Inconsistent(&'static str),
}
