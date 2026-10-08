//! Variant call evidence lookup failures shared by sample aggregation and reporting.

/// Why a variant call mapping could not be resolved to its call evidence.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum CallEvidenceError {
    /// The mapping refers past the end of the read calls.
    #[error("variant references missing call index {index}")]
    MissingCall {
        /// Call index (0-based).
        index: usize,
    },
    /// The mapping refers past the end of the quality records.
    #[error("variant references missing quality index {index}")]
    MissingQuality {
        /// Call index (0-based).
        index: usize,
    },
    /// The call and quality records at the index disagree on their own index.
    #[error("variant call index {index} does not match call/quality records")]
    IndexMismatch {
        /// Call index (0-based).
        index: usize,
    },
    /// The call has no primary-event peak evidence.
    #[error("variant call index {index} lacks primary-event peak evidence")]
    MissingPeakEvidence {
        /// Call index (0-based).
        index: usize,
    },
}
