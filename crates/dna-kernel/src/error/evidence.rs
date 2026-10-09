//! Modality-to-core read-evidence failures.

/// Why a modality's per-read evidence cannot enter the core.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum EvidenceError {
    /// The informative interval does not fit the read's calls.
    #[error("informative interval {start}..{end} does not fit {calls} calls")]
    InvalidInformativeInterval {
        /// Interval start (0-based, inclusive).
        start: usize,
        /// Interval end (0-based, exclusive).
        end: usize,
        /// Call count.
        calls: usize,
    },
    /// More support vetoes were declared than a veto set can hold.
    #[error("{count} support vetoes exceed the maximum {maximum}")]
    TooManyVetoes {
        /// Declared veto count.
        count: usize,
        /// Largest supported vocabulary.
        maximum: usize,
    },
    /// A call sets a veto the read did not declare.
    #[error("call {index} sets an undeclared support veto")]
    UndeclaredVeto {
        /// Offending call index (0-based).
        index: usize,
    },
    /// A modality reason label is malformed, repeated, or reserved by the core.
    #[error("invalid modality reason label {0:?}")]
    InvalidReason(&'static str),
}
