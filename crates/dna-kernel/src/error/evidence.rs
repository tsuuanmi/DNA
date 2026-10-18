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
    /// A call's verdict holds a reason the read's reason order does not
    /// list.
    #[error("call {index} holds an undeclared reason")]
    UndeclaredReason {
        /// Offending call index (0-based).
        index: usize,
    },
    /// A modality reason label is malformed, repeated, or reserved by the core.
    #[error("invalid modality reason label {0:?}")]
    InvalidReason(&'static str),
}
