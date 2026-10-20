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
    /// A verdict was given more reasons than it can hold.
    #[error("a verdict holds at most {maximum} reasons")]
    TooManyReasons {
        /// Largest number of reasons one verdict holds.
        maximum: usize,
    },
    /// A call's verdict holds a reason the read's reason order does not
    /// list.
    #[error("call {index} holds an undeclared reason")]
    UndeclaredReason {
        /// Offending call index (0-based).
        index: usize,
    },
    /// A phase mixture does not fit the read's calls (ADR-0082).
    #[error("phase mixture {index} does not fit the read's calls")]
    InvalidMixture {
        /// Offending mixture index (0-based).
        index: usize,
    },
    /// An unreliable span does not fit the read's calls, or its reason label
    /// is malformed.
    #[error("unreliable span {index} does not fit the read's calls")]
    InvalidUnreliable {
        /// Offending span index (0-based).
        index: usize,
    },
    /// A modality reason label is malformed, repeated, or reserved by the core.
    #[error("invalid modality reason label {0:?}")]
    InvalidReason(&'static str),
}
