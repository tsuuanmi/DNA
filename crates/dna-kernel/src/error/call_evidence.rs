//! Variant call evidence lookup failures of the core's call records.

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
    /// The call's primary base is unresolved (not A/C/G/T).
    #[error("variant call index {index} has an unresolved primary base")]
    UnresolvedCall {
        /// Call index (0-based).
        index: usize,
    },
    /// No call of the variant is resolved.
    #[error("variant has no resolved call")]
    NoResolvedCalls,
}
