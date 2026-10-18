//! Relative-quality scoring and trimming failures.

/// Why a read could not be scored or trimmed.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum QualityControlError {
    /// Penalties need at least one call and a positive window.
    #[error("quality penalties require calls and a positive window")]
    EmptyPenaltyInput,
    /// The callable span does not fit the read's calls.
    #[error("callable span {start}..{end} does not fit {calls} calls")]
    InvalidCallableSpan {
        /// Span start (0-based, inclusive).
        start: usize,
        /// Span end (0-based, exclusive).
        end: usize,
        /// Call count.
        calls: usize,
    },
    /// A penalty computation overflowed.
    #[error("{0}")]
    Overflow(&'static str),
}
