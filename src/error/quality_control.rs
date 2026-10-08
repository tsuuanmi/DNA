//! Relative-quality scoring and end-trimming failures.

/// Why a read could not be scored or trimmed to a sufficient retained interval.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum QualityControlError {
    /// Penalties need at least one call and a positive window.
    #[error("quality penalties require calls and a positive window")]
    EmptyPenaltyInput,
    /// The read has fewer calls than the minimum retained length.
    #[error("{calls} calls are fewer than minimum_retained_bases {minimum}")]
    TooFewCalls {
        /// Call count.
        calls: usize,
        /// Configured minimum retained bases.
        minimum: usize,
    },
    /// Trimming left fewer bases than the minimum retained length.
    #[error("retained interval {start}..{end} is shorter than minimum {minimum}")]
    RetainedTooShort {
        /// Trim start (0-based, inclusive).
        start: usize,
        /// Trim end (0-based, exclusive).
        end: usize,
        /// Configured minimum retained bases.
        minimum: usize,
    },
    /// A penalty computation overflowed.
    #[error("{0}")]
    Overflow(&'static str),
}
