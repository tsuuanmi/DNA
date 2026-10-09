//! Signal-derived read-callability failures.

/// Why a read's callability could not be derived.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum CallabilityError {
    /// Call spacing needs at least two positions.
    #[error("{positions} call positions are fewer than the minimum {minimum}")]
    TooFewPositions {
        /// Call position count.
        positions: usize,
        /// Smallest count the method handles.
        minimum: usize,
    },
    /// Per-position evidence and calls disagree in count.
    #[error("expected {expected} callability positions, found {found}")]
    PositionCountMismatch {
        /// Call count.
        expected: usize,
        /// Evidence record count.
        found: usize,
    },
    /// Too few calls are callable for the read to be analyzed.
    #[error("{callable} callable calls are fewer than minimum_callable_calls {minimum}")]
    TooFewCallableCalls {
        /// Unmasked call count.
        callable: usize,
        /// Configured minimum callable calls.
        minimum: usize,
    },
    /// Derived segments, mask, or callable span violate their invariants.
    #[error("inconsistent callability: {0}")]
    Inconsistent(&'static str),
}
