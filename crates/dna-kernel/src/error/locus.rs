//! Shared Sanger locus-window geometry failures.

/// Why symmetric sample windows could not be built around the trace loci.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum LocusWindowError {
    /// Window geometry needs at least two loci.
    #[error("at least two locus positions are required")]
    TooFewLoci,
    /// A neighboring-locus midpoint overflowed.
    #[error("locus-window midpoint overflow")]
    MidpointOverflow,
    /// A window is empty, exceeds the trace, or excludes its own locus.
    #[error("invalid locus window {start}..{end} for locus position {position}")]
    InvalidWindow {
        /// Window start sample (inclusive).
        start: usize,
        /// Window end sample (exclusive).
        end: usize,
        /// Locus sample position.
        position: usize,
    },
}
