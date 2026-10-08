//! Signal-derived base re-calling failures.

use super::LocusWindowError;

/// Why bases could not be re-called from the trace signal.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum BasecallingError {
    /// Locus windows could not be built.
    #[error("{0}")]
    LocusWindow(LocusWindowError),
    /// A selected peak lies outside its own call window.
    #[error("selected peak escaped call window {start}..{end} at call {call}")]
    PeakOutsideWindow {
        /// Window start sample (inclusive).
        start: usize,
        /// Window end sample (exclusive).
        end: usize,
        /// Call index (0-based).
        call: usize,
    },
}
