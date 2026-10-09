//! Observation-only signal-processing failures.

use super::LocusWindowError;

/// Why signal-quality evidence could not be derived from the trace.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum SignalError {
    /// Locus windows could not be built.
    #[error("{0}")]
    LocusWindow(LocusWindowError),
    /// The read has fewer calls than one feature window.
    #[error("{calls} calls are fewer than window_size_bases {window}")]
    TooFewCalls {
        /// Call count.
        calls: usize,
        /// Configured window size in bases.
        window: usize,
    },
    /// The trace has fewer loci than one feature window.
    #[error("{loci} loci are fewer than window_size_bases {window}")]
    TooFewLoci {
        /// Locus count.
        loci: usize,
        /// Configured window size in bases.
        window: usize,
    },
    /// Locus windows and loci disagree in count.
    #[error("expected {expected} locus windows, found {found}")]
    WindowCountMismatch {
        /// Locus count.
        expected: usize,
        /// Window count.
        found: usize,
    },
    /// Locus evidence and trace calls disagree in count.
    #[error("trace integrity expected {expected} loci, found {found}")]
    LocusCountMismatch {
        /// Trace call count.
        expected: usize,
        /// Locus evidence count.
        found: usize,
    },
    /// A call window spans an empty or out-of-trace sample interval.
    #[error(
        "invalid sample interval {sample_start}..{sample_end} for call window {call_start}..{call_end}"
    )]
    InvalidSampleInterval {
        /// First sample (inclusive).
        sample_start: usize,
        /// End sample (exclusive).
        sample_end: usize,
        /// First call (inclusive).
        call_start: usize,
        /// End call (exclusive).
        call_end: usize,
    },
    /// A locus context spans an empty or out-of-trace sample interval.
    #[error("invalid locus context sample interval {start}..{end}")]
    InvalidContextInterval {
        /// First sample (inclusive).
        start: usize,
        /// End sample (exclusive).
        end: usize,
    },
    /// A locus event window is empty, out of trace, or excludes its locus.
    #[error("invalid locus event window {start}..{end} at canonical locus {locus}")]
    InvalidEventWindow {
        /// First sample (inclusive).
        start: usize,
        /// End sample (exclusive).
        end: usize,
        /// Canonical locus sample.
        locus: usize,
    },
    /// Derived locus evidence violates its coordinate or profile invariants.
    #[error("inconsistent locus evidence at call {call}")]
    InconsistentLocus {
        /// Call index (0-based).
        call: usize,
    },
    /// A statistic was requested over too few or non-finite samples.
    #[error("{0}")]
    InsufficientData(&'static str),
}
