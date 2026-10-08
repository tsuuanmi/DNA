//! Evidence-profile alignment and orientation failures.

/// Why a read could not be aligned to a unique, sufficiently supported placement.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum AlignmentError {
    /// The query or the reference is empty.
    #[error("query and reference must both be non-empty")]
    EmptyInput,
    /// The query bases and their evidence profiles differ in length.
    #[error("query/profile length mismatch: {bases} query bases, {profiles} profiles")]
    ProfileLengthMismatch {
        /// Query base count.
        bases: usize,
        /// Evidence-profile count.
        profiles: usize,
    },
    /// The dynamic-programming matrix would exceed the allocation cap.
    #[error("alignment requires {cells} cells; cap is {maximum}")]
    TooManyCells {
        /// Required matrix cells.
        cells: usize,
        /// Largest accepted matrix.
        maximum: usize,
    },
    /// No traceback satisfied the alignment bounds.
    #[error("no valid bounded alignment traceback was found")]
    NoTraceback,
    /// Forward and reverse orientations score equally.
    #[error("forward and reverse evidence-profile scores are tied")]
    OrientationTie,
    /// The selected orientation has more than one best placement.
    #[error("selected orientation has multiple equally scoring placements")]
    AmbiguousPlacement,
    /// Too few callable columns support the selected placement.
    #[error("alignment has {found} callable columns; minimum is {minimum}")]
    TooFewCallableColumns {
        /// Callable columns in the placement.
        found: usize,
        /// Configured minimum.
        minimum: usize,
    },
    /// Callable identity of the selected placement is below the configured minimum.
    #[error("alignment callable identity {identity:.4} is below {minimum:.4}")]
    LowIdentity {
        /// Observed callable identity.
        identity: f64,
        /// Configured minimum identity.
        minimum: f64,
    },
    /// Signal loci and quality records disagree in count.
    #[error("signal/quality call count mismatch: {loci} loci, {qualities} quality records")]
    CallCountMismatch {
        /// Signal locus count.
        loci: usize,
        /// Quality record count.
        qualities: usize,
    },
    /// The trim interval does not fit the locus profiles.
    #[error("invalid trim interval {start}..{end} for {profiles} locus profiles")]
    InvalidTrim {
        /// Trim start (0-based, inclusive).
        start: usize,
        /// Trim end (0-based, exclusive).
        end: usize,
        /// Locus profile count.
        profiles: usize,
    },
    /// The retained sequence and its profiles differ in length.
    #[error("retained sequence/profile length mismatch: {bases} bases, {profiles} profiles")]
    RetainedLengthMismatch {
        /// Retained base count.
        bases: usize,
        /// Retained profile count.
        profiles: usize,
    },
    /// The canonical traceback rescored differently from the dynamic program.
    #[error("traceback score {observed} disagrees with DP score {expected}")]
    ScoreMismatch {
        /// Rescored traceback score.
        observed: i64,
        /// Dynamic-programming score.
        expected: i64,
    },
    /// A canonical column refers past the end of the query profiles.
    #[error("canonical alignment query index {index} is out of profile bounds")]
    ProfileIndexOutOfBounds {
        /// Offending query index.
        index: usize,
    },
    /// A size or score computation overflowed.
    #[error("{0}")]
    Overflow(&'static str),
    /// An internal alignment invariant was violated.
    #[error("{0}")]
    Inconsistent(&'static str),
}
