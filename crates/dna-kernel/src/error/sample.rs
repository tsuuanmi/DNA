//! Multi-read sample evidence failures.

/// Why reads could not be aggregated into consistent sample evidence.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum SampleError {
    /// The sample command received no trace.
    #[error("sample analysis requires at least one AB1 trace")]
    NoTraces,
    /// The sample identifier is not a safe file-name component.
    #[error(
        "sample id must be 1..=128 ASCII characters, start with an alphanumeric character, and contain only alphanumeric, '_', '.', or '-'"
    )]
    InvalidSampleId,
    /// Aggregation received no read observation.
    #[error("at least one read observation is required")]
    NoReads,
    /// Every read was rejected: too few callable calls, or no placement.
    #[error("all {rejected} reads were rejected and none could be analyzed")]
    NoAdmittedReads {
        /// Rejected read count.
        rejected: usize,
    },
    /// Reads were analyzed against different references.
    #[error("all reads must use the same reference identity")]
    MixedReference,
    /// Reads were analyzed with different scientific configurations.
    #[error("all reads must use the same scientific configuration identity")]
    MixedConfiguration,
    /// The same trace content was supplied more than once.
    #[error("duplicate input trace content cannot contribute twice")]
    DuplicateTrace,
    /// A read has no aligned reference segment.
    #[error("read {read} has no mapped reference coverage")]
    NoCoverage {
        /// Read content identity.
        read: String,
    },
    /// A read has an empty or reversed reference segment.
    #[error("read {read} contains empty or reversed reference coverage")]
    InvalidCoverage {
        /// Read content identity.
        read: String,
    },
    /// A read has overlapping reference segments.
    #[error("read {read} contains overlapping reference segments")]
    OverlappingSegments {
        /// Read content identity.
        read: String,
    },
    /// A read maps one reference coordinate twice.
    #[error("read {read} contains duplicate reference coordinate {position}")]
    DuplicateCoordinate {
        /// Read content identity.
        read: String,
        /// 1-based reference position.
        position: usize,
    },
    /// A read reports the same normalized variant twice.
    #[error("read {read} contains duplicate normalized variant identity")]
    DuplicateVariant {
        /// Read content identity.
        read: String,
    },
    /// Reads disagree on the reference base at one coordinate.
    #[error("reference base disagrees at position {position}")]
    ReferenceMismatch {
        /// 1-based reference position.
        position: usize,
    },
    /// Aggregated evidence refers to a read index that does not exist.
    #[error("{context} references missing read {index}")]
    MissingRead {
        /// Evidence kind holding the reference.
        context: &'static str,
        /// Read index.
        index: usize,
    },
    /// An eligible nucleotide observation has no profile evidence.
    #[error("eligible nucleotide observation for read {read} lacks profile evidence")]
    MissingProfile {
        /// Read index.
        read: usize,
    },
    /// A call has no locus evidence.
    #[error("call index {call} lacks matching locus evidence")]
    MissingLocusEvidence {
        /// Call index (0-based).
        call: usize,
    },
    /// A variant call could not be resolved to its call evidence.
    #[error("{0}")]
    CallEvidence(super::CallEvidenceError),
    /// A count or coordinate computation overflowed.
    #[error("{0}")]
    Overflow(&'static str),
    /// An internal aggregation invariant was violated.
    #[error("{0}")]
    Inconsistent(&'static str),
}
