//! Consensus-document failures.

/// Why a `dna.consensus/v1` document cannot be called.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum ConsensusDocumentError {
    /// The file size is zero or exceeds the accepted maximum.
    #[error("consensus file size {bytes} is outside 1..={maximum} bytes")]
    FileSize {
        /// Observed file size.
        bytes: u64,
        /// Largest accepted file size.
        maximum: usize,
    },
    /// The document is not a supported consensus contract.
    #[error("unsupported schema version {found:?}; expected \"dna.consensus/v1\"")]
    UnsupportedSchema {
        /// The declared schema version.
        found: String,
    },
    /// The document belongs to another sample.
    #[error("document sample {found:?} is not the requested sample {expected:?}")]
    SampleMismatch {
        /// The requested sample identifier.
        expected: String,
        /// The document's sample identifier.
        found: String,
    },
    /// The document was built against another reference sequence.
    #[error("document reference does not match the supplied reference")]
    ReferenceMismatch,
    /// A consensus document was given together with other sequence inputs.
    #[error("a consensus document is called on its own, without other sequence files")]
    Mixed,
    /// A segment or site is outside the reference, out of order, or does not
    /// agree with the segment sequence.
    #[error("{0}")]
    Inconsistent(&'static str),
}
