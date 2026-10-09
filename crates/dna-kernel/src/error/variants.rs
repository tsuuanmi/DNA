//! Variants-document failures.

/// Why a `dna.variants/v1` document cannot be used.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum VariantsError {
    /// The file size is zero or exceeds the accepted maximum.
    #[error("variants file size {bytes} is outside 1..={maximum} bytes")]
    FileSize {
        /// Observed file size.
        bytes: u64,
        /// Largest accepted file size.
        maximum: usize,
    },
    /// The document is not a supported variants contract.
    #[error("unsupported schema version {found:?}; expected \"dna.variants/v1\"")]
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
    /// The document was called against another reference sequence.
    #[error("document reference does not match the supplied reference")]
    ReferenceMismatch,
    /// The target profile declares no notation to derive.
    #[error("the target profile declares no notation")]
    NoNotation,
    /// A variant has an unknown kind or a non-canonical allele.
    #[error("read {read:?} has an invalid variant at position {position}")]
    InvalidVariant {
        /// The read's name.
        read: String,
        /// The variant's 1-based position.
        position: usize,
    },
}
