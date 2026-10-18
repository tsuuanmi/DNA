//! Variant extraction, anchoring, and eligibility failures.

/// Why primary-sequence differences could not be called from an alignment.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum VariantError {
    /// An anchored reference allele disagrees with the reference sequence.
    #[error(
        "variant reference allele disagrees with the supplied reference at position {position}"
    )]
    ReferenceAlleleMismatch {
        /// 1-based reference position of the disagreeing base.
        position: usize,
    },
    /// A reported variant has no supporting trace call.
    #[error("{kind} at position {position} has no supporting calls")]
    NoSupportingCalls {
        /// Variant type label.
        kind: &'static str,
        /// 1-based variant position.
        position: usize,
    },
    /// A variant call mapping refers past the end of the read calls.
    #[error("variant filter references missing call index {index}")]
    MissingCall {
        /// Call index (0-based).
        index: usize,
    },
    /// A coordinate computation overflowed.
    #[error("{0}")]
    Overflow(&'static str),
    /// An internal alignment-to-variant invariant was violated.
    #[error("{0}")]
    Inconsistent(&'static str),
}
