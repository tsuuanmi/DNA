//! Haplotype-preserving variant representation failures.

/// Why called variants could not be converted to, applied as, or rendered
/// from reference sequence edits.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum RepresentationError {
    /// A variant is placed on a different contig than the reference.
    #[error("called variant contig does not match the supplied reference")]
    ContigMismatch,
    /// A variant position is zero.
    #[error("called variant position must be one-based")]
    PositionNotOneBased,
    /// An allele is empty or contains symbols other than uppercase DNA bases.
    #[error("called variant {allele} allele must contain uppercase DNA bases")]
    InvalidAllele {
        /// `"reference"` or `"alternate"`.
        allele: &'static str,
    },
    /// The variant kind disagrees with its allele lengths.
    #[error("called variant kind disagrees with its reference/alternate alleles")]
    KindMismatch,
    /// A variant spans the end/start seam of the reference coordinates.
    #[error("origin-spanning variants are not represented across the canonical seam")]
    OriginSpanning,
    /// A reference allele disagrees with the reference sequence.
    #[error("called variant reference allele disagrees with the supplied reference")]
    ReferenceAlleleMismatch,
    /// A variant's alleles are identical.
    #[error("called variant does not change the reference sequence")]
    NoChange,
    /// An edit lies outside the reference.
    #[error("variant edit lies outside the supplied reference")]
    EditOutsideReference,
    /// Two edits overlap on the reference.
    #[error("variant edits overlap on the reference")]
    OverlappingEdits,
    /// Two insertions share one reference boundary, so their order is undefined.
    #[error("multiple insertion edits share one reference boundary")]
    SharedInsertionBoundary,
    /// A rendered edit is neither an insertion nor a deletion nor an SNV.
    #[error("variant representation produced an unsupported replacement edit")]
    UnsupportedReplacement,
    /// A rendered deletion lies outside the reference.
    #[error("variant deletion lies outside the reference")]
    DeletionOutsideReference,
    /// An indel anchor lies outside the reference.
    #[error("variant representation anchor lies outside the reference")]
    AnchorOutsideReference,
    /// A coordinate computation overflowed or underflowed.
    #[error("{0}")]
    Overflow(&'static str),
}

/// Why post-calling normalization could not be applied.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum NormalizationError {
    /// The called-variant reference identity differs from the supplied reference.
    #[error("called variants do not match the supplied reference identity")]
    ReferenceIdentityMismatch,
    /// Normalized edits would change the reconstructed haplotype.
    #[error("normalized edits changed the reconstructed haplotype")]
    HaplotypeChanged,
    /// Variants could not be represented as sequence edits.
    #[error("{0}")]
    Representation(RepresentationError),
}

/// Why target nomenclature could not be applied.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum NomenclatureError {
    /// The nomenclature input reference identity differs from the supplied reference.
    #[error("nomenclature input does not match the supplied reference identity")]
    ReferenceIdentityMismatch,
    /// Normalized variants do not reconstruct the supplied alternate haplotype.
    #[error("normalized variants do not reproduce the supplied alternate haplotype")]
    InconsistentInput,
    /// The profile windows would change the reconstructed haplotype.
    #[error("nomenclature changed the reconstructed haplotype")]
    HaplotypeChanged,
    /// The reference does not carry a profile window's sequence.
    #[error("{window} reference window does not match the profile sequence")]
    ReferenceMotifMismatch {
        /// Window name declared by the profile, for example `"HVS-II"`.
        window: String,
    },
    /// An edit straddles a window boundary.
    #[error("variant edit crosses the {window} nomenclature window")]
    WindowCrossing {
        /// Window name declared by the profile, for example `"HVS-II"`.
        window: String,
    },
    /// A window rule would change the local haplotype.
    #[error("{window} representation changed the local haplotype")]
    LocalHaplotypeChanged {
        /// Window name declared by the profile, for example `"HVS-II"`.
        window: String,
    },
    /// Variants could not be represented as sequence edits.
    #[error("{0}")]
    Representation(RepresentationError),
}

impl From<RepresentationError> for NormalizationError {
    fn from(error: RepresentationError) -> Self {
        Self::Representation(error)
    }
}

impl From<RepresentationError> for NomenclatureError {
    fn from(error: RepresentationError) -> Self {
        Self::Representation(error)
    }
}
