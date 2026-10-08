//! Canonical called-variant contracts shared by every capability.
//!
//! These public types are the cross-modality boundary of ADR-0058: Variant
//! Analysis produces them, and normalization and nomenclature consume them.
//! They carry no evidence; the crate-internal `model::variant` keeps the
//! call mappings and eligibility from which they are projected.

use crate::model::variant as internal;

/// Stable reference identity carried by analysis results.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReferenceIdentity {
    /// FASTA record identifier of the reference.
    pub name: String,
    /// SHA-256 of the normalized reference sequence.
    pub sha256: String,
}

/// Supported called-variant type.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[non_exhaustive]
pub enum VariantKind {
    /// Single-nucleotide substitution.
    Snv,
    /// Anchored insertion of one or more bases.
    Ins,
    /// Anchored deletion of one or more bases.
    Del,
}

/// One reportable evidence-backed primary-sequence difference.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Variant {
    /// Reference record identifier the variant is placed on.
    pub contig: String,
    /// 1-based position of the first reference allele base.
    pub position_1based: usize,
    /// Reference allele, including the anchor base for indels.
    pub reference: String,
    /// Alternate allele, including the anchor base for indels.
    pub alternate: String,
    /// Variant type.
    pub kind: VariantKind,
}

/// Cross-modality boundary for evidence-backed variants against one reference.
///
/// This type does not imply right/left alignment, nomenclature, VCF
/// normalization, or another representation policy.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CalledVariantSet {
    /// Reference every variant in the set is placed on.
    pub reference: ReferenceIdentity,
    /// Evidence-backed called variants.
    pub variants: Vec<Variant>,
}

impl From<&internal::Variant> for Variant {
    /// Projects an internal called variant into the public variant boundary.
    fn from(variant: &internal::Variant) -> Self {
        Self {
            contig: variant.contig.clone(),
            position_1based: variant.position_1based,
            reference: variant.reference.clone(),
            alternate: variant.alternate.clone(),
            kind: match variant.kind {
                internal::VariantKind::Snv => VariantKind::Snv,
                internal::VariantKind::Ins => VariantKind::Ins,
                internal::VariantKind::Del => VariantKind::Del,
            },
        }
    }
}
