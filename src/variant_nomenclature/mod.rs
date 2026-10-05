//! Optional target-specific variant nomenclature boundary.
//!
//! Nomenclature is distinct from evidence-backed variant calling and generic
//! sequence-equivalent normalization. This module currently exposes only the
//! immutable context seam that future target policies may consume.

use crate::variant_analysis::{ReferenceIdentity, Variant};
use crate::variant_normalization::VariantNormalizationResult;

/// Immutable context available to a target-specific nomenclature policy.
///
/// This borrowed view keeps reference identity, exact source calls, the complete
/// reconstructed alternate haplotype, and normalized variants together without
/// copying or rewriting any representation.
#[derive(Debug, Clone, Copy)]
pub struct NomenclatureInput<'a> {
    pub reference: &'a ReferenceIdentity,
    pub source_variants: &'a [Variant],
    pub alternate_sequence: &'a str,
    pub normalized_variants: &'a [Variant],
}

/// Borrows the complete nomenclature context from one normalization result.
#[must_use]
pub fn from_normalization(normalized: &VariantNormalizationResult) -> NomenclatureInput<'_> {
    NomenclatureInput {
        reference: &normalized.reference,
        source_variants: &normalized.source_variants,
        alternate_sequence: &normalized.alternate_sequence,
        normalized_variants: &normalized.normalized_variants,
    }
}
