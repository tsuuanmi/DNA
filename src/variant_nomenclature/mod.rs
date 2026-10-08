//! Optional target-specific variant nomenclature boundary.
//!
//! Nomenclature is distinct from evidence-backed variant calling and generic
//! sequence-equivalent normalization. The module exposes immutable context for
//! target policies and currently implements the validated human-mtDNA HVS-II
//! 309/315 representation under `mtdna`.

use crate::variant_analysis::{ReferenceIdentity, Variant};
use crate::variant_normalization::VariantNormalizationResult;

pub mod mtdna;

/// Immutable context available to a target-specific nomenclature policy.
///
/// This borrowed view keeps reference identity, exact source calls, the complete
/// reconstructed alternate haplotype, and normalized variants together without
/// copying or rewriting any representation.
#[derive(Debug, Clone, Copy)]
pub struct NomenclatureInput<'a> {
    /// Identity of the reference the haplotype was reconstructed against.
    pub reference: &'a ReferenceIdentity,
    /// Exact called variants supplied to normalization.
    pub source_variants: &'a [Variant],
    /// Reconstructed alternate haplotype.
    pub alternate_sequence: &'a str,
    /// Variants after the normalization policy.
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

/// Source, normalized, and target-represented views of one unchanged haplotype.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct VariantNomenclatureResult {
    /// Identity of the reference the haplotype was reconstructed against.
    pub reference: ReferenceIdentity,
    /// Exact called variants supplied to normalization.
    pub source_variants: Vec<Variant>,
    /// Reconstructed alternate haplotype, unchanged by nomenclature.
    pub alternate_sequence: String,
    /// Variants after the normalization policy.
    pub normalized_variants: Vec<Variant>,
    /// Variants expressed in the target nomenclature.
    pub represented_variants: Vec<Variant>,
}
