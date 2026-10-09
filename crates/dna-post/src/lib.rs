//! DNA post-calling plugins (ADR-0069): haplotype-preserving representation,
//! normalization, target nomenclature, and notation-convention conformance
//! over called variants.
//!
//! It depends only on the kernel.

pub mod conformance;
pub mod variant_nomenclature;
pub mod variant_normalization;
pub mod variant_representation;
