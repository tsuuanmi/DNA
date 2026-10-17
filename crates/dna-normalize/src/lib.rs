//! DNA post-calling normalization (ADR-0069, ADR-0077): haplotype-preserving
//! representation, sequence-equivalent normalization, and target nomenclature
//! over called variants, in the lab's reviewer convention.
//!
//! It depends only on the kernel.

pub mod variant_nomenclature;
pub mod variant_normalization;
pub mod variant_representation;
