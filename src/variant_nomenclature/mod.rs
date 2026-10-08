//! Optional target-specific variant nomenclature boundary.
//!
//! Nomenclature is distinct from evidence-backed variant calling and generic
//! sequence-equivalent normalization. Target knowledge — the reference windows
//! and their ordered representation rules — comes from a
//! [`Profile`]; this module owns only the generic
//! engine and its haplotype-preservation proof.

mod windows;

use std::path::Path;

use crate::error::{NomenclatureError, Result};
use crate::model::reference::Reference;
use crate::profile::Profile;
use crate::reference;
use crate::variant_analysis::{ReferenceIdentity, Variant};
use crate::variant_normalization::VariantNormalizationResult;
use crate::variant_representation::{apply_edits, render_edits, variants_to_edits};

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

/// Applies a target profile's nomenclature windows to one normalized haplotype.
///
/// Inside each profile window, sequence-equivalent haplotypes are expressed by
/// the first window rule whose candidate reproduces the window haplotype (for
/// the human-mtDNA profile, forensic forms such as `309.1C 315.1C`,
/// `513A 523DEL 524DEL`, `16183C 16184A 16189C`). Normalized variants outside
/// the windows, and window haplotypes no rule represents, are retained
/// unchanged. A profile without windows represents the normalized variants.
///
/// # Errors
///
/// Returns [`Error`](crate::error::Error) when the reference cannot be loaded,
/// is not the profile's reference, does not match the input reference identity
/// or a profile window sequence, an edit crosses a window boundary, or the input
/// variants do not reconstruct the supplied haplotype.
pub fn apply(
    reference_path: &Path,
    profile: &Profile,
    input: NomenclatureInput<'_>,
) -> Result<VariantNomenclatureResult> {
    let reference = reference::load(reference_path, profile.topology)?;
    profile.require_reference(&reference)?;
    apply_with(&reference, profile, input)
}

/// Applies a profile's windows against an already loaded and profile-checked
/// reference; see [`apply`].
pub(crate) fn apply_with(
    reference: &Reference,
    profile: &Profile,
    input: NomenclatureInput<'_>,
) -> Result<VariantNomenclatureResult> {
    if input.reference.name != reference.name || input.reference.sha256 != reference.sequence_sha256
    {
        return Err(NomenclatureError::ReferenceIdentityMismatch.into());
    }

    let normalized_edits = variants_to_edits(
        &reference.name,
        &reference.sequence,
        input.normalized_variants,
    )
    .map_err(NomenclatureError::from)?;
    if apply_edits(&reference.sequence, &normalized_edits).map_err(NomenclatureError::from)?
        != input.alternate_sequence
    {
        return Err(NomenclatureError::InconsistentInput.into());
    }

    let represented_edits =
        windows::represent_windows(&profile.windows, &reference.sequence, &normalized_edits)?;
    if apply_edits(&reference.sequence, &represented_edits).map_err(NomenclatureError::from)?
        != input.alternate_sequence
    {
        return Err(NomenclatureError::HaplotypeChanged.into());
    }

    let represented_variants =
        render_edits(&reference.name, &reference.sequence, &represented_edits)
            .map_err(NomenclatureError::from)?;

    Ok(VariantNomenclatureResult {
        reference: input.reference.clone(),
        source_variants: input.source_variants.to_vec(),
        alternate_sequence: input.alternate_sequence.to_owned(),
        normalized_variants: input.normalized_variants.to_vec(),
        represented_variants,
    })
}
