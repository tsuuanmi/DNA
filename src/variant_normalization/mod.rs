//! Optional post-calling variant normalization.
//!
//! Variant calling establishes evidence-backed biological differences. This
//! capability may then choose another sequence-equivalent representation under
//! an explicit policy while preserving the source calls and complete haplotype.

use crate::plugin::{Contract, PluginDescriptor, PluginFamily};
mod right;

use std::path::Path;

use crate::error::{NormalizationError, Result};
use crate::model::reference::{Reference, ReferenceTopology};
use crate::profile::IndelPlacement;
use crate::reference;
use crate::variant::{CalledVariantSet, ReferenceIdentity, Variant};
use crate::variant_representation::{apply_edits, render_edits, variants_to_edits};

/// Explicit post-calling normalization policy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum NormalizationPolicy {
    /// 3'/right-most sequence-equivalent indel placement.
    ///
    /// The FASTA boundaries define a fixed coordinate seam; equivalent events
    /// are never rotated across the end/start boundary of a circular reference.
    RightAligned,
}

impl NormalizationPolicy {
    /// The policy implementing a target profile's indel placement.
    pub(crate) fn for_placement(placement: IndelPlacement) -> Self {
        match placement {
            IndelPlacement::Right => Self::RightAligned,
        }
    }
}

/// Source and normalized representations of one unchanged called haplotype.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct VariantNormalizationResult {
    /// Identity of the reference used to reconstruct and normalize this haplotype.
    pub reference: ReferenceIdentity,
    /// Exact called variants supplied to normalization.
    pub source_variants: Vec<Variant>,
    /// Reconstructed alternate haplotype before representation movement.
    pub alternate_sequence: String,
    /// Variants after the selected representation policy.
    pub normalized_variants: Vec<Variant>,
}

/// Applies an explicit representation policy to evidence-backed called variants.
///
/// The supplied reference must have the same name and sequence identity carried
/// by the called-variant set. Normalization never changes the reconstructed
/// haplotype and preserves the original source variants in the result.
///
/// # Errors
///
/// Returns [`Error`](crate::error::Error) when the reference cannot be loaded,
/// does not match the called-variant reference identity, or a called variant
/// is inconsistent with the reference sequence.
pub fn normalize(
    reference_path: &Path,
    called: &CalledVariantSet,
    policy: NormalizationPolicy,
) -> Result<VariantNormalizationResult> {
    // Normalization never moves an edit across the FASTA seam, so the reference
    // topology does not affect it.
    let reference = reference::load(reference_path, ReferenceTopology::Linear)?;
    normalize_with(&reference, called, policy)
}

/// Normalizes against an already loaded reference; see [`normalize`].
pub(crate) fn normalize_with(
    reference: &Reference,
    called: &CalledVariantSet,
    policy: NormalizationPolicy,
) -> Result<VariantNormalizationResult> {
    if called.reference.name != reference.name
        || called.reference.sha256 != reference.sequence_sha256
    {
        return Err(NormalizationError::ReferenceIdentityMismatch.into());
    }

    let source_variants = called.variants.clone();
    let source_edits = variants_to_edits(&reference.name, &reference.sequence, &source_variants)
        .map_err(NormalizationError::from)?;
    let alternate_sequence =
        apply_edits(&reference.sequence, &source_edits).map_err(NormalizationError::from)?;

    let normalized_edits = match policy {
        NormalizationPolicy::RightAligned => {
            right::right_align(&reference.sequence, &alternate_sequence, &source_edits)
                .map_err(NormalizationError::from)?
        }
    };

    if apply_edits(&reference.sequence, &normalized_edits).map_err(NormalizationError::from)?
        != alternate_sequence
    {
        return Err(NormalizationError::HaplotypeChanged.into());
    }

    let normalized_variants = render_edits(&reference.name, &reference.sequence, &normalized_edits)
        .map_err(NormalizationError::from)?;

    Ok(VariantNormalizationResult {
        reference: called.reference.clone(),
        source_variants,
        alternate_sequence,
        normalized_variants,
    })
}

/// Haplotype-preserving variant normalization.
pub(crate) const PLUGIN: PluginDescriptor = PluginDescriptor {
    id: "normalization",
    family: PluginFamily::PostCalling,
    version: 1,
    provides: &[Contract::NormalizedVariants],
    requires: &[Contract::CalledVariants],
    config_sections: &[],
};
