//! Human-mtDNA target-specific representation policies.

use std::path::Path;

use crate::error::{NomenclatureError, Result};
use crate::model::reference::ReferenceTopology;
use crate::reference;
use crate::variant_representation::{
    SequenceEdit, apply_edits, render_edits, sort_edits, variants_to_edits,
};

use super::{NomenclatureInput, VariantNomenclatureResult};

const HV2_START_0BASED: usize = 302;
const HV2_END_0BASED_EXCLUSIVE: usize = 315;
const HV2_REFERENCE: &str = "CCCCCCCTCCCCC";
const HV2_ANCHOR_INDEX: usize = 7;
const HV2_RIGHT_RUN_LENGTH: usize = 5;

/// Applies the validated human-mtDNA HVS-II 309/315 poly-C representation.
///
/// This first nomenclature rule recognizes sequence-equivalent movement of the
/// T310 anchor inside the rCRS C-runs and expresses that haplotype as run-length
/// changes at the 309 and 315 boundaries. Other normalized variants are retained
/// unchanged.
///
/// # Errors
///
/// Returns [`Error`](crate::error::Error) when the reference cannot be loaded,
/// does not match the input reference identity, or the input variants do not
/// reconstruct the supplied haplotype.
pub fn apply_hv2_polyc(
    reference_path: &Path,
    input: NomenclatureInput<'_>,
) -> Result<VariantNomenclatureResult> {
    let reference = reference::load(reference_path, ReferenceTopology::Circular)?;
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

    let represented_edits = represent_hv2_polyc(&reference.sequence, &normalized_edits)?;
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

fn represent_hv2_polyc(reference: &str, normalized: &[SequenceEdit]) -> Result<Vec<SequenceEdit>> {
    if reference.get(HV2_START_0BASED..HV2_END_0BASED_EXCLUSIVE) != Some(HV2_REFERENCE) {
        return Err(NomenclatureError::ReferenceMotifMismatch.into());
    }

    let mut outside = Vec::new();
    let mut inside = Vec::new();
    for edit in normalized {
        if intersects_hv2(edit) {
            if !contained_in_hv2(edit) {
                return Err(NomenclatureError::WindowCrossing.into());
            }
            inside.push(SequenceEdit {
                start: edit.start - HV2_START_0BASED,
                end: edit.end - HV2_START_0BASED,
                alternate: edit.alternate.clone(),
            });
        } else {
            outside.push(edit.clone());
        }
    }

    if inside.is_empty() {
        return Ok(normalized.to_vec());
    }

    let local_alternate = apply_edits(HV2_REFERENCE, &inside).map_err(NomenclatureError::from)?;
    let Some(local_represented) = anchored_run_length_representation(&local_alternate) else {
        return Ok(normalized.to_vec());
    };

    if apply_edits(HV2_REFERENCE, &local_represented).map_err(NomenclatureError::from)?
        != local_alternate
    {
        return Err(NomenclatureError::LocalHaplotypeChanged.into());
    }

    outside.extend(local_represented.into_iter().map(|edit| SequenceEdit {
        start: edit.start + HV2_START_0BASED,
        end: edit.end + HV2_START_0BASED,
        alternate: edit.alternate,
    }));
    sort_edits(&mut outside);
    Ok(outside)
}

fn anchored_run_length_representation(alternate: &str) -> Option<Vec<SequenceEdit>> {
    if alternate.bytes().any(|base| !matches!(base, b'C' | b'T'))
        || alternate.bytes().filter(|base| *base == b'T').count() != 1
    {
        return None;
    }

    let alternate_anchor = alternate.find('T')?;
    let alternate_right_run = alternate.len().checked_sub(alternate_anchor + 1)?;
    let mut edits = Vec::new();

    if alternate_anchor < HV2_ANCHOR_INDEX {
        edits.push(SequenceEdit {
            start: alternate_anchor,
            end: HV2_ANCHOR_INDEX,
            alternate: String::new(),
        });
    } else if alternate_anchor > HV2_ANCHOR_INDEX {
        edits.push(SequenceEdit {
            start: HV2_ANCHOR_INDEX,
            end: HV2_ANCHOR_INDEX,
            alternate: "C".repeat(alternate_anchor - HV2_ANCHOR_INDEX),
        });
    }

    if alternate_right_run < HV2_RIGHT_RUN_LENGTH {
        edits.push(SequenceEdit {
            start: HV2_REFERENCE.len() - (HV2_RIGHT_RUN_LENGTH - alternate_right_run),
            end: HV2_REFERENCE.len(),
            alternate: String::new(),
        });
    } else if alternate_right_run > HV2_RIGHT_RUN_LENGTH {
        edits.push(SequenceEdit {
            start: HV2_REFERENCE.len(),
            end: HV2_REFERENCE.len(),
            alternate: "C".repeat(alternate_right_run - HV2_RIGHT_RUN_LENGTH),
        });
    }

    if edits.is_empty() {
        None
    } else {
        sort_edits(&mut edits);
        Some(edits)
    }
}

fn intersects_hv2(edit: &SequenceEdit) -> bool {
    if edit.start == edit.end {
        return edit.start > HV2_START_0BASED && edit.start <= HV2_END_0BASED_EXCLUSIVE;
    }

    edit.start < HV2_END_0BASED_EXCLUSIVE && edit.end > HV2_START_0BASED
}

fn contained_in_hv2(edit: &SequenceEdit) -> bool {
    if edit.start == edit.end {
        return edit.start > HV2_START_0BASED && edit.start <= HV2_END_0BASED_EXCLUSIVE;
    }

    edit.start >= HV2_START_0BASED && edit.end <= HV2_END_0BASED_EXCLUSIVE
}
