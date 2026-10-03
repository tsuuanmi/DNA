//! Optional post-calling variant normalization.
//!
//! Variant calling establishes evidence-backed biological differences. This
//! capability may then choose another sequence-equivalent representation under
//! an explicit policy while preserving the source calls and complete haplotype.

use std::collections::BTreeSet;
use std::path::Path;

use crate::error::{Error, Result};
use crate::model::reference::ReferenceTopology;
use crate::reference;
use crate::variant_analysis::{CalledVariantSet, Variant, VariantKind};

/// Explicit post-calling normalization policy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NormalizationPolicy {
    /// Human-mtDNA 3'/right-most sequence-equivalent indel placement.
    ///
    /// The FASTA boundaries define the fixed rCRS coordinate seam; equivalent
    /// events are never rotated across the end/start boundary.
    MtDnaRightAligned,
}

/// Source and normalized representations of one unchanged called haplotype.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VariantNormalizationResult {
    /// Exact called variants supplied to normalization.
    pub source_variants: Vec<Variant>,
    /// Reconstructed alternate haplotype before representation movement.
    pub alternate_sequence: String,
    /// Variants after the selected representation policy.
    pub normalized_variants: Vec<Variant>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct SequenceEdit {
    start: usize,
    end: usize,
    alternate: String,
}

/// Applies an explicit representation policy to evidence-backed called variants.
///
/// The supplied reference must have the same name and sequence identity carried
/// by the called-variant set. Normalization never changes the reconstructed
/// haplotype and preserves the original source variants in the result.
pub fn normalize(
    reference_path: &Path,
    called: &CalledVariantSet,
    policy: NormalizationPolicy,
) -> Result<VariantNormalizationResult> {
    let reference = reference::load(reference_path, ReferenceTopology::Circular)?;
    if called.reference.name != reference.name
        || called.reference.sha256 != reference.sequence_sha256
    {
        return Err(normalization_error(
            "called variants do not match the supplied reference identity",
        ));
    }

    let source_variants = called.variants.clone();
    let source_edits = source_variants
        .iter()
        .map(|variant| variant_to_edit(&reference.name, &reference.sequence, variant))
        .collect::<Result<Vec<_>>>()?;
    let alternate_sequence = apply_edits(&reference.sequence, &source_edits)?;

    let normalized_edits = match policy {
        NormalizationPolicy::MtDnaRightAligned => {
            right_align_equivalent_indels(&reference.sequence, &alternate_sequence, &source_edits)?
        }
    };

    if apply_edits(&reference.sequence, &normalized_edits)? != alternate_sequence {
        return Err(normalization_error(
            "normalized edits changed the reconstructed haplotype",
        ));
    }

    let normalized_variants =
        render_edits(&reference.name, &reference.sequence, &normalized_edits)?;

    Ok(VariantNormalizationResult {
        source_variants,
        alternate_sequence,
        normalized_variants,
    })
}

fn variant_to_edit(contig: &str, reference: &str, variant: &Variant) -> Result<SequenceEdit> {
    if variant.contig != contig {
        return Err(normalization_error(
            "called variant contig does not match the supplied reference",
        ));
    }
    if variant.position_1based == 0 {
        return Err(normalization_error(
            "called variant position must be one-based",
        ));
    }
    validate_allele("reference", &variant.reference)?;
    validate_allele("alternate", &variant.alternate)?;
    validate_kind(variant)?;

    let start = variant.position_1based - 1;
    let allele_end = start
        .checked_add(variant.reference.len())
        .ok_or_else(|| normalization_error("called variant reference span overflow"))?;
    if allele_end > reference.len() {
        return Err(normalization_error(
            "origin-spanning source variants are not normalized across the canonical seam",
        ));
    }
    if reference.as_bytes().get(start..allele_end) != Some(variant.reference.as_bytes()) {
        return Err(normalization_error(
            "called variant reference allele disagrees with the supplied reference",
        ));
    }

    let reference_bytes = variant.reference.as_bytes();
    let alternate_bytes = variant.alternate.as_bytes();
    let mut prefix = 0;
    while prefix < reference_bytes.len()
        && prefix < alternate_bytes.len()
        && reference_bytes[prefix] == alternate_bytes[prefix]
    {
        prefix += 1;
    }

    let mut suffix = 0;
    while suffix < reference_bytes.len().saturating_sub(prefix)
        && suffix < alternate_bytes.len().saturating_sub(prefix)
        && reference_bytes[reference_bytes.len() - 1 - suffix]
            == alternate_bytes[alternate_bytes.len() - 1 - suffix]
    {
        suffix += 1;
    }

    let edit_start = start
        .checked_add(prefix)
        .ok_or_else(|| normalization_error("called variant edit start overflow"))?;
    let edit_end = allele_end
        .checked_sub(suffix)
        .ok_or_else(|| normalization_error("called variant edit end underflow"))?;
    let alternate_end = alternate_bytes
        .len()
        .checked_sub(suffix)
        .ok_or_else(|| normalization_error("called variant alternate span underflow"))?;
    let edit_alternate = String::from_utf8(alternate_bytes[prefix..alternate_end].to_vec())
        .map_err(|_| normalization_error("called variant alternate allele is not ASCII DNA"))?;

    if edit_start == edit_end && edit_alternate.is_empty() {
        return Err(normalization_error(
            "called variant does not change the reference sequence",
        ));
    }

    Ok(SequenceEdit {
        start: edit_start,
        end: edit_end,
        alternate: edit_alternate,
    })
}

fn validate_kind(variant: &Variant) -> Result<()> {
    let valid = match variant.kind {
        VariantKind::Snv => variant.reference.len() == 1 && variant.alternate.len() == 1,
        VariantKind::Ins => variant.alternate.len() > variant.reference.len(),
        VariantKind::Del => variant.reference.len() > variant.alternate.len(),
    };
    if valid {
        Ok(())
    } else {
        Err(normalization_error(
            "called variant kind disagrees with its reference/alternate alleles",
        ))
    }
}

fn validate_allele(label: &str, allele: &str) -> Result<()> {
    if allele.is_empty()
        || !allele
            .bytes()
            .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T' | b'N'))
    {
        return Err(normalization_error(format!(
            "called variant {label} allele must contain uppercase DNA bases"
        )));
    }
    Ok(())
}

fn apply_edits(reference: &str, edits: &[SequenceEdit]) -> Result<String> {
    let mut ordered = edits.to_vec();
    ordered.sort_by(|left, right| {
        (left.start, left.end, &left.alternate).cmp(&(right.start, right.end, &right.alternate))
    });

    let mut output = String::with_capacity(reference.len());
    let mut cursor = 0;
    let mut insertion_anchors = BTreeSet::new();

    for edit in ordered {
        if edit.start > edit.end || edit.end > reference.len() {
            return Err(normalization_error(
                "variant edit lies outside the supplied reference",
            ));
        }
        if edit.start < cursor {
            return Err(normalization_error(
                "called variant edits overlap on the reference",
            ));
        }
        if edit.start == edit.end && !insertion_anchors.insert(edit.start) {
            return Err(normalization_error(
                "multiple insertion edits share one reference boundary",
            ));
        }

        output.push_str(&reference[cursor..edit.start]);
        output.push_str(&edit.alternate);
        cursor = edit.end;
    }
    output.push_str(&reference[cursor..]);
    Ok(output)
}

fn right_align_equivalent_indels(
    reference: &str,
    alternate_sequence: &str,
    source_edits: &[SequenceEdit],
) -> Result<Vec<SequenceEdit>> {
    let mut normalized = source_edits.to_vec();

    loop {
        let mut order = (0..normalized.len()).collect::<Vec<_>>();
        order.sort_by_key(|index| std::cmp::Reverse(normalized[*index].start));

        let mut moved = false;
        for index in order {
            let Some(shifted) = shift_right_once(reference, &normalized[index]) else {
                continue;
            };
            let mut candidate = normalized.clone();
            candidate[index] = shifted;
            if apply_edits(reference, &candidate)? == alternate_sequence {
                normalized = candidate;
                moved = true;
                break;
            }
        }
        if !moved {
            break;
        }
    }

    normalized.sort_by(|left, right| {
        (left.start, left.end, &left.alternate).cmp(&(right.start, right.end, &right.alternate))
    });
    Ok(normalized)
}

fn shift_right_once(reference: &str, edit: &SequenceEdit) -> Option<SequenceEdit> {
    let reference_bytes = reference.as_bytes();

    if edit.start == edit.end {
        let first_inserted = *edit.alternate.as_bytes().first()?;
        let crossed = *reference_bytes.get(edit.start)?;
        if first_inserted != crossed {
            return None;
        }
        let mut alternate = edit.alternate.as_bytes()[1..].to_vec();
        alternate.push(crossed);
        return Some(SequenceEdit {
            start: edit.start + 1,
            end: edit.end + 1,
            alternate: String::from_utf8(alternate).ok()?,
        });
    }

    if edit.alternate.is_empty() {
        let first_deleted = *reference_bytes.get(edit.start)?;
        let following = *reference_bytes.get(edit.end)?;
        if first_deleted != following {
            return None;
        }
        return Some(SequenceEdit {
            start: edit.start + 1,
            end: edit.end + 1,
            alternate: String::new(),
        });
    }

    None
}

fn render_edits(contig: &str, reference: &str, edits: &[SequenceEdit]) -> Result<Vec<Variant>> {
    let mut variants = edits
        .iter()
        .map(|edit| render_edit(contig, reference, edit))
        .collect::<Result<Vec<_>>>()?;
    variants.sort_by(|left, right| {
        (
            &left.contig,
            left.position_1based,
            &left.reference,
            &left.alternate,
        )
            .cmp(&(
                &right.contig,
                right.position_1based,
                &right.reference,
                &right.alternate,
            ))
    });
    Ok(variants)
}

fn render_edit(contig: &str, reference: &str, edit: &SequenceEdit) -> Result<Variant> {
    if edit.start == edit.end {
        if edit.start > 0 {
            let anchor = reference_base(reference, edit.start - 1)?;
            return Ok(Variant {
                contig: contig.to_owned(),
                position_1based: edit.start,
                reference: anchor.to_string(),
                alternate: format!("{anchor}{}", edit.alternate),
                kind: VariantKind::Ins,
            });
        }

        let anchor = reference_base(reference, 0)?;
        return Ok(Variant {
            contig: contig.to_owned(),
            position_1based: 1,
            reference: anchor.to_string(),
            alternate: format!("{}{anchor}", edit.alternate),
            kind: VariantKind::Ins,
        });
    }

    if edit.alternate.is_empty() {
        let deleted = reference
            .get(edit.start..edit.end)
            .ok_or_else(|| normalization_error("normalized deletion lies outside the reference"))?;
        if edit.start > 0 {
            let anchor = reference_base(reference, edit.start - 1)?;
            return Ok(Variant {
                contig: contig.to_owned(),
                position_1based: edit.start,
                reference: format!("{anchor}{deleted}"),
                alternate: anchor.to_string(),
                kind: VariantKind::Del,
            });
        }

        let anchor = reference_base(reference, edit.end)?;
        return Ok(Variant {
            contig: contig.to_owned(),
            position_1based: 1,
            reference: format!("{deleted}{anchor}"),
            alternate: anchor.to_string(),
            kind: VariantKind::Del,
        });
    }

    if edit.end - edit.start == 1 && edit.alternate.len() == 1 {
        return Ok(Variant {
            contig: contig.to_owned(),
            position_1based: edit.start + 1,
            reference: reference_base(reference, edit.start)?.to_string(),
            alternate: edit.alternate.clone(),
            kind: VariantKind::Snv,
        });
    }

    Err(normalization_error(
        "normalization produced an unsupported replacement edit",
    ))
}

fn reference_base(reference: &str, index: usize) -> Result<char> {
    reference
        .as_bytes()
        .get(index)
        .copied()
        .map(char::from)
        .ok_or_else(|| normalization_error("normalization anchor lies outside the reference"))
}

fn normalization_error(message: impl Into<String>) -> Error {
    Error::VariantNormalization(message.into())
}
