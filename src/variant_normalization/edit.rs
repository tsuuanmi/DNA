//! Source called variants to minimal sequence edits and haplotype reconstruction.

use std::collections::BTreeSet;

use crate::error::Result;
use crate::variant_analysis::{Variant, VariantKind};

use super::normalization_error;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct SequenceEdit {
    pub(crate) start: usize,
    pub(crate) end: usize,
    pub(crate) alternate: String,
}

pub(crate) fn variants_to_edits(
    contig: &str,
    reference: &str,
    variants: &[Variant],
) -> Result<Vec<SequenceEdit>> {
    variants
        .iter()
        .map(|variant| variant_to_edit(contig, reference, variant))
        .collect()
}

pub(crate) fn apply_edits(reference: &str, edits: &[SequenceEdit]) -> Result<String> {
    let mut ordered = edits.to_vec();
    sort_edits(&mut ordered);

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

pub(crate) fn sort_edits(edits: &mut [SequenceEdit]) {
    edits.sort_by(|left, right| {
        (left.start, left.end, &left.alternate).cmp(&(right.start, right.end, &right.alternate))
    });
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
