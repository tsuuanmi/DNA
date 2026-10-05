//! Shared crate-internal mechanics for sequence-edit and anchored variant representations.

use std::collections::BTreeSet;

use crate::variant_analysis::{Variant, VariantKind};

pub(crate) type RepresentationResult<T> = std::result::Result<T, String>;

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
) -> RepresentationResult<Vec<SequenceEdit>> {
    variants
        .iter()
        .map(|variant| variant_to_edit(contig, reference, variant))
        .collect()
}

pub(crate) fn apply_edits(reference: &str, edits: &[SequenceEdit]) -> RepresentationResult<String> {
    let mut ordered = edits.to_vec();
    sort_edits(&mut ordered);

    let mut output = String::with_capacity(reference.len());
    let mut cursor = 0;
    let mut insertion_anchors = BTreeSet::new();

    for edit in ordered {
        if edit.start > edit.end || edit.end > reference.len() {
            return Err("variant edit lies outside the supplied reference".into());
        }
        if edit.start < cursor {
            return Err("called variant edits overlap on the reference".into());
        }
        if edit.start == edit.end && !insertion_anchors.insert(edit.start) {
            return Err("multiple insertion edits share one reference boundary".into());
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

pub(crate) fn render_edits(
    contig: &str,
    reference: &str,
    edits: &[SequenceEdit],
) -> RepresentationResult<Vec<Variant>> {
    let mut variants = edits
        .iter()
        .map(|edit| render_edit(contig, reference, edit))
        .collect::<RepresentationResult<Vec<_>>>()?;
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

fn variant_to_edit(
    contig: &str,
    reference: &str,
    variant: &Variant,
) -> RepresentationResult<SequenceEdit> {
    if variant.contig != contig {
        return Err("called variant contig does not match the supplied reference".into());
    }
    if variant.position_1based == 0 {
        return Err("called variant position must be one-based".into());
    }
    validate_allele("reference", &variant.reference)?;
    validate_allele("alternate", &variant.alternate)?;
    validate_kind(variant)?;

    let start = variant.position_1based - 1;
    let allele_end = start
        .checked_add(variant.reference.len())
        .ok_or_else(|| "called variant reference span overflow".to_owned())?;
    if allele_end > reference.len() {
        return Err(
            "origin-spanning source variants are not normalized across the canonical seam".into(),
        );
    }
    if reference.as_bytes().get(start..allele_end) != Some(variant.reference.as_bytes()) {
        return Err("called variant reference allele disagrees with the supplied reference".into());
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
        .ok_or_else(|| "called variant edit start overflow".to_owned())?;
    let edit_end = allele_end
        .checked_sub(suffix)
        .ok_or_else(|| "called variant edit end underflow".to_owned())?;
    let alternate_end = alternate_bytes
        .len()
        .checked_sub(suffix)
        .ok_or_else(|| "called variant alternate span underflow".to_owned())?;
    let edit_alternate = String::from_utf8(alternate_bytes[prefix..alternate_end].to_vec())
        .map_err(|_| "called variant alternate allele is not ASCII DNA".to_owned())?;

    if edit_start == edit_end && edit_alternate.is_empty() {
        return Err("called variant does not change the reference sequence".into());
    }

    Ok(SequenceEdit {
        start: edit_start,
        end: edit_end,
        alternate: edit_alternate,
    })
}

fn validate_kind(variant: &Variant) -> RepresentationResult<()> {
    let valid = match variant.kind {
        VariantKind::Snv => variant.reference.len() == 1 && variant.alternate.len() == 1,
        VariantKind::Ins => variant.alternate.len() > variant.reference.len(),
        VariantKind::Del => variant.reference.len() > variant.alternate.len(),
    };
    if valid {
        Ok(())
    } else {
        Err("called variant kind disagrees with its reference/alternate alleles".into())
    }
}

fn validate_allele(label: &str, allele: &str) -> RepresentationResult<()> {
    if allele.is_empty()
        || !allele
            .bytes()
            .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T' | b'N'))
    {
        return Err(format!(
            "called variant {label} allele must contain uppercase DNA bases"
        ));
    }
    Ok(())
}

fn render_edit(
    contig: &str,
    reference: &str,
    edit: &SequenceEdit,
) -> RepresentationResult<Variant> {
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
            .ok_or_else(|| "normalized deletion lies outside the reference".to_owned())?;
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

    Err("normalization produced an unsupported replacement edit".into())
}

fn reference_base(reference: &str, index: usize) -> RepresentationResult<char> {
    reference
        .as_bytes()
        .get(index)
        .copied()
        .map(char::from)
        .ok_or_else(|| "normalization anchor lies outside the reference".to_owned())
}
