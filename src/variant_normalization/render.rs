//! Minimal sequence edits back to anchored public variant representations.

use crate::error::Result;
use crate::variant_analysis::{Variant, VariantKind};

use super::edit::SequenceEdit;
use super::normalization_error;

pub(super) fn render_edits(
    contig: &str,
    reference: &str,
    edits: &[SequenceEdit],
) -> Result<Vec<Variant>> {
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
