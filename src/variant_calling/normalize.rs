//! Minimal indel representation preserving canonical alignment placement.

use crate::error::{Error, Result};
use crate::model::coordinate::reference_one_based;
use crate::model::reference::{Reference, ReferenceTopology};
use crate::model::variant::{Variant, VariantCallMapping, VariantKind};
use crate::variant_calling::mapping;

pub(crate) fn snv(
    reference: &Reference,
    position: usize,
    alternate: char,
    calls: Vec<VariantCallMapping>,
) -> Result<Variant> {
    mapping::validate_snv(&calls, position)?;
    let reference_base = reference
        .sequence
        .as_bytes()
        .get(position)
        .copied()
        .ok_or_else(|| Error::Variant("SNV reference position is out of bounds".into()))?;
    validated(
        reference,
        Variant {
            contig: reference.name.clone(),
            position_1based: reference_one_based(position)?,
            reference: char::from(reference_base).to_string(),
            alternate: alternate.to_string(),
            kind: VariantKind::Snv,
            calls,
        },
    )
}

pub(crate) fn insertion(
    reference: &Reference,
    previous_reference: Option<usize>,
    next_reference: Option<usize>,
    inserted: String,
    calls: Vec<VariantCallMapping>,
) -> Result<Variant> {
    if inserted.is_empty() {
        return Err(Error::Variant("insertion allele is empty".into()));
    }
    mapping::validate_insertion(&calls)?;
    let anchor = observed_anchor(reference, previous_reference, next_reference)?;
    if let Some(anchor) = anchor {
        return build_insertion(reference, anchor, &inserted, calls);
    }

    let right_anchor = next_reference
        .ok_or_else(|| Error::Variant("leading insertion lacks a right anchor".into()))?;
    let base = reference_base(reference, right_anchor)?;
    validated(
        reference,
        Variant {
            contig: reference.name.clone(),
            position_1based: reference_one_based(right_anchor)?,
            reference: base.to_string(),
            alternate: format!("{inserted}{base}"),
            kind: VariantKind::Ins,
            calls,
        },
    )
}

pub(crate) fn deletion(
    reference: &Reference,
    previous_reference: Option<usize>,
    first_deleted_reference: usize,
    next_reference: Option<usize>,
    deleted: String,
    calls: Vec<VariantCallMapping>,
) -> Result<Variant> {
    if deleted.is_empty() {
        return Err(Error::Variant("deletion allele is empty".into()));
    }
    mapping::validate_deletion(&calls)?;
    let anchor = observed_anchor(reference, previous_reference, Some(first_deleted_reference))?;
    if let Some(anchor) = anchor {
        return build_deletion(reference, anchor, &deleted, calls);
    }

    let right_anchor = next_reference
        .ok_or_else(|| Error::Variant("leading deletion lacks a right anchor".into()))?;
    let base = reference_base(reference, right_anchor)?;
    validated(
        reference,
        Variant {
            contig: reference.name.clone(),
            position_1based: 1,
            reference: format!("{deleted}{base}"),
            alternate: base.to_string(),
            kind: VariantKind::Del,
            calls,
        },
    )
}

fn observed_anchor(
    reference: &Reference,
    previous_reference: Option<usize>,
    event_reference: Option<usize>,
) -> Result<Option<usize>> {
    if let Some(anchor) = previous_reference {
        reference_base(reference, anchor)?;
        return Ok(Some(anchor));
    }
    let Some(position) = event_reference else {
        return Ok(None);
    };
    reference_base(reference, position)?;
    Ok(match reference.topology {
        ReferenceTopology::Linear => position.checked_sub(1),
        ReferenceTopology::Circular => Some((position + reference.len() - 1) % reference.len()),
    })
}

fn build_insertion(
    reference: &Reference,
    anchor: usize,
    inserted: &str,
    calls: Vec<VariantCallMapping>,
) -> Result<Variant> {
    let base = reference_base(reference, anchor)?;
    validated(
        reference,
        Variant {
            contig: reference.name.clone(),
            position_1based: reference_one_based(anchor)?,
            reference: base.to_string(),
            alternate: format!("{base}{inserted}"),
            kind: VariantKind::Ins,
            calls,
        },
    )
}

fn build_deletion(
    reference: &Reference,
    anchor: usize,
    deleted: &str,
    calls: Vec<VariantCallMapping>,
) -> Result<Variant> {
    let base = reference_base(reference, anchor)?;
    validated(
        reference,
        Variant {
            contig: reference.name.clone(),
            position_1based: reference_one_based(anchor)?,
            reference: format!("{base}{deleted}"),
            alternate: base.to_string(),
            kind: VariantKind::Del,
            calls,
        },
    )
}

fn validated(reference: &Reference, variant: Variant) -> Result<Variant> {
    let start = variant
        .position_1based
        .checked_sub(1)
        .ok_or_else(|| Error::Variant("variant position must be one-based".into()))?;
    for (offset, observed) in variant.reference.bytes().enumerate() {
        let unwrapped = start
            .checked_add(offset)
            .ok_or_else(|| Error::Variant("variant reference span overflow".into()))?;
        let position = match reference.topology {
            ReferenceTopology::Linear => unwrapped,
            ReferenceTopology::Circular => unwrapped % reference.len(),
        };
        let expected = reference
            .sequence
            .as_bytes()
            .get(position)
            .copied()
            .ok_or_else(|| Error::Variant("variant reference span is out of bounds".into()))?;
        if observed != expected {
            return Err(Error::Variant(format!(
                "variant reference allele disagrees with the supplied reference at position {}",
                position + 1
            )));
        }
    }
    Ok(variant)
}

fn reference_base(reference: &Reference, position: usize) -> Result<char> {
    reference
        .sequence
        .as_bytes()
        .get(position)
        .copied()
        .map(char::from)
        .ok_or_else(|| Error::Variant("indel anchor is outside the reference".into()))
}

#[cfg(test)]
mod tests {
    use crate::model::variant::VariantCallRole;

    use super::*;

    fn reference(sequence: &str, topology: ReferenceTopology) -> Reference {
        Reference {
            name: "ref".into(),
            sequence: sequence.into(),
            topology,
            sequence_sha256: String::new(),
        }
    }

    fn inserted_calls() -> Vec<VariantCallMapping> {
        vec![VariantCallMapping {
            role: VariantCallRole::Supporting,
            call_index_0based: 0,
            reference_position_0based: None,
        }]
    }

    fn deletion_flanks() -> Vec<VariantCallMapping> {
        vec![VariantCallMapping {
            role: VariantCallRole::Flanking,
            call_index_0based: 0,
            reference_position_0based: Some(0),
        }]
    }

    #[test]
    fn preserves_canonical_homopolymer_insertion_anchor() -> Result<()> {
        let variant = insertion(
            &reference("CAAAAG", ReferenceTopology::Linear),
            Some(4),
            Some(5),
            "A".into(),
            inserted_calls(),
        )?;
        assert_eq!(variant.position_1based, 5);
        assert_eq!(variant.reference, "A");
        assert_eq!(variant.alternate, "AA");
        Ok(())
    }

    #[test]
    fn normalization_preserves_observed_call_positions() -> Result<()> {
        let calls = vec![
            VariantCallMapping {
                role: VariantCallRole::Flanking,
                call_index_0based: 4,
                reference_position_0based: Some(4),
            },
            VariantCallMapping {
                role: VariantCallRole::Flanking,
                call_index_0based: 5,
                reference_position_0based: Some(6),
            },
        ];
        let variant = deletion(
            &reference("CAAAAAG", ReferenceTopology::Linear),
            Some(4),
            5,
            Some(6),
            "A".into(),
            calls.clone(),
        )?;
        assert_eq!(variant.position_1based, 5);
        assert_eq!(variant.reference, "AA");
        assert_eq!(variant.alternate, "A");
        assert_eq!(variant.calls, calls);
        Ok(())
    }

    #[test]
    fn circular_representation_preserves_origin_seam_anchor() -> Result<()> {
        let variant = deletion(
            &reference("AAAA", ReferenceTopology::Circular),
            Some(3),
            0,
            Some(1),
            "A".into(),
            deletion_flanks(),
        )?;
        assert_eq!(variant.position_1based, 4);
        assert_eq!(variant.reference, "AA");
        assert_eq!(variant.alternate, "A");
        Ok(())
    }

    #[test]
    fn circular_insertion_preserves_observed_seam_anchor() -> Result<()> {
        let reference = reference("AACAA", ReferenceTopology::Circular);
        let variant = insertion(&reference, Some(4), Some(0), "A".into(), inserted_calls())?;
        assert_eq!(variant.position_1based, 5);
        assert_eq!(variant.reference, "A");
        assert_eq!(variant.alternate, "AA");
        Ok(())
    }

    #[test]
    fn rejects_reference_alleles_that_disagree_with_the_reference() {
        let result = deletion(
            &reference("TTCG", ReferenceTopology::Linear),
            Some(1),
            2,
            Some(3),
            "A".into(),
            deletion_flanks(),
        );
        assert!(matches!(result, Err(Error::Variant(message)) if message.contains("disagrees")));
    }

    #[test]
    fn derives_internal_linear_predecessor_without_an_aligned_left_flank() -> Result<()> {
        let reference = reference("TTCG", ReferenceTopology::Linear);
        let variant = deletion(&reference, None, 2, Some(3), "C".into(), deletion_flanks())?;
        assert_eq!(variant.position_1based, 2);
        assert_eq!(variant.reference, "TC");
        assert_eq!(variant.alternate, "T");
        Ok(())
    }

    #[test]
    fn derives_non_origin_circular_predecessor_without_an_aligned_left_flank() -> Result<()> {
        let reference = reference("ACGT", ReferenceTopology::Circular);
        let variant = deletion(&reference, None, 2, Some(3), "G".into(), deletion_flanks())?;
        assert_eq!(variant.position_1based, 2);
        assert_eq!(variant.reference, "CG");
        assert_eq!(variant.alternate, "C");
        Ok(())
    }
}
