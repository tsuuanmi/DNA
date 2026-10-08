//! Per-read human-mtDNA representation for the sample notation view.
//!
//! ADR-0060 composes representation per read before sample reconciliation:
//! each read's eligible calls are right-aligned and then given the
//! control-region representation, so sequence-equivalent descriptions from
//! different reads converge on one notation.

use crate::error::{Error, NomenclatureError, Result};
use crate::model::read_observation::ReadObservation;
use crate::model::reference::Reference;
use crate::model::variant::Variant as CalledVariant;
use crate::report::ReadRepresentation;
use crate::variant_analysis::{self, CalledVariantSet, ReferenceIdentity, Variant};
use crate::variant_nomenclature::from_normalization;
use crate::variant_nomenclature::mtdna::{control_region_with, is_rcrs};
use crate::variant_normalization::{NormalizationPolicy, normalize_with};

/// Represents every read's eligible calls, or `None` when the reference is not
/// the rCRS the human-mtDNA policies are validated against.
pub(crate) fn represent(
    reads: &[ReadObservation],
    reference: &Reference,
) -> Result<Option<Vec<ReadRepresentation>>> {
    if !is_rcrs(reference) {
        return Ok(None);
    }
    reads
        .iter()
        .map(|read| {
            Ok(ReadRepresentation {
                input_sha256: read.input_sha256.clone(),
                variants: represent_read(reference, &read.variants.reported)?,
            })
        })
        .collect::<Result<_>>()
        .map(Some)
}

/// Right-aligns one read's calls, then applies the control-region policy. An edit that
/// straddles a validated window cannot take that representation
/// (SRS-NOM-007), so the read keeps its right-aligned form.
fn represent_read(reference: &Reference, reported: &[CalledVariant]) -> Result<Vec<Variant>> {
    let called = CalledVariantSet {
        reference: ReferenceIdentity {
            name: reference.name.clone(),
            sha256: reference.sequence_sha256.clone(),
        },
        variants: reported
            .iter()
            .map(variant_analysis::project_variant)
            .collect(),
    };
    let normalized = normalize_with(reference, &called, NormalizationPolicy::MtDnaRightAligned)?;
    match control_region_with(reference, from_normalization(&normalized)) {
        Ok(represented) => Ok(represented.represented_variants),
        Err(Error::VariantNomenclature(NomenclatureError::WindowCrossing { .. })) => {
            Ok(normalized.normalized_variants)
        }
        Err(error) => Err(error),
    }
}

#[cfg(test)]
mod tests {
    use std::cmp::Ordering;
    use std::path::Path;

    use crate::model::reference::ReferenceTopology;
    use crate::model::variant::{Variant as CalledVariant, VariantKind};
    use crate::reference;
    use crate::variant_analysis;

    use super::*;

    type TestResult = std::result::Result<(), Box<dyn std::error::Error>>;

    fn rcrs() -> Result<Reference> {
        reference::load(
            &Path::new(env!("CARGO_MANIFEST_DIR")).join("references/rCRS.fasta"),
            ReferenceTopology::Circular,
        )
    }

    fn called(
        position: usize,
        reference: &str,
        alternate: &str,
        kind: VariantKind,
    ) -> CalledVariant {
        CalledVariant {
            contig: "rCRS".into(),
            position_1based: position,
            reference: reference.into(),
            alternate: alternate.into(),
            kind,
            calls: Vec::new(),
        }
    }

    fn public(position: usize, reference: &str, alternate: &str) -> Variant {
        Variant {
            contig: "rCRS".into(),
            position_1based: position,
            reference: reference.into(),
            alternate: alternate.into(),
            kind: match reference.len().cmp(&alternate.len()) {
                Ordering::Less => variant_analysis::VariantKind::Ins,
                Ordering::Greater => variant_analysis::VariantKind::Del,
                Ordering::Equal => variant_analysis::VariantKind::Snv,
            },
        }
    }

    #[test]
    fn converges_a_left_run_c_insertion_on_309() -> TestResult {
        let reference = rcrs()?;
        let represented = represent_read(&reference, &[called(303, "C", "CC", VariantKind::Ins)])?;
        assert_eq!(represented, [public(309, "C", "CC")]);
        Ok(())
    }

    #[test]
    fn keeps_the_normalized_form_for_an_edit_crossing_the_hv2_window() -> TestResult {
        let reference = rcrs()?;
        let crossing = public(301, "AAC", "A");
        let normalized = normalize_with(
            &reference,
            &CalledVariantSet {
                reference: ReferenceIdentity {
                    name: reference.name.clone(),
                    sha256: reference.sequence_sha256.clone(),
                },
                variants: vec![crossing.clone()],
            },
            NormalizationPolicy::MtDnaRightAligned,
        )?;
        assert!(matches!(
            control_region_with(&reference, from_normalization(&normalized)),
            Err(Error::VariantNomenclature(
                NomenclatureError::WindowCrossing { window: "HVS-II" }
            ))
        ));

        let represented = represent_read(&reference, &[called(301, "AAC", "A", VariantKind::Del)])?;
        assert_eq!(represented, [crossing]);
        Ok(())
    }

    #[test]
    fn applies_only_to_the_rcrs_reference() -> TestResult {
        let reference = rcrs()?;
        assert!(is_rcrs(&reference));
        let mut other = reference.clone();
        other.sequence_sha256 = "0".repeat(64);
        assert!(represent(&[], &other)?.is_none());
        assert!(represent(&[], &reference)?.is_some());
        Ok(())
    }
}
