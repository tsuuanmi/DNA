//! Per-read profile representation for the notation view of sample, call,
//! and notation documents.
//!
//! ADR-0060 composes representation per read before sample reconciliation:
//! each read's eligible calls are normalized and then given the profile's
//! window representation, so sequence-equivalent descriptions from different
//! reads converge on one notation.

use crate::report::{ReadRepresentation, SampleNotation};
use dna_core::model::called_read::CalledRead;
use dna_kernel::error::{Error, NomenclatureError, Result};
use dna_kernel::model::reference::Reference;
use dna_kernel::profile::{Notation, Profile};
use dna_kernel::variant::{CalledVariantSet, ReferenceIdentity, Variant};
use dna_post::variant_nomenclature::{self, from_normalization};
use dna_post::variant_normalization::{NormalizationPolicy, normalize_with};

/// Represents every called read's eligible calls, or `None` when the profile
/// declares no notation. The reference has already been checked against the
/// profile.
pub(crate) fn represent(
    reads: &[CalledRead],
    reference: &Reference,
    profile: &Profile,
) -> Result<Option<SampleNotation>> {
    represent_variants(
        reads.iter().map(|read| {
            (
                read.input_sha256.clone(),
                read.variants.reported.iter().map(Variant::from).collect(),
            )
        }),
        reference,
        profile,
    )
}

/// Represents each read's eligible variants, given with its content identity.
pub(crate) fn represent_variants(
    reads: impl IntoIterator<Item = (String, Vec<Variant>)>,
    reference: &Reference,
    profile: &Profile,
) -> Result<Option<SampleNotation>> {
    let Some(Notation {
        indel_placement,
        style,
    }) = profile.notation
    else {
        return Ok(None);
    };
    let policy = NormalizationPolicy::for_placement(indel_placement);
    let reads = reads
        .into_iter()
        .map(|(input_sha256, variants)| {
            Ok(ReadRepresentation {
                input_sha256,
                variants: represent_read(reference, profile, policy, variants)?,
            })
        })
        .collect::<Result<_>>()?;
    Ok(Some(SampleNotation { style, reads }))
}

/// Normalizes one read's calls, then applies the profile windows. An edit that
/// straddles a window cannot take that representation (SRS-NOM-007), so the
/// read keeps its normalized form.
fn represent_read(
    reference: &Reference,
    profile: &Profile,
    policy: NormalizationPolicy,
    variants: Vec<Variant>,
) -> Result<Vec<Variant>> {
    let called = CalledVariantSet {
        reference: ReferenceIdentity {
            name: reference.name.clone(),
            sha256: reference.sequence_sha256.clone(),
        },
        variants,
    };
    let normalized = normalize_with(reference, &called, policy)?;
    match variant_nomenclature::apply_with(reference, profile, from_normalization(&normalized)) {
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

    use dna_core::model::variant::{Variant as CalledVariant, VariantKind};
    use dna_kernel::model::reference::ReferenceTopology;
    use dna_kernel::profile::fixtures::human_mtdna;
    use dna_kernel::reference;
    use dna_kernel::variant;

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
                Ordering::Less => variant::VariantKind::Ins,
                Ordering::Greater => variant::VariantKind::Del,
                Ordering::Equal => variant::VariantKind::Snv,
            },
        }
    }

    #[test]
    fn converges_a_left_run_c_insertion_on_309() -> TestResult {
        let reference = rcrs()?;
        let represented = represent_read(
            &reference,
            &human_mtdna()?,
            NormalizationPolicy::RightAligned,
            vec![Variant::from(&called(303, "C", "CC", VariantKind::Ins))],
        )?;
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
            NormalizationPolicy::RightAligned,
        )?;
        let profile = human_mtdna()?;
        assert!(matches!(
            variant_nomenclature::apply_with(&reference, &profile, from_normalization(&normalized)),
            Err(Error::VariantNomenclature(
                NomenclatureError::WindowCrossing { window }
            )) if window == "HVS-II"
        ));

        let represented = represent_read(
            &reference,
            &profile,
            NormalizationPolicy::RightAligned,
            vec![Variant::from(&called(301, "AAC", "A", VariantKind::Del))],
        )?;
        assert_eq!(represented, [crossing]);
        Ok(())
    }

    #[test]
    fn produces_notation_only_when_the_profile_declares_it() -> TestResult {
        let reference = rcrs()?;
        let mut profile = human_mtdna()?;
        assert!(represent(&[], &reference, &profile)?.is_some());
        profile.notation = None;
        assert!(represent(&[], &reference, &profile)?.is_none());
        Ok(())
    }
}
