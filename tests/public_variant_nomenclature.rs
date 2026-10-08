//! Public `variant_nomenclature` capability and the mtDNA HV2 poly-C rule.

mod support;

use std::path::Path;

use sha2::{Digest, Sha256};
use tempfile::tempdir;

use dna::variant_analysis::{CalledVariantSet, ReferenceIdentity, Variant, VariantKind};
use dna::variant_nomenclature;
use dna::variant_normalization::{self, NormalizationPolicy, VariantNormalizationResult};
use support::write_reference;

fn sha256(sequence: &str) -> String {
    format!("{:x}", Sha256::digest(sequence.as_bytes()))
}

fn called(sequence: &str, variants: Vec<Variant>) -> CalledVariantSet {
    CalledVariantSet {
        reference: ReferenceIdentity {
            name: "synthetic".into(),
            sha256: sha256(sequence),
        },
        variants,
    }
}

fn variant(position_1based: usize, reference: &str, alternate: &str, kind: VariantKind) -> Variant {
    Variant {
        contig: "synthetic".into(),
        position_1based,
        reference: reference.into(),
        alternate: alternate.into(),
        kind,
    }
}

fn normalize(
    reference_path: &Path,
    called: &CalledVariantSet,
) -> Result<VariantNormalizationResult, dna::error::Error> {
    variant_normalization::normalize(
        reference_path,
        called,
        NormalizationPolicy::MtDnaRightAligned,
    )
}

fn hv2_reference() -> String {
    format!("{}{}{}", "A".repeat(302), "CCCCCCCTCCCCC", "G".repeat(5))
}

#[test]
fn normalization_result_exposes_complete_read_only_nomenclature_context()
-> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let reference_path = directory.path().join("reference.fa");
    let sequence = "CAAAAG";
    write_reference(&reference_path, sequence)?;

    let source = variant(1, "C", "CA", VariantKind::Ins);
    let called = called(sequence, vec![source]);
    let normalized = normalize(&reference_path, &called)?;
    let input = variant_nomenclature::from_normalization(&normalized);

    assert_eq!(input.reference, &called.reference);
    assert_eq!(input.source_variants, normalized.source_variants.as_slice());
    assert_eq!(input.alternate_sequence, normalized.alternate_sequence);
    assert_eq!(
        input.normalized_variants,
        normalized.normalized_variants.as_slice()
    );
    Ok(())
}

#[test]
fn hv2_anchor_shift_is_represented_as_run_length_change() -> Result<(), Box<dyn std::error::Error>>
{
    let directory = tempdir()?;
    let reference_path = directory.path().join("reference.fa");
    let sequence = hv2_reference();
    write_reference(&reference_path, &sequence)?;

    let called = called(
        &sequence,
        vec![
            variant(309, "C", "T", VariantKind::Snv),
            variant(310, "T", "C", VariantKind::Snv),
        ],
    );
    let normalized = normalize(&reference_path, &called)?;
    assert_eq!(normalized.normalized_variants, called.variants);

    let input = variant_nomenclature::from_normalization(&normalized);
    let result = variant_nomenclature::mtdna::apply_hv2_polyc(&reference_path, input)?;

    assert_eq!(
        result.represented_variants,
        vec![
            variant(308, "CC", "C", VariantKind::Del),
            variant(315, "C", "CC", VariantKind::Ins),
        ]
    );
    assert_eq!(result.alternate_sequence, normalized.alternate_sequence);
    Ok(())
}

#[test]
fn hv2_multiple_c_insertions_keep_309_and_315_run_boundaries()
-> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let reference_path = directory.path().join("reference.fa");
    let sequence = hv2_reference();
    write_reference(&reference_path, &sequence)?;

    let called = called(
        &sequence,
        vec![
            variant(309, "C", "CCC", VariantKind::Ins),
            variant(315, "C", "CC", VariantKind::Ins),
        ],
    );
    let normalized = normalize(&reference_path, &called)?;
    let input = variant_nomenclature::from_normalization(&normalized);
    let result = variant_nomenclature::mtdna::apply_hv2_polyc(&reference_path, input)?;

    assert_eq!(
        result.represented_variants,
        vec![
            variant(309, "C", "CCC", VariantKind::Ins),
            variant(315, "C", "CC", VariantKind::Ins),
        ]
    );
    assert_eq!(result.alternate_sequence, normalized.alternate_sequence);
    Ok(())
}

#[test]
fn hv2_nomenclature_preserves_variants_outside_the_window() -> Result<(), Box<dyn std::error::Error>>
{
    let directory = tempdir()?;
    let reference_path = directory.path().join("reference.fa");
    let sequence = hv2_reference();
    write_reference(&reference_path, &sequence)?;

    let called = called(
        &sequence,
        vec![
            variant(200, "A", "G", VariantKind::Snv),
            variant(309, "C", "T", VariantKind::Snv),
            variant(310, "T", "C", VariantKind::Snv),
        ],
    );
    let normalized = normalize(&reference_path, &called)?;
    let input = variant_nomenclature::from_normalization(&normalized);
    let result = variant_nomenclature::mtdna::apply_hv2_polyc(&reference_path, input)?;

    assert_eq!(
        result.represented_variants,
        vec![
            variant(200, "A", "G", VariantKind::Snv),
            variant(308, "CC", "C", VariantKind::Del),
            variant(315, "C", "CC", VariantKind::Ins),
        ]
    );
    assert_eq!(result.alternate_sequence, normalized.alternate_sequence);
    Ok(())
}
