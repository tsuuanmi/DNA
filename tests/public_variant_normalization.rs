//! Public `variant_normalization` capability and its policies.

mod support;

use std::path::Path;

use sha2::{Digest, Sha256};
use tempfile::tempdir;

use dna::variant_analysis::{CalledVariantSet, ReferenceIdentity, Variant, VariantKind};
use dna::variant_normalization::{self, NormalizationPolicy};
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
) -> Result<variant_normalization::VariantNormalizationResult, dna::error::Error> {
    variant_normalization::normalize(
        reference_path,
        called,
        NormalizationPolicy::MtDnaRightAligned,
    )
}

#[test]
fn right_alignment_leaves_snv_unchanged() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let reference_path = directory.path().join("reference.fa");
    let sequence = "CAAAAG";
    write_reference(&reference_path, sequence)?;
    let source = variant(6, "G", "T", VariantKind::Snv);
    let called = called(sequence, vec![source.clone()]);

    let result = normalize(&reference_path, &called)?;

    assert_eq!(result.source_variants, vec![source.clone()]);
    assert_eq!(result.normalized_variants, vec![source]);
    Ok(())
}

#[test]
fn right_aligns_homopolymer_deletion_without_changing_haplotype()
-> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let reference_path = directory.path().join("reference.fa");
    let sequence = "CAAAAG";
    write_reference(&reference_path, sequence)?;
    let source = variant(1, "CA", "C", VariantKind::Del);
    let called = called(sequence, vec![source.clone()]);

    let result = normalize(&reference_path, &called)?;

    assert_eq!(result.source_variants, vec![source]);
    assert_eq!(
        result.normalized_variants,
        vec![variant(4, "AA", "A", VariantKind::Del)]
    );
    Ok(())
}

#[test]
fn right_aligns_homopolymer_insertion_without_changing_haplotype()
-> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let reference_path = directory.path().join("reference.fa");
    let sequence = "CAAAAG";
    write_reference(&reference_path, sequence)?;
    let source = variant(1, "C", "CA", VariantKind::Ins);
    let called = called(sequence, vec![source.clone()]);

    let result = normalize(&reference_path, &called)?;

    assert_eq!(result.source_variants, vec![source]);
    assert_eq!(
        result.normalized_variants,
        vec![variant(5, "A", "AA", VariantKind::Ins)]
    );
    Ok(())
}

#[test]
fn right_aligns_tandem_repeat_deletion() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let reference_path = directory.path().join("reference.fa");
    let sequence = "CATATG";
    write_reference(&reference_path, sequence)?;
    let source = variant(1, "CAT", "C", VariantKind::Del);
    let called = called(sequence, vec![source.clone()]);

    let result = normalize(&reference_path, &called)?;

    assert_eq!(result.source_variants, vec![source]);
    assert_eq!(result.alternate_sequence, "CATG");
    assert_eq!(
        result.normalized_variants,
        vec![variant(3, "TAT", "T", VariantKind::Del)]
    );
    Ok(())
}

#[test]
fn right_aligns_tandem_repeat_insertion() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let reference_path = directory.path().join("reference.fa");
    let sequence = "CATATG";
    write_reference(&reference_path, sequence)?;
    let source = variant(1, "C", "CAT", VariantKind::Ins);
    let called = called(sequence, vec![source.clone()]);

    let result = normalize(&reference_path, &called)?;

    assert_eq!(result.source_variants, vec![source]);
    assert_eq!(result.alternate_sequence, "CATATATG");
    assert_eq!(
        result.normalized_variants,
        vec![variant(5, "T", "TAT", VariantKind::Ins)]
    );
    Ok(())
}

#[test]
fn normalization_preserves_phase_across_nearby_edits() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let reference_path = directory.path().join("reference.fa");
    let sequence = "CAAAAG";
    write_reference(&reference_path, sequence)?;
    let insertion = variant(1, "C", "CA", VariantKind::Ins);
    let snv = variant(5, "A", "T", VariantKind::Snv);
    let called = called(sequence, vec![insertion.clone(), snv.clone()]);

    let result = normalize(&reference_path, &called)?;

    assert_eq!(result.source_variants, vec![insertion, snv]);
    assert_eq!(
        result.normalized_variants,
        vec![
            variant(4, "A", "AA", VariantKind::Ins),
            variant(5, "A", "T", VariantKind::Snv),
        ]
    );
    Ok(())
}

#[test]
fn right_alignment_stops_at_the_mtdna_coordinate_seam() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let reference_path = directory.path().join("reference.fa");
    let sequence = "AAAA";
    write_reference(&reference_path, sequence)?;
    let source = variant(1, "AA", "A", VariantKind::Del);
    let called = called(sequence, vec![source]);

    let result = normalize(&reference_path, &called)?;

    assert_eq!(
        result.normalized_variants,
        vec![variant(3, "AA", "A", VariantKind::Del)]
    );
    Ok(())
}

#[test]
fn rejects_reference_identity_mismatch() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let reference_path = directory.path().join("reference.fa");
    let sequence = "CAAAAG";
    write_reference(&reference_path, sequence)?;
    let mut called = called(sequence, vec![variant(6, "G", "T", VariantKind::Snv)]);
    called.reference.sha256 = "wrong".into();

    let Err(error) = normalize(&reference_path, &called) else {
        return Err("expected reference identity mismatch".into());
    };

    assert!(error.to_string().contains("reference identity"));
    Ok(())
}

#[test]
fn rejects_source_allele_that_disagrees_with_reference() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let reference_path = directory.path().join("reference.fa");
    let sequence = "CAAAAG";
    write_reference(&reference_path, sequence)?;
    let called = called(sequence, vec![variant(2, "G", "T", VariantKind::Snv)]);

    let Err(error) = normalize(&reference_path, &called) else {
        return Err("expected reference allele mismatch".into());
    };

    assert!(error.to_string().contains("reference allele"));
    Ok(())
}
