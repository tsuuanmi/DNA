//! Public `variant_normalization` capability and its policies.

pub mod support;

use std::path::Path;

use sha2::{Digest, Sha256};
use tempfile::tempdir;

use dna::error::{KernelError, NormalizationError, RepresentationError};
use dna::variant::{CalledVariantSet, ReferenceIdentity, Variant, VariantKind};
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
) -> Result<variant_normalization::VariantNormalizationResult, KernelError> {
    variant_normalization::normalize(
        reference_path,
        called,
        NormalizationPolicy::RightAligned { max_trade: 20 },
    )
}

// SRS-VN-005 SRS-VN-008 SRS-API-007
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

// SRS-VN-006
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

// SRS-VN-006
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

// SRS-VN-006
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

// SRS-VN-006
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

// SRS-VN-005 SRS-VN-007
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

// SRS-VN-006 SRS-VN-007
#[test]
fn an_insertion_crosses_a_substituted_base_of_its_own() -> Result<(), Box<dyn std::error::Error>> {
    // As 16172C with an inserted C at 16171: the sample's C run is 4-6 plus
    // the inserted C, so the insertion belongs after 6 (16174.1C), although
    // the reference base it crosses at 4 is T.
    let directory = tempdir()?;
    let reference_path = directory.path().join("reference.fa");
    let sequence = "GAATCCAG";
    write_reference(&reference_path, sequence)?;
    let insertion = variant(3, "A", "AC", VariantKind::Ins);
    let snv = variant(4, "T", "C", VariantKind::Snv);
    let called = called(sequence, vec![insertion, snv.clone()]);

    let result = normalize(&reference_path, &called)?;

    assert_eq!(
        result.normalized_variants,
        vec![snv, variant(6, "C", "CC", VariantKind::Ins)]
    );
    Ok(())
}

// SRS-VN-007 SRS-VN-008
#[test]
fn a_deletion_moves_3_prime_across_substitutions_it_replaces()
-> Result<(), Box<dyn std::error::Error>> {
    // As 16293DEL 16297C 16298C over CCCACCCTT: the sample's C run is eight
    // bases, written with the deletion at its 3' end: 16293C 16297C 16298DEL.
    let directory = tempdir()?;
    let reference_path = directory.path().join("reference.fa");
    let sequence = "CCCACCCTTAA";
    write_reference(&reference_path, sequence)?;
    let called = called(
        sequence,
        vec![
            variant(3, "CA", "C", VariantKind::Del),
            variant(8, "T", "C", VariantKind::Snv),
            variant(9, "T", "C", VariantKind::Snv),
        ],
    );

    let result = normalize(&reference_path, &called)?;

    assert_eq!(
        result.normalized_variants,
        vec![
            variant(4, "A", "C", VariantKind::Snv),
            variant(8, "T", "C", VariantKind::Snv),
            variant(8, "TT", "T", VariantKind::Del),
        ]
    );
    Ok(())
}

// SRS-VN-006
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

// SRS-VN-002
#[test]
fn rejects_reference_identity_mismatch() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let reference_path = directory.path().join("reference.fa");
    let sequence = "CAAAAG";
    write_reference(&reference_path, sequence)?;
    let mut called = called(sequence, vec![variant(6, "G", "T", VariantKind::Snv)]);
    called.reference.sha256 = "wrong".into();

    assert!(matches!(
        normalize(&reference_path, &called),
        Err(KernelError::VariantNormalization(
            NormalizationError::ReferenceIdentityMismatch
        ))
    ));
    Ok(())
}

#[test]
fn rejects_source_allele_that_disagrees_with_reference() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let reference_path = directory.path().join("reference.fa");
    let sequence = "CAAAAG";
    write_reference(&reference_path, sequence)?;
    let called = called(sequence, vec![variant(2, "G", "T", VariantKind::Snv)]);

    assert!(matches!(
        normalize(&reference_path, &called),
        Err(KernelError::VariantNormalization(
            NormalizationError::Representation(RepresentationError::ReferenceAlleleMismatch)
        ))
    ));
    Ok(())
}

// SRS-VN-003 SRS-VN-004
#[test]
fn source_variants_are_applied_jointly_into_one_alternate_haplotype()
-> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let reference_path = directory.path().join("reference.fa");
    let sequence = "CAAAAG";
    write_reference(&reference_path, sequence)?;
    let called = called(
        sequence,
        vec![
            variant(1, "C", "CA", VariantKind::Ins),
            variant(5, "A", "T", VariantKind::Snv),
        ],
    );

    let result = normalize(&reference_path, &called)?;

    // C, inserted A, A2..A4, substituted T, G: both edits in one haplotype.
    assert_eq!(result.alternate_sequence, "CAAAATG");
    // The moved insertion and the substitution still give the same haplotype.
    let moved = variant(4, "A", "AA", VariantKind::Ins);
    assert!(result.normalized_variants.contains(&moved));
    assert_eq!(result.normalized_variants.len(), 2);
    Ok(())
}

// SRS-VN-009
#[test]
fn normalized_variants_are_emitted_in_position_order() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let reference_path = directory.path().join("reference.fa");
    let sequence = "CAAAAGTTCC";
    write_reference(&reference_path, sequence)?;
    let late = variant(9, "C", "A", VariantKind::Snv);
    let early = variant(6, "G", "T", VariantKind::Snv);
    let called = called(sequence, vec![late.clone(), early.clone()]);

    let result = normalize(&reference_path, &called)?;

    assert_eq!(result.normalized_variants, vec![early, late]);
    Ok(())
}
