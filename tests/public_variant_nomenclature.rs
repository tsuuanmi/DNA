mod support;

use sha2::{Digest, Sha256};
use tempfile::tempdir;

use dna::variant_analysis::{
    CalledVariantSet, ReferenceIdentity, Variant, VariantKind,
};
use dna::variant_nomenclature;
use dna::variant_normalization::{self, NormalizationPolicy};
use support::write_reference;

fn sha256(sequence: &str) -> String {
    format!("{:x}", Sha256::digest(sequence.as_bytes()))
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

#[test]
fn normalization_result_exposes_complete_read_only_nomenclature_context()
-> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let reference_path = directory.path().join("reference.fa");
    let sequence = "CAAAAG";
    write_reference(&reference_path, sequence)?;

    let source = variant(1, "C", "CA", VariantKind::Ins);
    let called = CalledVariantSet {
        reference: ReferenceIdentity {
            name: "synthetic".into(),
            sha256: sha256(sequence),
        },
        variants: vec![source],
    };

    let normalized = variant_normalization::normalize(
        &reference_path,
        &called,
        NormalizationPolicy::MtDnaRightAligned,
    )?;
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
