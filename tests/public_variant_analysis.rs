//! Public `variant_analysis` capability without CLI side effects.

pub mod support;

use std::fs;

use tempfile::tempdir;

use dna::variant::VariantKind;
use dna::variant_analysis;
use support::{write_abif, write_config, write_reference};

const QUERY: &str = "ACGTCAGTACGATCGTACCTGAGTACGA";

#[test]
fn sanger_analysis_returns_canonical_variants_without_cli_side_effects()
-> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");

    write_abif(&trace, QUERY)?;
    let mut reference_query = QUERY.as_bytes().to_vec();
    reference_query[10] = if reference_query[10] == b'A' {
        b'C'
    } else {
        b'A'
    };
    write_reference(
        &reference,
        &format!("TTTT{}CCCC", String::from_utf8(reference_query)?),
    )?;
    write_config(&config, "linear")?;

    let result = variant_analysis::analyze_sanger(&trace, &reference, &config)?;

    assert_eq!(result.reference.name, "synthetic");
    assert_eq!(result.variants.len(), 1);
    let variant = &result.variants[0];
    assert_eq!(variant.kind, VariantKind::Snv);
    assert_eq!(variant.position_1based, 15);
    assert_eq!(variant.alternate.len(), 1);
    assert_ne!(variant.reference, variant.alternate);
    assert!(!result.input_sha256.is_empty());
    assert!(!result.reference.sha256.is_empty());
    assert!(!result.configuration_sha256.is_empty());
    assert!(!result.reference_segments.is_empty());

    let called = result.called_variants();
    assert_eq!(called.reference, result.reference);
    assert_eq!(called.variants, result.variants);

    assert!(!directory.path().join("results").exists());
    assert!(!directory.path().join("logs").exists());
    assert!(!fs::exists(directory.path().join("trace.log"))?);

    Ok(())
}

#[test]
fn sanger_analysis_accepts_abif_filename_without_extension_coupling()
-> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let trace = directory.path().join("trace.abif");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");

    write_abif(&trace, QUERY)?;
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;

    let result = variant_analysis::analyze_sanger(&trace, &reference, &config)?;

    assert_eq!(result.reference.name, "synthetic");
    assert!(result.variants.is_empty());
    assert!(!result.input_sha256.is_empty());

    Ok(())
}
