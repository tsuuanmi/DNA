mod support;

use std::fs;

use tempfile::tempdir;

use dna::variant_analysis::{self, SangerAnalyzer};
use support::{write_abif, write_config, write_reference};

const QUERY: &str = "ACGTCAGTACGATCGTACCTGAGTACGA";

fn assert_send_sync<T: Send + Sync>() {}

#[test]
fn prepared_sanger_analyzer_is_send_sync() {
    assert_send_sync::<SangerAnalyzer>();
}

#[test]
fn prepared_sanger_analyzer_matches_one_shot_analysis() -> Result<(), Box<dyn std::error::Error>> {
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

    let analyzer = SangerAnalyzer::load(&reference, &config)?;
    let prepared = analyzer.analyze(&trace)?;
    let one_shot = variant_analysis::analyze_sanger(&trace, &reference, &config)?;

    assert_eq!(prepared, one_shot);
    assert_eq!(analyzer.reference_identity(), &prepared.reference);
    assert_eq!(
        analyzer.configuration_sha256(),
        prepared.configuration_sha256
    );

    assert!(!directory.path().join("results").exists());
    assert!(!directory.path().join("logs").exists());
    assert!(!fs::exists(directory.path().join("trace.log"))?);

    Ok(())
}

#[test]
fn one_prepared_analyzer_can_process_multiple_traces() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let trace_a = directory.path().join("a.ab1");
    let trace_b = directory.path().join("b.ab1");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");

    write_abif(&trace_a, QUERY)?;
    write_abif(&trace_b, QUERY)?;
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;

    let analyzer = SangerAnalyzer::load(&reference, &config)?;
    let left = analyzer.analyze(&trace_a)?;
    let right = analyzer.analyze(&trace_b)?;

    assert!(left.variants.is_empty());
    assert!(right.variants.is_empty());
    assert_eq!(left.reference, right.reference);
    assert_eq!(left.configuration_sha256, right.configuration_sha256);
    assert_ne!(left.input_sha256, "");
    assert_ne!(right.input_sha256, "");

    Ok(())
}
