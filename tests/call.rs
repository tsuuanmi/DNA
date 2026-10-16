//! CLI contract for the core-only call over reviewed consensus sequences.

pub mod support;

use std::fs;
use std::path::Path;

use assert_cmd::Command;
use predicates::prelude::*;
use serde_json::{Value, json};
use tempfile::tempdir;

use support::{write_config, write_reference};

const SAMPLE_ID: &str = "sample-1";

fn dna_binary() -> String {
    std::env::var("CARGO_BIN_EXE_dna")
        .unwrap_or_else(|error| panic!("Cargo must provide the DNA test binary: {error}"))
}

/// A deterministic reference without repeats long enough to confuse placement.
fn reference_sequence(length: usize) -> String {
    let mut state: u32 = 0x2545_f491;
    (0..length)
        .map(|_| {
            state = state.wrapping_mul(1_103_515_245).wrapping_add(12_345);
            ['A', 'C', 'G', 'T'][((state >> 16) & 3) as usize]
        })
        .collect()
}

fn substitute(sequence: &str, index: usize) -> String {
    let mut bases = sequence.chars().collect::<Vec<_>>();
    bases[index] = if bases[index] == 'A' { 'C' } else { 'A' };
    bases.into_iter().collect()
}

/// Writes a strict configuration whose read-end margin is `margin`.
fn write_config_with_margin(path: &Path, margin: usize) -> Result<(), Box<dyn std::error::Error>> {
    write_config(path, "linear")?;
    let text = fs::read_to_string(path)?;
    fs::write(
        path,
        text.replace("read_end_margin=0", &format!("read_end_margin={margin}")),
    )?;
    Ok(())
}

fn run(
    sequences: &Path,
    reference: &Path,
    config: &Path,
    workdir: &Path,
) -> assert_cmd::assert::Assert {
    Command::new(dna_binary())
        .current_dir(workdir)
        .env("DNA_CONFIG", config)
        .env("DNA_LOG_DIR", workdir.join("logs"))
        .arg("call")
        .arg(SAMPLE_ID)
        .arg(sequences)
        .arg("--reference")
        .arg(reference)
        .assert()
}

#[test]
fn calls_variants_from_consensus_sequences_with_the_core_alone()
-> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let sequences = directory.path().join("consensus.fasta");
    let sequence = reference_sequence(160);
    write_reference(&reference, &sequence)?;
    write_config_with_margin(&config, 12)?;

    // The first record changes its very first base and one inner base; the
    // second inserts one base. Reviewed ends are vouched for, so the first
    // base is callable despite the read-end margin.
    let first = substitute(&sequence[..70], 0);
    let first = substitute(&first, 35);
    let second = format!("{}A{}", &sequence[80..120], &sequence[120..160]);
    fs::write(
        &sequences,
        format!(">HV1\n{first}\n>HV2 reviewed\n{}\n", second.to_lowercase()),
    )?;

    run(&sequences, &reference, &config, directory.path())
        .success()
        .stdout(predicate::str::is_empty());

    let value: Value = serde_json::from_slice(&fs::read(
        directory
            .path()
            .join("results")
            .join(format!("{SAMPLE_ID}.variants.json")),
    )?)?;
    assert_eq!(value["schema_version"], "dna.variants/v1");
    assert_eq!(value["sample_id"], SAMPLE_ID);
    assert_eq!(
        value["provenance"]["plugins"],
        json!([
            {"id": "sequence", "family": "modality", "version": 1},
            {"id": "core", "family": "core", "version": 4},
        ])
    );
    assert!(value.get("notation").is_none());
    let reads = value["reads"].as_array().ok_or("reads must be an array")?;
    assert_eq!(reads.len(), 2);
    assert_eq!(reads[0]["name"], "HV1");
    assert_eq!(reads[0]["alignment"]["orientation"], "forward");
    assert_eq!(
        reads[0]["alignment"]["reference_segments"],
        json!([{"start": 0, "end": 70}])
    );
    let positions = |read: &Value| -> Vec<(u64, bool)> {
        read["variants"]
            .as_array()
            .map(|variants| {
                variants
                    .iter()
                    .filter_map(|variant| {
                        Some((
                            variant["position"].as_u64()?,
                            variant["eligible"].as_bool()?,
                        ))
                    })
                    .collect()
            })
            .unwrap_or_default()
    };
    assert_eq!(positions(&reads[0]), [(1, true), (36, true)]);
    assert_eq!(reads[1]["name"], "HV2");
    assert_eq!(reads[1]["variants"][0]["kind"], "INS");
    assert_eq!(reads[1]["variants"][0]["eligible"], true);
    assert_eq!(reads[1]["variants"].as_array().map(Vec::len), Some(1));
    Ok(())
}

#[test]
fn rejects_sequences_with_unsupported_symbols() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let sequences = directory.path().join("consensus.fasta");
    let sequence = reference_sequence(120);
    write_reference(&reference, &sequence)?;
    write_config_with_margin(&config, 0)?;
    fs::write(
        &sequences,
        format!(">HV1\n{}-{}\n", &sequence[..40], &sequence[41..80]),
    )?;

    run(&sequences, &reference, &config, directory.path())
        .failure()
        .stderr(predicate::str::contains(
            "invalid sequence FASTA: unsupported sequence symbol '-' in \"HV1\"",
        ));
    assert!(!directory.path().join("results").exists());
    Ok(())
}

#[test]
fn unresolved_bases_split_the_callable_reference_segments() -> Result<(), Box<dyn std::error::Error>>
{
    let directory = tempdir()?;
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let sequences = directory.path().join("consensus.fasta");
    let sequence = reference_sequence(160);
    write_reference(&reference, &sequence)?;
    write_config_with_margin(&config, 12)?;

    // An `N` run over reference 40..48, an ambiguity code at 60, and an
    // inner substitution two bases after the `N` run.
    let mut bases = substitute(&sequence[..100], 50).chars().collect::<Vec<_>>();
    bases[40..48].fill('N');
    bases[60] = 'y';
    let record = bases.into_iter().collect::<String>();
    fs::write(&sequences, format!(">HV1\n{record}\n"))?;

    run(&sequences, &reference, &config, directory.path()).success();

    let value: Value = serde_json::from_slice(&fs::read(
        directory
            .path()
            .join("results")
            .join(format!("{SAMPLE_ID}.variants.json")),
    )?)?;
    let alignment = &value["reads"][0]["alignment"];
    assert_eq!(
        alignment["reference_segments"],
        json!([{"start": 0, "end": 100}])
    );
    assert_eq!(
        alignment["callable_reference_segments"],
        json!([
            {"start": 0, "end": 40},
            {"start": 48, "end": 60},
            {"start": 61, "end": 100},
        ])
    );
    assert_eq!(alignment["unresolved_bases"], 9);
    assert_eq!(alignment["masked_bases"], 0);
    // An unresolved base masks no call: the substitution beside the `N` run
    // stays eligible, and the ambiguity code is no variant.
    let variants = value["reads"][0]["variants"]
        .as_array()
        .ok_or("variants must be an array")?;
    assert_eq!(variants.len(), 1);
    assert_eq!(variants[0]["position"], 51);
    assert_eq!(variants[0]["eligible"], true);
    Ok(())
}

#[test]
fn rejects_a_sequence_without_a_resolved_base() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let sequences = directory.path().join("consensus.fasta");
    write_reference(&reference, &reference_sequence(120))?;
    write_config_with_margin(&config, 0)?;
    // Nothing places a read without a resolved base: it fails typed.
    fs::write(&sequences, format!(">HV1\n{}\n", "N".repeat(60)))?;

    run(&sequences, &reference, &config, directory.path())
        .failure()
        .stderr(predicate::str::contains("error: alignment failed:"));
    assert!(!directory.path().join("results").exists());
    Ok(())
}
