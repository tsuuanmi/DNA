//! CLI contract for post-calling notation of a variants document.

pub mod support;

use std::fs;
use std::path::{Path, PathBuf};

use assert_cmd::Command;
use predicates::prelude::*;
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use tempfile::tempdir;

use support::{human_mtdna_profile, write_config_with_profile};

const SAMPLE_ID: &str = "sample-1";

fn dna_binary() -> String {
    std::env::var("CARGO_BIN_EXE_dna")
        .unwrap_or_else(|error| panic!("Cargo must provide the DNA test binary: {error}"))
}

fn rcrs_path() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("references/rCRS.fasta")
}

fn rcrs() -> Result<String, Box<dyn std::error::Error>> {
    Ok(fs::read_to_string(rcrs_path())?
        .lines()
        .filter(|line| !line.starts_with('>'))
        .collect())
}

fn dna(workdir: &Path, config: &Path, args: &[&str]) -> assert_cmd::assert::Assert {
    Command::new(dna_binary())
        .current_dir(workdir)
        .env("DNA_CONFIG", config)
        .env("DNA_LOG_DIR", workdir.join("logs"))
        .args(args)
        .arg("--reference")
        .arg(rcrs_path())
        .assert()
}

/// HVS-II 73-340 of a reviewed haplotype written `73G 309.1C 309.2C 314G
/// 315.1C` by EMPOP convention; DNA's window representation writes the right
/// C run as `313.1G` instead, with the same sequence.
fn hvs2(reference: &str) -> String {
    let mut region = String::new();
    for position in 73..=340 {
        region.push(match position {
            73 | 314 => 'G',
            _ => char::from(reference.as_bytes()[position - 1]),
        });
        match position {
            309 => region.push_str("CC"),
            315 => region.push('C'),
            _ => {}
        }
    }
    region
}

// SRS-NOM-016 SRS-NOM-018 SRS-IN-010 SRS-IN-014 SRS-OUT-001 SRS-OUT-002 SRS-OUT-011
#[test]
fn derives_the_in_process_notation_and_reports_conformance()
-> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let config = directory.path().join("dna.toml");
    write_config_with_profile(&config, &human_mtdna_profile())?;
    let sequences = directory.path().join("consensus.fasta");
    fs::write(&sequences, format!(">HV2\n{}\n", hvs2(&rcrs()?)))?;
    let sequences = sequences.to_string_lossy().into_owned();

    dna(directory.path(), &config, &["call", SAMPLE_ID, &sequences]).success();
    let variants_path = directory
        .path()
        .join("results")
        .join(format!("{SAMPLE_ID}.variants.json"));
    let variants_bytes = fs::read(&variants_path)?;
    let variants: Value = serde_json::from_slice(&variants_bytes)?;

    let document = variants_path.to_string_lossy().into_owned();
    dna(
        directory.path(),
        &config,
        &["notation", SAMPLE_ID, &document],
    )
    .success()
    .stdout(predicate::str::is_empty());
    let notation: Value = serde_json::from_slice(&fs::read(
        directory
            .path()
            .join("results")
            .join(format!("{SAMPLE_ID}.notation.json")),
    )?)?;

    assert_eq!(notation["schema_version"], "dna.notation/v1");
    assert_eq!(notation["notation"], variants["notation"]);
    assert_eq!(
        notation["notation"]["calls"].as_array().map(|calls| calls
            .iter()
            .map(|call| call["call"].clone())
            .collect::<Vec<_>>()),
        Some(vec![
            json!("73G"),
            json!("309.1C"),
            json!("309.2C"),
            json!("313.1G")
        ])
    );
    assert_eq!(
        notation["provenance"]["plugins"],
        json!([
            {"id": "normalization", "family": "post_calling", "version": 3},
            {"id": "nomenclature", "family": "post_calling", "version": 3},
            {"id": "conformance", "family": "post_calling", "version": 1},
            {"id": "empop", "family": "post_calling", "version": 2},
        ])
    );
    assert_eq!(
        notation["empop"],
        json!({
            "ranges": [[73, 340]],
            "ranges_text": "73-340",
            "mitotype": "73G 309.1C 309.2C 313.1G",
            "cut": [],
        })
    );
    assert_eq!(
        notation["provenance"]["source"],
        json!({
            "schema_version": "dna.variants/v1",
            "sha256": format!("{:x}", Sha256::digest(&variants_bytes)),
        })
    );
    assert_eq!(
        notation["conformance"],
        json!({
            "rules": ["insertion_at_run_end", "insertion_matches_run"],
            "findings": [{"rule": "insertion_matches_run", "read": "HV2", "calls": ["313.1G"]}],
        })
    );
    Ok(())
}

// SRS-NOM-018
#[test]
fn an_unresolved_base_splits_the_empop_range() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let config = directory.path().join("dna.toml");
    write_config_with_profile(&config, &human_mtdna_profile())?;
    let sequences = directory.path().join("consensus.fasta");
    // HVS-II with an `N` at 200 (before the 309 insertions, so the index
    // in the record is the reference position less 73).
    let mut record = hvs2(&rcrs()?);
    record.replace_range(200 - 73..200 - 72, "N");
    fs::write(&sequences, format!(">HV2\n{record}\n"))?;
    let sequences = sequences.to_string_lossy().into_owned();
    dna(directory.path(), &config, &["call", SAMPLE_ID, &sequences]).success();
    let document = directory
        .path()
        .join("results")
        .join(format!("{SAMPLE_ID}.variants.json"))
        .to_string_lossy()
        .into_owned();
    dna(
        directory.path(),
        &config,
        &["notation", SAMPLE_ID, &document],
    )
    .success();
    let notation: Value = serde_json::from_slice(&fs::read(
        directory
            .path()
            .join("results")
            .join(format!("{SAMPLE_ID}.notation.json")),
    )?)?;
    assert_eq!(notation["empop"]["ranges_text"], "73-199 201-340");
    assert_eq!(
        notation["empop"]["cut"],
        json!([{"start": 200, "end": 200, "reason": "unresolved"}])
    );
    assert_eq!(notation["empop"]["mitotype"], "73G 309.1C 309.2C 313.1G");
    Ok(())
}

// SRS-IN-006 SRS-IN-014
#[test]
fn rejects_a_document_of_another_sample() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let config = directory.path().join("dna.toml");
    write_config_with_profile(&config, &human_mtdna_profile())?;
    let sequences = directory.path().join("consensus.fasta");
    fs::write(&sequences, format!(">HV2\n{}\n", hvs2(&rcrs()?)))?;
    let sequences = sequences.to_string_lossy().into_owned();
    dna(directory.path(), &config, &["call", SAMPLE_ID, &sequences]).success();
    let document = directory
        .path()
        .join("results")
        .join(format!("{SAMPLE_ID}.variants.json"))
        .to_string_lossy()
        .into_owned();

    dna(directory.path(), &config, &["notation", "sample-2", &document])
        .failure()
        .stderr(predicate::str::contains(
            "invalid variants document: document sample \"sample-1\" is not the requested sample \"sample-2\"",
        ));
    assert!(
        !directory
            .path()
            .join("results/sample-2.notation.json")
            .exists()
    );
    Ok(())
}
