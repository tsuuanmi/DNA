//! CLI contract for reference-free read evidence (`dna.read_evidence/v1`,
//! PROP-0011): deterministic JSON that matches the basecalls of the same trace.

pub mod support;

use std::fs;
use std::path::{Path, PathBuf};

use assert_cmd::Command;
use predicates::prelude::*;
use serde_json::Value;
use tempfile::tempdir;

use support::{basecall_output_path, write_abif, write_config};

const QUERY: &str = "ACGTCAGTACGATCGTACCTGAGTACGA";

fn dna_binary() -> String {
    std::env::var("CARGO_BIN_EXE_dna")
        .unwrap_or_else(|error| panic!("Cargo must provide the DNA test binary: {error}"))
}

fn evidence_output_path(
    workdir: &Path,
    trace: &Path,
) -> Result<PathBuf, Box<dyn std::error::Error>> {
    let stem = trace
        .file_stem()
        .and_then(|stem| stem.to_str())
        .ok_or("trace stem is not UTF-8")?;
    Ok(workdir
        .join("results")
        .join(format!("{stem}.evidence.json")))
}

#[test]
fn writes_deterministic_evidence_that_matches_the_basecalls()
-> Result<(), Box<dyn std::error::Error>> {
    let first = tempdir()?;
    let second = tempdir()?;
    for directory in [first.path(), second.path()] {
        let trace = directory.join("trace.ab1");
        let config = directory.join("dna.toml");
        write_abif(&trace, QUERY)?;
        write_config(&config, "linear")?;
        run("evidence", &trace, &config, directory)?
            .success()
            .stdout(predicate::str::is_empty())
            .stderr(predicate::str::is_empty());
    }
    let first_trace = first.path().join("trace.ab1");
    let bytes = fs::read(evidence_output_path(first.path(), &first_trace)?)?;
    assert_eq!(
        bytes,
        fs::read(evidence_output_path(
            second.path(),
            &second.path().join("trace.ab1")
        )?)?
    );
    let value: Value = serde_json::from_slice(&bytes)?;
    assert_eq!(value["schema_version"], "dna.read_evidence/v1");
    assert_eq!(
        value["provenance"]["plugins"],
        serde_json::json!([{"id": "sanger", "family": "modality", "version": 12}])
    );
    let calls = value["calls"].as_array().ok_or("calls must be an array")?;
    assert_eq!(
        value["read"]["call_count"].as_u64(),
        u64::try_from(calls.len()).ok()
    );
    for (index, call) in calls.iter().enumerate() {
        assert_eq!(call["index"].as_u64(), u64::try_from(index).ok());
    }

    // The evidence's bases and callability are the basecalls' of the same trace.
    let config = first.path().join("dna.toml");
    run("basecall", &first_trace, &config, first.path())?.success();
    let basecalls: Value =
        serde_json::from_slice(&fs::read(basecall_output_path(first.path(), &first_trace))?)?;
    let bases = calls
        .iter()
        .map(|call| call["base"].as_str().unwrap_or_default())
        .collect::<String>();
    assert_eq!(Some(bases.as_str()), basecalls["read"]["primary"].as_str());
    assert_eq!(
        value["read"]["callability"],
        basecalls["read"]["callability"]
    );
    assert_eq!(
        value["provenance"]["configuration_sha256"],
        basecalls["provenance"]["configuration_sha256"]
    );
    Ok(())
}

#[test]
fn refuses_to_overwrite_an_existing_result() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let config = directory.path().join("dna.toml");
    write_abif(&trace, QUERY)?;
    write_config(&config, "linear")?;
    run("evidence", &trace, &config, directory.path())?.success();
    let output = evidence_output_path(directory.path(), &trace)?;
    let before = fs::read(&output)?;
    run("evidence", &trace, &config, directory.path())?.failure();
    assert_eq!(fs::read(&output)?, before);
    Ok(())
}

fn run(
    command: &str,
    trace: &Path,
    config: &Path,
    workdir: &Path,
) -> Result<assert_cmd::assert::Assert, Box<dyn std::error::Error>> {
    Ok(Command::new(dna_binary())
        .current_dir(workdir)
        .env("DNA_CONFIG", config)
        .env("DNA_LOG_DIR", workdir.join("logs"))
        .args([command, trace.to_str().ok_or("trace path is not UTF-8")?])
        .assert())
}
