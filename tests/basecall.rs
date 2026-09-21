mod support;

use std::collections::BTreeSet;
use std::fs;
use std::path::Path;

use assert_cmd::Command;
use predicates::prelude::*;
use serde_json::Value;
use tempfile::tempdir;

use support::{
    analysis_output_path, basecall_output_path, write_abif, write_config, write_reference,
};

const QUERY: &str = "ACGTCAGTACGATCGTACCTGAGTACGA";

fn dna_binary() -> String {
    std::env::var("CARGO_BIN_EXE_dna")
        .unwrap_or_else(|error| panic!("Cargo must provide the DNA test binary: {error}"))
}

#[test]
fn writes_deterministic_reference_free_json() -> Result<(), Box<dyn std::error::Error>> {
    let first = tempdir()?;
    let second = tempdir()?;
    for directory in [first.path(), second.path()] {
        let trace = directory.join("trace.ab1");
        let config = directory.join("dna.toml");
        write_abif(&trace, QUERY)?;
        write_config(&config, "linear")?;
        run(&trace, &config, directory)?
            .success()
            .stdout(predicate::str::is_empty())
            .stderr(predicate::str::is_empty());
        assert!(!directory.join("reference.fa").exists());
    }

    let first_bytes = fs::read(basecall_output_path(
        first.path(),
        &first.path().join("trace.ab1"),
    ))?;
    let second_bytes = fs::read(basecall_output_path(
        second.path(),
        &second.path().join("trace.ab1"),
    ))?;
    assert_eq!(first_bytes, second_bytes);
    let value: Value = serde_json::from_slice(&first_bytes)?;
    assert_eq!(value["schema_version"], "dna.basecalls/v2");
    assert_object_keys(
        &value,
        &[
            "schema_version",
            "provenance",
            "read",
            "signal_quality",
            "warnings",
        ],
    );
    assert_object_keys(&value["provenance"], &["input", "configuration_sha256"]);
    assert_object_keys(
        &value["read"],
        &["call_count", "primary", "ambiguity", "retained", "trim"],
    );
    assert_object_keys(
        &value["warnings"],
        &[
            "unresolved_primary_calls",
            "multi_channel_unresolved_calls",
            "vendor_disagreements",
            "ploc_vendor_length_mismatches",
            "clipped_channel_samples",
        ],
    );
    assert_eq!(value["read"]["call_count"], QUERY.len());
    assert_eq!(value["read"]["primary"], QUERY);
    assert_eq!(value["read"]["ambiguity"], QUERY);
    assert_eq!(value["read"]["retained"], QUERY);
    assert_eq!(value["read"]["trim"]["start"], 0);
    assert_eq!(value["read"]["trim"]["end"], QUERY.len());
    assert!(value["signal_quality"]["noisy_regions"].is_array());
    assert_eq!(
        value["signal_quality"]["integrity"]["ploc_count"],
        QUERY.len()
    );
    assert_eq!(
        value["signal_quality"]["integrity"]["vendor_primary_count"],
        QUERY.len()
    );
    assert_eq!(
        value["signal_quality"]["integrity"]["vendor_quality_count"],
        QUERY.len()
    );
    assert_eq!(
        value["signal_quality"]["integrity"]["clipped_channel_samples"],
        0
    );
    assert!(value.get("reference").is_none());
    assert!(value.get("alignment").is_none());
    assert!(value.get("variants").is_none());

    let log = fs::read_to_string(first.path().join("logs/trace.log"))?;
    let mut search_start = 0;
    for event in [
        "event=basecall_started",
        "event=basecall_inputs_loaded",
        "event=basecalling_completed",
        "event=signal_processing_completed",
        "event=quality_control_completed",
        "event=basecall_ready_for_publication",
    ] {
        let offset = log[search_start..]
            .find(event)
            .ok_or_else(|| format!("missing ordered log event {event}"))?;
        search_start += offset + event.len();
    }
    assert!(!log.contains(QUERY));
    assert!(!log.contains("\"schema_version\""));
    Ok(())
}

#[test]
fn malformed_abif_leaves_no_basecall_output() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let config = directory.path().join("dna.toml");
    fs::write(&trace, b"not an ABIF file")?;
    write_config(&config, "linear")?;

    run(&trace, &config, directory.path())?
        .failure()
        .stderr(predicate::str::contains("invalid ABIF input"));
    assert!(!basecall_output_path(directory.path(), &trace).exists());
    let log = fs::read_to_string(directory.path().join("logs/trace.log"))?;
    assert!(log.contains("event=basecall_failed stage=input_loading"));
    Ok(())
}

#[test]
fn refuses_to_overwrite_basecall_output() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let config = directory.path().join("dna.toml");
    write_abif(&trace, QUERY)?;
    write_config(&config, "linear")?;
    let output = basecall_output_path(directory.path(), &trace);
    fs::create_dir_all(output.parent().ok_or("output has no parent")?)?;
    fs::write(&output, b"owned")?;

    run(&trace, &config, directory.path())?
        .failure()
        .stderr(predicate::str::contains("target already exists"));
    assert_eq!(fs::read(output)?, b"owned");
    Ok(())
}

#[test]
fn basecall_and_analysis_outputs_coexist() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    write_abif(&trace, QUERY)?;
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;

    run(&trace, &config, directory.path())?.success();
    analyze(&trace, &reference, &config, directory.path())?.success();
    assert!(basecall_output_path(directory.path(), &trace).exists());
    assert!(analysis_output_path(directory.path(), &trace).exists());
    Ok(())
}

fn run(
    trace: &Path,
    config: &Path,
    workdir: &Path,
) -> Result<assert_cmd::assert::Assert, Box<dyn std::error::Error>> {
    let mut command = Command::new(dna_binary());
    Ok(command
        .current_dir(workdir)
        .env("DNA_CONFIG", config)
        .env("DNA_LOG_DIR", workdir.join("logs"))
        .args(["basecall", trace.to_str().ok_or("trace path is not UTF-8")?])
        .assert())
}

fn analyze(
    trace: &Path,
    reference: &Path,
    config: &Path,
    workdir: &Path,
) -> Result<assert_cmd::assert::Assert, Box<dyn std::error::Error>> {
    let mut command = Command::new(dna_binary());
    Ok(command
        .current_dir(workdir)
        .env("DNA_CONFIG", config)
        .env("DNA_LOG_DIR", workdir.join("logs"))
        .args([
            "analyze",
            trace.to_str().ok_or("trace path is not UTF-8")?,
            "--reference",
            reference.to_str().ok_or("reference path is not UTF-8")?,
        ])
        .assert())
}

fn assert_object_keys(value: &Value, expected: &[&str]) {
    let actual = value
        .as_object()
        .map(|object| object.keys().map(String::as_str).collect::<BTreeSet<_>>())
        .unwrap_or_default();
    let expected = expected.iter().copied().collect::<BTreeSet<_>>();
    assert_eq!(actual, expected);
}
