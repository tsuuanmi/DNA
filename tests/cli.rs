use assert_cmd::Command;
use predicates::prelude::*;
use tempfile::tempdir;

fn dna_binary() -> String {
    std::env::var("CARGO_BIN_EXE_dna")
        .unwrap_or_else(|error| panic!("Cargo must provide the DNA test binary: {error}"))
}

#[test]
fn help_succeeds() {
    let mut command = Command::new(dna_binary());
    command
        .arg("--help")
        .assert()
        .success()
        .stdout(predicate::str::contains("Process Sanger sequencing traces"));
}

#[test]
fn basecall_requires_a_trace() {
    let mut command = Command::new(dna_binary());
    command
        .arg("basecall")
        .assert()
        .failure()
        .stderr(predicate::str::contains("Usage:"));
}

#[test]
fn basecall_rejects_reference_option() {
    let mut command = Command::new(dna_binary());
    command
        .args(["basecall", "trace.ab1", "--reference", "reference.fa"])
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "unexpected argument '--reference'",
        ));
}

#[test]
fn analyze_requires_arguments() {
    let mut command = Command::new(dna_binary());
    command
        .arg("analyze")
        .assert()
        .failure()
        .stderr(predicate::str::contains("Usage:"));
}

#[test]
fn rejects_removed_out_prefix_option() {
    let mut command = Command::new(dna_binary());
    command
        .args([
            "analyze",
            "trace.ab1",
            "--reference",
            "reference.fa",
            "--out-prefix",
            "legacy",
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "unexpected argument '--out-prefix'",
        ));
}

#[test]
fn missing_trace_fails_explicitly() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let mut command = Command::new(dna_binary());
    let log_directory = directory.path().join("custom-logs");
    command
        .current_dir(directory.path())
        .env("DNA_LOG_DIR", &log_directory)
        .args([
            "analyze",
            "missing.ab1",
            "--reference",
            "references/rCRS.fasta",
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains("failed to read AB1 file"));
    let log = std::fs::read_to_string(log_directory.join("missing.log"))?;
    assert!(log.contains("event=analysis_started"));
    assert!(log.contains("event=analysis_failed stage=input_loading"));
    assert!(log.contains("failed to read AB1 file"));
    Ok(())
}

#[test]
fn sample_requires_at_least_one_trace() {
    let mut command = Command::new(dna_binary());
    command
        .args(["sample", "sample-1", "--reference", "reference.fa"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("Usage:"));
}

#[test]
fn sample_rejects_unsafe_sample_id_before_io() {
    let mut command = Command::new(dna_binary());
    command
        .args([
            "sample",
            "../sample",
            "missing.ab1",
            "--reference",
            "missing.fa",
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "sample id must be 1..=128 ASCII characters",
        ));
}
