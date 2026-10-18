//! CLI contract for a run signal: deterministic JSON over a run's traces.

pub mod support;

use std::fs;
use std::path::Path;

use assert_cmd::Command;
use predicates::prelude::*;
use serde_json::Value;
use tempfile::tempdir;

use support::{write_abif, write_abif_with_secondary_signal, write_config};

const QUERY: &str = "ACGTCAGTACGATCGTACCTGAGTACGATTGCAGCTAGCATGCATCGATGCA";

fn dna_binary() -> String {
    std::env::var("CARGO_BIN_EXE_dna")
        .unwrap_or_else(|error| panic!("Cargo must provide the DNA test binary: {error}"))
}

fn run(
    workdir: &Path,
    config: &Path,
    traces: &[&Path],
) -> Result<assert_cmd::assert::Assert, Box<dyn std::error::Error>> {
    let mut command = Command::new(dna_binary());
    command
        .current_dir(workdir)
        .env("DNA_CONFIG", config)
        .env("DNA_LOG_DIR", workdir.join("logs"))
        .args(["signal", "RUN1"]);
    for trace in traces {
        command.arg(trace.to_str().ok_or("trace path is not UTF-8")?);
    }
    Ok(command.assert())
}

#[test]
fn writes_a_deterministic_run_signal_whatever_the_trace_order()
-> Result<(), Box<dyn std::error::Error>> {
    let mut outputs = Vec::new();
    for reversed in [false, true] {
        let directory = tempdir()?;
        let config = directory.path().join("dna.toml");
        write_config(&config, "linear")?;
        let first = directory.path().join("a.ab1");
        let second = directory.path().join("b.ab1");
        let short = directory.path().join("c.ab1");
        write_abif(&first, QUERY)?;
        write_abif_with_secondary_signal(&second, QUERY, 20, b'T', 300)?;
        write_abif(&short, &QUERY[..9])?;
        let mut traces = vec![first.as_path(), second.as_path(), short.as_path()];
        if reversed {
            traces.reverse();
        }
        run(directory.path(), &config, &traces)?
            .success()
            .stdout(predicate::str::is_empty());
        outputs.push(fs::read(directory.path().join("results/RUN1.signal.json"))?);
    }
    assert_eq!(outputs[0], outputs[1]);
    let value: Value = serde_json::from_slice(&outputs[0])?;
    assert_eq!(value["schema_version"], "dna.signal/v1");
    assert_eq!(value["run_id"], "RUN1");
    assert_eq!(value["traces"].as_array().map(Vec::len), Some(3));
    assert_eq!(value["reads"], 2);
    assert_eq!(value["rejected_reads"], 1);
    let total = value["keys"]
        .as_array()
        .ok_or("keys must be an array")?
        .iter()
        .flat_map(|key| key["shares"].as_array().cloned().unwrap_or_default())
        .map(|pair| pair[1].as_u64().unwrap_or(0))
        .sum::<u64>();
    // Each selected call counts under its three other channels.
    assert_eq!(value["calls"].as_u64().map(|calls| 3 * calls), Some(total));
    Ok(())
}

#[test]
fn refuses_an_invalid_run_identifier() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let config = directory.path().join("dna.toml");
    write_config(&config, "linear")?;
    let trace = directory.path().join("a.ab1");
    write_abif(&trace, QUERY)?;
    Command::new(dna_binary())
        .current_dir(directory.path())
        .env("DNA_CONFIG", &config)
        .env("DNA_LOG_DIR", directory.path().join("logs"))
        .args([
            "signal",
            "../run",
            trace.to_str().ok_or("trace path is not UTF-8")?,
        ])
        .assert()
        .failure();
    assert!(!directory.path().join("results").exists());
    Ok(())
}
