//! CLI contract for reference-free read evidence (`dna.read_evidence/v1`,
//! PROP-0011): deterministic JSON that matches the basecalls of the same
//! trace, and from which `analyze --evidence` reproduces the trace's analysis
//! byte for byte.

pub mod support;

use std::fs;
use std::path::{Path, PathBuf};

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

// SRS-API-017 SRS-OUT-001 SRS-OUT-007
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
        serde_json::json!([{"id": "sanger", "family": "modality", "version": 14}])
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

// SRS-API-017 SRS-OUT-006
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

/// The reverse complement of `sequence`.
fn reverse_complement(sequence: &str) -> String {
    sequence
        .chars()
        .rev()
        .map(|base| match base {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            _ => 'A',
        })
        .collect()
}

// SRS-API-018 SRS-IN-010
#[test]
fn analyzes_a_read_from_its_evidence_exactly_as_from_its_trace()
-> Result<(), Box<dyn std::error::Error>> {
    // A read that matches its reference, one with a substitution, and the
    // same read on the reverse strand.
    let substituted = format!("{}T{}", &QUERY[..10], &QUERY[11..]);
    let cases = [
        ("matching", QUERY.to_owned(), QUERY.to_owned()),
        ("forward", QUERY.to_owned(), substituted.clone()),
        ("reverse", reverse_complement(QUERY), substituted),
    ];
    for (name, read, reference_sequence) in cases {
        let from_trace = tempdir()?;
        let from_evidence = tempdir()?;
        for directory in [from_trace.path(), from_evidence.path()] {
            write_abif(&directory.join(format!("{name}.ab1")), &read)?;
            write_reference(
                &directory.join("reference.fa"),
                &format!("TTTT{reference_sequence}CCCC"),
            )?;
            write_config(&directory.join("dna.toml"), "linear")?;
        }
        let trace = from_trace.path().join(format!("{name}.ab1"));
        analyze(
            &[trace.as_os_str()],
            from_trace.path(),
            &from_trace.path().join("dna.toml"),
        )
        .success();

        let directory = from_evidence.path();
        let config = directory.join("dna.toml");
        let evidence_trace = directory.join(format!("{name}.ab1"));
        run("evidence", &evidence_trace, &config, directory)?.success();
        let evidence = evidence_output_path(directory, &evidence_trace)?;
        analyze(
            &["--evidence".as_ref(), evidence.as_os_str()],
            directory,
            &config,
        )
        .success()
        .stdout(predicate::str::is_empty())
        .stderr(predicate::str::is_empty());

        // Published under the trace's stem, byte-identical.
        let expected = fs::read(analysis_output_path(from_trace.path(), &trace))?;
        let actual = fs::read(analysis_output_path(directory, &evidence_trace))?;
        assert_eq!(actual, expected, "{name}");
        let value: Value = serde_json::from_slice(&actual)?;
        let variants = value["variants"].as_array().map_or(0, Vec::len);
        assert_eq!(variants > 0, name != "matching", "{name}");
    }
    Ok(())
}

// SRS-API-018
#[test]
fn refuses_evidence_it_cannot_read_back() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let directory = directory.path();
    let trace = directory.join("trace.ab1");
    let config = directory.join("dna.toml");
    write_abif(&trace, QUERY)?;
    write_reference(&directory.join("reference.fa"), &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;
    run("evidence", &trace, &config, directory)?.success();
    let document: Value =
        serde_json::from_slice(&fs::read(evidence_output_path(directory, &trace)?)?)?;

    let mut schema = document.clone();
    schema["schema_version"] = "dna.read_evidence/v0".into();
    let mut reason = document.clone();
    reason["read"]["reasons"][0] = "made_up_reason".into();
    let mut order = document.clone();
    order["calls"][0]["index"] = 1.into();
    order["calls"][1]["index"] = 0.into();
    let mut configuration = document;
    configuration["provenance"]["configuration_sha256"] = "0".repeat(64).into();
    for (name, edited, message) in [
        ("schema", schema, "unsupported schema version"),
        ("reason", reason, "unknown reason label \"made_up_reason\""),
        ("order", order, "call at position 0 has index 1"),
        (
            "configuration",
            configuration,
            "is not the loaded configuration",
        ),
    ] {
        let path = directory.join(format!("{name}.evidence.json"));
        fs::write(&path, serde_json::to_vec(&edited)?)?;
        analyze(
            &["--evidence".as_ref(), path.as_os_str()],
            directory,
            &config,
        )
        .failure()
        .stderr(
            predicate::str::contains("invalid read evidence")
                .and(predicate::str::contains(message)),
        );
        assert!(
            !directory
                .join("results")
                .join(format!("{name}.json"))
                .exists()
        );
    }

    // A trace and its evidence together, or neither, is no analysis.
    let evidence = evidence_output_path(directory, &trace)?;
    analyze(
        &[
            trace.as_os_str(),
            "--evidence".as_ref(),
            evidence.as_os_str(),
        ],
        directory,
        &config,
    )
    .failure()
    .stderr(predicate::str::contains("cannot be used with"));
    analyze(&[], directory, &config).failure();
    assert!(!directory.join("results").join("trace.json").exists());
    Ok(())
}

/// Runs `analyze` with `arguments` and the reference `reference.fa` in `workdir`.
fn analyze(
    arguments: &[&std::ffi::OsStr],
    workdir: &Path,
    config: &Path,
) -> assert_cmd::assert::Assert {
    Command::new(dna_binary())
        .current_dir(workdir)
        .env("DNA_CONFIG", config)
        .env("DNA_LOG_DIR", workdir.join("logs"))
        .arg("analyze")
        .args(arguments)
        .arg("--reference")
        .arg(workdir.join("reference.fa"))
        .assert()
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
