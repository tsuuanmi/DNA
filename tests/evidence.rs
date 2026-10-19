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
    analysis_output_path, basecall_output_path, write_abif, write_abif_with_incoherent_doubles,
    write_abif_with_secondary_signal, write_config, write_reference,
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
        serde_json::json!([{"id": "sanger", "family": "modality", "version": 13}])
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

/// Runs `dna` with `arguments` in `workdir` under `config`.
fn dna(
    arguments: &[&std::ffi::OsStr],
    workdir: &Path,
    config: &Path,
) -> assert_cmd::assert::Assert {
    Command::new(dna_binary())
        .current_dir(workdir)
        .env("DNA_CONFIG", config)
        .env("DNA_LOG_DIR", workdir.join("logs"))
        .args(arguments)
        .assert()
}

/// A sample's inputs: a reference and a configuration whose profile declares
/// notation, and five reads: a forward read with a substitution and a second
/// peak, the same read on the reverse strand, a read with a substitution
/// only, a read too short for one signal window, and a read with too few
/// callable calls.
fn sample_inputs(
    directory: &Path,
) -> Result<(PathBuf, PathBuf, Vec<PathBuf>), Box<dyn std::error::Error>> {
    let reference = directory.join("reference.fa");
    let config = directory.join("dna.toml");
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;
    let profile = config.with_extension("profile.toml");
    let mut text = fs::read_to_string(&profile)?;
    text.push_str(
        "[normalization]\nindel_placement='right'\nmax_trade=20\n[notation]\nstyle='per_base_decimal'\n",
    );
    fs::write(&profile, text)?;
    let substituted = format!("{}A{}", &QUERY[..10], &QUERY[11..]);
    let traces = ["forward", "reverse", "plain", "short", "doubled"]
        .map(|name| directory.join(format!("read-{name}.ab1")));
    write_abif_with_secondary_signal(&traces[0], &substituted, 14, b'T', 400)?;
    write_abif(&traces[1], &reverse_complement(&substituted))?;
    write_abif(&traces[2], &substituted)?;
    write_abif(&traces[3], &QUERY[..5])?;
    write_abif_with_incoherent_doubles(&traces[4], QUERY, 0..QUERY.len(), 500)?;
    Ok((reference, config, traces.to_vec()))
}

// SRS-API-019 SRS-OUT-007
#[test]
fn sample_commands_from_evidence_write_what_the_traces_write()
-> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let directory = directory.path();
    let (reference, config, traces) = sample_inputs(directory)?;
    let workdir = |name: &str| -> Result<PathBuf, Box<dyn std::error::Error>> {
        let path = directory.join(name);
        fs::create_dir(&path)?;
        Ok(path)
    };

    // Every read's evidence, a rejected read's included.
    let recorded = workdir("recorded")?;
    let mut documents = Vec::new();
    for trace in &traces {
        run("evidence", trace, &config, &recorded)?.success();
        documents.push(evidence_output_path(&recorded, trace)?);
    }
    let short: Value = serde_json::from_slice(&fs::read(&documents[3])?)?;
    assert_eq!(
        short["rejection"],
        serde_json::json!({"reason": "too_short", "calls": 5, "minimum_calls": 10})
    );
    let doubled: Value = serde_json::from_slice(&fs::read(&documents[4])?)?;
    assert_eq!(
        doubled["rejection"]["reason"],
        "callable_calls_below_minimum"
    );
    assert!(doubled["callability"].is_object() && doubled.get("calls").is_none());

    // The run signal of the sample's run, which the review takes out again.
    let run_id = "run-1";
    let mut arguments = vec!["signal".as_ref(), run_id.as_ref()];
    arguments.extend(traces.iter().map(|trace| trace.as_os_str()));
    dna(&arguments, &recorded, &config).success();
    let signal = recorded
        .join("results")
        .join(format!("{run_id}.signal.json"));

    // Every read, then the rejected reads alone: a review without an
    // admitted read.
    for (command, outputs, with_signal, reads) in [
        ("sample", &["sample.json"][..], false, 0..5),
        (
            "consensus",
            &["consensus.json", "consensus.fasta"][..],
            false,
            0..5,
        ),
        ("review", &["review.json"][..], false, 0..5),
        ("review", &["review.json"][..], true, 0..5),
        ("review", &["review.json"][..], false, 3..5),
    ] {
        let name = format!("{command}-{with_signal}-{}", reads.start);
        let from_traces = workdir(&format!("{name}-traces"))?;
        let from_evidence = workdir(&format!("{name}-evidence"))?;
        let mut tail = vec!["--reference".as_ref(), reference.as_os_str()];
        if with_signal {
            tail.extend(["--signal".as_ref(), signal.as_os_str()]);
        }
        let mut arguments = vec![command.as_ref(), "sample-1".as_ref()];
        arguments.extend(traces[reads.clone()].iter().map(|trace| trace.as_os_str()));
        arguments.extend(&tail);
        dna(&arguments, &from_traces, &config).success();
        let mut arguments = vec![command.as_ref(), "sample-1".as_ref(), "--evidence".as_ref()];
        arguments.extend(documents[reads].iter().map(|document| document.as_os_str()));
        arguments.extend(&tail);
        dna(&arguments, &from_evidence, &config)
            .success()
            .stdout(predicate::str::is_empty())
            .stderr(predicate::str::is_empty());
        for output in outputs {
            let file = format!("sample-1.{output}");
            let expected = fs::read(from_traces.join("results").join(&file))?;
            let actual = fs::read(from_evidence.join("results").join(&file))?;
            assert_eq!(actual, expected, "{name}: {file}");
        }
    }

    // The documents name the reads after their traces and list the rejected
    // reads with their causes.
    let sample: Value = serde_json::from_slice(&fs::read(
        directory
            .join("sample-false-0-evidence")
            .join("results")
            .join("sample-1.sample.json"),
    )?)?;
    let names = sample["reads"]
        .as_array()
        .ok_or("reads must be an array")?
        .iter()
        .filter_map(|read| read["name"].as_str())
        .collect::<std::collections::BTreeSet<_>>();
    assert_eq!(names, ["read-forward", "read-plain", "read-reverse"].into());
    let rejected = sample["rejected_reads"]
        .as_array()
        .ok_or("rejected reads must be an array")?
        .iter()
        .map(|read| (read["name"].as_str(), read["rejection"]["reason"].as_str()))
        .collect::<std::collections::BTreeSet<_>>();
    assert_eq!(
        rejected,
        [
            (Some("read-doubled"), Some("callable_calls_below_minimum")),
            (Some("read-short"), Some("too_short")),
        ]
        .into()
    );
    let review: Value = serde_json::from_slice(&fs::read(
        directory
            .join("review-true-0-evidence")
            .join("results")
            .join("sample-1.review.json"),
    )?)?;
    assert_eq!(review["provenance"]["signal"]["run_id"], run_id);
    let rejected: Value = serde_json::from_slice(&fs::read(
        directory
            .join("review-false-3-evidence")
            .join("results")
            .join("sample-1.review.json"),
    )?)?;
    assert_eq!(rejected["status"], "rejected");
    Ok(())
}

// SRS-API-017 SRS-API-018
#[test]
fn a_rejected_read_fails_its_analysis_as_its_trace_does() -> Result<(), Box<dyn std::error::Error>>
{
    let directory = tempdir()?;
    let directory = directory.path();
    let (reference, config, traces) = sample_inputs(directory)?;
    for trace in &traces[3..] {
        let stem = trace
            .file_stem()
            .and_then(|stem| stem.to_str())
            .ok_or("trace stem is not UTF-8")?;
        let from_trace = directory.join(format!("{stem}-trace"));
        let from_evidence = directory.join(format!("{stem}-evidence"));
        fs::create_dir(&from_trace)?;
        fs::create_dir(&from_evidence)?;
        let expected = dna(
            &[
                "analyze".as_ref(),
                trace.as_os_str(),
                "--reference".as_ref(),
                reference.as_os_str(),
            ],
            &from_trace,
            &config,
        )
        .failure()
        .get_output()
        .stderr
        .clone();
        run("evidence", trace, &config, &from_evidence)?.success();
        let document = evidence_output_path(&from_evidence, trace)?;
        let actual = dna(
            &[
                "analyze".as_ref(),
                "--evidence".as_ref(),
                document.as_os_str(),
                "--reference".as_ref(),
                reference.as_os_str(),
            ],
            &from_evidence,
            &config,
        )
        .failure()
        .get_output()
        .stderr
        .clone();
        assert_eq!(actual, expected, "{stem}");
        assert!(
            !from_evidence
                .join("results")
                .join(format!("{stem}.json"))
                .exists()
        );
    }
    Ok(())
}

// SRS-API-019 SRS-IN-001
#[test]
fn sample_commands_take_traces_or_evidence_not_both() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let directory = directory.path();
    let (reference, config, traces) = sample_inputs(directory)?;
    run("evidence", &traces[0], &config, directory)?.success();
    let document = evidence_output_path(directory, &traces[0])?;
    for command in ["sample", "consensus", "review"] {
        dna(
            &[
                command.as_ref(),
                "sample-1".as_ref(),
                traces[1].as_os_str(),
                "--evidence".as_ref(),
                document.as_os_str(),
                "--reference".as_ref(),
                reference.as_os_str(),
            ],
            directory,
            &config,
        )
        .failure()
        .stderr(predicate::str::contains("cannot be used with"));
        dna(
            &[
                command.as_ref(),
                "sample-1".as_ref(),
                "--reference".as_ref(),
                reference.as_os_str(),
            ],
            directory,
            &config,
        )
        .failure();
    }
    // A trace given as evidence is no read evidence document.
    dna(
        &[
            "sample".as_ref(),
            "sample-1".as_ref(),
            "--evidence".as_ref(),
            traces[1].as_os_str(),
            "--reference".as_ref(),
            reference.as_os_str(),
        ],
        directory,
        &config,
    )
    .failure()
    .stderr(predicate::str::contains("invalid read evidence"));
    assert!(
        !directory
            .join("results")
            .join("sample-1.sample.json")
            .exists()
    );
    Ok(())
}
