//! CLI contract for the sample consensus (PROP-0003).

pub mod support;

use std::fs;
use std::path::{Path, PathBuf};

use assert_cmd::Command;
use serde_json::Value;
use tempfile::tempdir;

use support::{write_abif, write_config, write_reference};

const QUERY: &str = "ACGTCAGTACGATCGTACCTGAGTACGA";
const SAMPLE_ID: &str = "sample-1";

fn dna_binary() -> String {
    std::env::var("CARGO_BIN_EXE_dna")
        .unwrap_or_else(|error| panic!("Cargo must provide the DNA test binary: {error}"))
}

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

/// `QUERY` with its 11th base (1-based position 15 on the reference) as `A`.
fn alternate() -> String {
    let mut bases = QUERY.as_bytes().to_vec();
    bases[10] = b'A';
    String::from_utf8_lossy(&bases).into_owned()
}

fn run(
    traces: &[&PathBuf],
    reference: &Path,
    config: &Path,
    workdir: &Path,
) -> assert_cmd::assert::Assert {
    let mut command = Command::new(dna_binary());
    command
        .current_dir(workdir)
        .env("DNA_CONFIG", config)
        .arg("consensus")
        .arg(SAMPLE_ID);
    for trace in traces {
        command.arg(trace);
    }
    command.arg("--reference").arg(reference).assert()
}

fn document(workdir: &Path) -> Result<Value, Box<dyn std::error::Error>> {
    let path = workdir
        .join("results")
        .join(format!("{SAMPLE_ID}.consensus.json"));
    Ok(serde_json::from_slice(&fs::read(path)?)?)
}

/// Two reads that agree on a substitution give it in the consensus sequence,
/// its FASTA, and one called site naming both reads.
#[test]
fn agreeing_reads_call_their_difference() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let forward = directory.path().join("read-forward.ab1");
    let reverse = directory.path().join("read-reverse.ab1");
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;
    write_abif(&forward, &alternate())?;
    write_abif(&reverse, &reverse_complement(&alternate()))?;

    run(&[&forward, &reverse], &reference, &config, directory.path()).success();

    let value = document(directory.path())?;
    assert_eq!(value["schema_version"], "dna.consensus/v1");
    assert_eq!(value["provenance"]["plugins"][2]["id"], "consensus");
    let segments = value["segments"]
        .as_array()
        .ok_or("segments must be an array")?;
    assert_eq!(segments.len(), 1);
    assert_eq!(segments[0]["sequence"], alternate());
    assert_eq!(segments[0]["reference"]["start"], 4);
    // The substituted A joins the A after it: the site is that whole run.
    let sites = value["sites"].as_array().ok_or("sites must be an array")?;
    assert_eq!(sites.len(), 1);
    assert_eq!(
        (&sites[0]["start"], &sites[0]["end"]),
        (&15.into(), &16.into())
    );
    assert_eq!(sites[0]["reference"], "GA");
    assert_eq!(sites[0]["call"], "AA");
    assert_eq!(sites[0]["state"], "called");
    assert_eq!(
        sites[0]["runs"],
        serde_json::json!([{"base": "A", "length": 2, "length_evidence": "in_phase"}])
    );
    assert_eq!(value["summary"]["reference_frame_runs"], 0);
    assert_eq!(
        sites[0]["supporting_reads"],
        serde_json::json!(["read-forward", "read-reverse"])
    );
    let fasta = fs::read_to_string(
        directory
            .path()
            .join("results")
            .join(format!("{SAMPLE_ID}.consensus.fasta")),
    )?;
    assert_eq!(fasta, format!(">{SAMPLE_ID}_5-32\n{}\n", alternate()));
    Ok(())
}

/// A read with fewer calls than one signal window is listed as rejected and
/// the consensus is built from the other reads (SRS-CONS-001, SRS-SAMPLE-030).
#[test]
fn a_read_shorter_than_the_signal_window_is_rejected() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let forward = directory.path().join("read-forward.ab1");
    let reverse = directory.path().join("read-reverse.ab1");
    let short = directory.path().join("read-short.ab1");
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;
    write_abif(&forward, &alternate())?;
    write_abif(&reverse, &reverse_complement(&alternate()))?;
    write_abif(&short, &QUERY[..5])?;

    run(
        &[&forward, &short, &reverse],
        &reference,
        &config,
        directory.path(),
    )
    .success();

    let value = document(directory.path())?;
    assert_eq!(value["reads"].as_array().map(Vec::len), Some(2));
    let rejected = value["rejected_reads"]
        .as_array()
        .ok_or("rejected_reads must be an array")?;
    assert_eq!(rejected.len(), 1);
    assert_eq!(rejected[0]["name"], "read-short");
    assert_eq!(
        rejected[0]["rejection"],
        serde_json::json!({"reason": "too_short", "calls": 5, "minimum_calls": 10})
    );
    assert_eq!(value["sites"][0]["state"], "called");
    Ok(())
}

/// One read against one read keeps the reference base and marks the site
/// contested.
#[test]
fn a_lone_difference_against_a_reference_read_is_contested()
-> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let forward = directory.path().join("read-forward.ab1");
    let reverse = directory.path().join("read-reverse.ab1");
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;
    write_abif(&forward, &alternate())?;
    write_abif(&reverse, &reverse_complement(QUERY))?;

    run(&[&forward, &reverse], &reference, &config, directory.path()).success();

    let value = document(directory.path())?;
    assert_eq!(value["segments"][0]["sequence"], QUERY);
    assert_eq!(value["sites"][0]["state"], "contested");
    assert_eq!(value["sites"][0]["call"], "GA");
    assert_eq!(
        value["sites"][0]["opposing_reads"],
        serde_json::json!(["read-forward"])
    );
    assert_eq!(value["summary"]["contested_sites"], 1);
    Ok(())
}

/// A log line without what differs between runs: the timestamp, the run id,
/// and timings.
#[cfg(target_os = "linux")]
fn stable_record(line: &str) -> String {
    let mut parts = line.splitn(3, " | ");
    let (_time, level, rest) = (parts.next(), parts.next(), parts.next());
    let fields = rest
        .unwrap_or_default()
        .split(' ')
        .filter(|field| !field.starts_with("run_id="))
        .map(|field| match field.split_once('=') {
            Some((name, _)) if name.ends_with("elapsed_ms") => format!("{name}=_"),
            _ => field.to_owned(),
        })
        .collect::<Vec<_>>();
    format!("{} | {}", level.unwrap_or_default(), fields.join(" "))
}

/// Reads on one CPU and on every available CPU give the same documents and the
/// same log, apart from timestamps, run ids, and timings (SRS-SAMPLE-001,
/// SRS-OUT-008).
#[cfg(target_os = "linux")]
#[test]
fn one_worker_and_many_publish_the_same_documents_and_log() -> Result<(), Box<dyn std::error::Error>>
{
    let directory = tempdir()?;
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;
    let mut traces = Vec::new();
    for (name, sequence) in [
        ("read-forward", alternate()),
        ("read-foreign", "AGCATGAGTCCATGCTAGCATGACTGCA".to_owned()),
        ("read-reverse", reverse_complement(&alternate())),
        ("read-reference", QUERY.to_owned()),
        ("read-reference-reverse", reverse_complement(QUERY)),
    ] {
        let trace = directory.path().join(format!("{name}.ab1"));
        write_abif(&trace, &sequence)?;
        traces.push(trace);
    }
    let mut outputs = Vec::new();
    for pinned in [true, false] {
        let workdir = directory.path().join(if pinned { "one" } else { "many" });
        fs::create_dir(&workdir)?;
        let mut command = if pinned {
            let mut taskset = Command::new("taskset");
            taskset.arg("-c").arg("0").arg(dna_binary());
            taskset
        } else {
            Command::new(dna_binary())
        };
        command
            .current_dir(&workdir)
            .env("DNA_CONFIG", &config)
            .arg("consensus")
            .arg(SAMPLE_ID)
            .args(&traces)
            .arg("--reference")
            .arg(&reference)
            .assert()
            .success();
        let results = workdir.join("results");
        let log = fs::read_to_string(workdir.join("logs").join(format!("{SAMPLE_ID}.log")))?;
        outputs.push((
            fs::read(results.join(format!("{SAMPLE_ID}.consensus.json")))?,
            fs::read(results.join(format!("{SAMPLE_ID}.consensus.fasta")))?,
            log.lines().map(stable_record).collect::<Vec<_>>(),
        ));
    }
    assert_eq!(outputs[0], outputs[1]);
    assert!(
        outputs[0]
            .2
            .iter()
            .any(|line| line.contains("read_index=4"))
    );
    Ok(())
}

/// `call` takes the consensus document on its own, of the requested sample
/// (PROP-0006).
#[test]
fn call_takes_a_consensus_document_alone_and_of_its_sample()
-> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let forward = directory.path().join("read-forward.ab1");
    let reverse = directory.path().join("read-reverse.ab1");
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;
    write_abif(&forward, &alternate())?;
    write_abif(&reverse, &reverse_complement(&alternate()))?;
    run(&[&forward, &reverse], &reference, &config, directory.path()).success();
    let results = directory.path().join("results");
    let document = results.join(format!("{SAMPLE_ID}.consensus.json"));
    let fasta = results.join(format!("{SAMPLE_ID}.consensus.fasta"));
    let call = |sample: &str, inputs: &[&PathBuf]| {
        let mut command = Command::new(dna_binary());
        command
            .current_dir(directory.path())
            .env("DNA_CONFIG", &config)
            .arg("call")
            .arg(sample);
        for input in inputs {
            command.arg(input);
        }
        command.arg("--reference").arg(&reference).assert()
    };
    call(SAMPLE_ID, &[&document, &fasta])
        .failure()
        .stderr(predicates::str::contains(
            "a consensus document is called on its own",
        ));
    call("sample-2", &[&document])
        .failure()
        .stderr(predicates::str::contains(
            "document sample \"sample-1\" is not the requested sample \"sample-2\"",
        ));
    call(SAMPLE_ID, &[&document]).success();
    let variants: Value = serde_json::from_slice(&fs::read(
        results.join(format!("{SAMPLE_ID}.variants.json")),
    )?)?;
    assert_eq!(variants["provenance"]["plugins"][0]["id"], "decided");
    assert_eq!(variants["reads"][0]["name"], format!("{SAMPLE_ID}_5-32"));
    Ok(())
}
