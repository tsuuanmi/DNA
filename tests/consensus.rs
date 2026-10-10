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
