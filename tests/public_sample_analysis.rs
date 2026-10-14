//! Public `sample_analysis` capability: one placement pass yields the sample
//! evidence and the consensus, which match the `sample` and `consensus`
//! commands.

pub mod support;

use std::fs;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

use assert_cmd::Command;
use serde_json::Value;
use tempfile::tempdir;

use dna::sample_analysis::{self, RejectionCause};
use support::{write_abif, write_config, write_reference};

const QUERY: &str = "ACGTCAGTACGATCGTACCTGAGTACGA";
const UNPLACEABLE: &str = "AGCATGAGTCCATGCTAGCATGACTGCA";
const SAMPLE_ID: &str = "sample-1";
const FOUR: NonZeroUsize = match NonZeroUsize::new(4) {
    Some(workers) => workers,
    None => NonZeroUsize::MIN,
};

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

/// Two agreeing reads carrying `15A` and one read the core cannot place.
fn write_sample(directory: &Path) -> Result<[PathBuf; 5], Box<dyn std::error::Error>> {
    let reference = directory.join("reference.fa");
    let config = directory.join("dna.toml");
    let forward = directory.join("read-forward.ab1");
    let foreign = directory.join("read-foreign.ab1");
    let reverse = directory.join("read-reverse.ab1");
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;
    write_abif(&forward, &alternate())?;
    write_abif(&foreign, UNPLACEABLE)?;
    write_abif(&reverse, &reverse_complement(&alternate()))?;
    Ok([reference, config, forward, foreign, reverse])
}

fn run_command(
    command: &str,
    traces: &[PathBuf],
    reference: &Path,
    config: &Path,
    workdir: &Path,
) -> Result<Value, Box<dyn std::error::Error>> {
    let binary = std::env::var("CARGO_BIN_EXE_dna")?;
    Command::new(binary)
        .current_dir(workdir)
        .env("DNA_CONFIG", config)
        .arg(command)
        .arg(SAMPLE_ID)
        .args(traces)
        .arg("--reference")
        .arg(reference)
        .assert()
        .success();
    let name = if command == "sample" {
        format!("{SAMPLE_ID}.sample.json")
    } else {
        format!("{SAMPLE_ID}.{command}.json")
    };
    Ok(serde_json::from_slice(&fs::read(
        workdir.join("results").join(name),
    )?)?)
}

#[test]
fn one_pass_returns_sample_evidence_and_consensus_without_cli_side_effects()
-> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let [reference, config, forward, foreign, reverse] = write_sample(directory.path())?;

    let result =
        sample_analysis::analyze_sanger(&[forward, foreign, reverse], &reference, &config, FOUR)?;

    assert_eq!(result.reference.name, "synthetic");
    assert_eq!(result.evidence.reference_sha256, result.reference.sha256);
    assert!(!result.evidence.configuration_sha256.is_empty());
    assert_eq!(result.evidence.reads.len(), 2);
    assert_eq!(result.reads.len(), 2);
    assert!(
        result
            .evidence
            .reads
            .iter()
            .all(|read| result.reads.contains_key(&read.input_sha256))
    );
    assert_eq!(result.rejected.len(), 1);
    let rejected = result.rejected.values().next().ok_or("one rejected read")?;
    assert_eq!(rejected.sanger.input_name, "read-foreign.ab1");
    assert!(matches!(rejected.cause, RejectionCause::Placement(_)));
    assert_eq!(result.consensus.segments.len(), 1);
    assert_eq!(result.consensus.segments[0].sequence, alternate());
    assert_eq!(result.consensus.sites.len(), 1);

    assert!(!directory.path().join("results").exists());
    assert!(!directory.path().join("logs").exists());
    Ok(())
}

/// The library pass and the two commands describe the same reads, variants,
/// and consensus.
#[test]
fn one_pass_matches_the_sample_and_consensus_commands() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let workdir = directory.path().join("work");
    fs::create_dir(&workdir)?;
    let [reference, config, forward, foreign, reverse] = write_sample(directory.path())?;
    let traces = [forward, foreign, reverse];

    let result = sample_analysis::analyze_sanger(&traces, &reference, &config, FOUR)?;
    let sample = run_command("sample", &traces, &reference, &config, &workdir)?;
    let consensus = run_command("consensus", &traces, &reference, &config, &workdir)?;

    let names = |value: &Value| -> Vec<String> {
        value
            .as_array()
            .map(|reads| {
                reads
                    .iter()
                    .filter_map(|read| read["sha256"].as_str().map(str::to_owned))
                    .collect()
            })
            .unwrap_or_default()
    };
    let admitted = result
        .evidence
        .reads
        .iter()
        .map(|read| read.input_sha256.clone())
        .collect::<Vec<_>>();
    assert_eq!(names(&sample["reads"]), admitted);
    assert_eq!(names(&consensus["reads"]), admitted);
    let rejected = result.rejected.keys().cloned().collect::<Vec<_>>();
    assert_eq!(names(&sample["rejected_reads"]), rejected);
    assert_eq!(names(&consensus["rejected_reads"]), rejected);
    assert_eq!(
        sample["variants"].as_array().map(Vec::len),
        Some(result.evidence.variants.len())
    );
    assert_eq!(
        consensus["segments"][0]["sequence"],
        result.consensus.segments[0].sequence.as_str()
    );
    assert_eq!(
        consensus["summary"]["called_positions"],
        result.consensus.summary.called_positions
    );
    assert_eq!(
        consensus["sites"].as_array().map(Vec::len),
        Some(result.consensus.sites.len())
    );
    Ok(())
}

/// One worker and several give the same result.
#[test]
fn the_result_does_not_depend_on_the_number_of_workers() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let [reference, config, forward, foreign, reverse] = write_sample(directory.path())?;
    let traces = [forward, foreign, reverse];

    let one = sample_analysis::analyze_sanger(&traces, &reference, &config, NonZeroUsize::MIN)?;
    let many = sample_analysis::analyze_sanger(&traces, &reference, &config, FOUR)?;

    assert_eq!(format!("{one:?}"), format!("{many:?}"));
    Ok(())
}
