//! Public `SampleAnalysis::call`: the `consensus` then `call` workflow in
//! process, which calls the consensus's decided edits (PROP-0006).

pub mod support;

use std::fs;
use std::path::{Path, PathBuf};

use assert_cmd::Command;
use serde_json::Value;
use tempfile::tempdir;

use std::num::NonZeroUsize;

use dna::sample_analysis;
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

/// Writes the synthetic configuration with a profile that declares notation.
fn write_config_with_notation(path: &Path) -> Result<(), Box<dyn std::error::Error>> {
    write_config(path, "linear")?;
    let profile = path.with_extension("profile.toml");
    let text = fs::read_to_string(&profile)?;
    fs::write(
        &profile,
        format!(
            "{text}[normalization]\nindel_placement='right'\nmax_trade=20\n[notation]\nstyle='per_base_decimal'\n"
        ),
    )?;
    Ok(())
}

/// Runs one CLI command in `workdir` and returns its published JSON document.
fn cli(
    workdir: &Path,
    config: &Path,
    arguments: &[&Path],
    command: &str,
    document: &str,
) -> Result<Value, Box<dyn std::error::Error>> {
    let mut process = Command::new(dna_binary());
    process
        .current_dir(workdir)
        .env("DNA_CONFIG", config)
        .env("DNA_LOG_DIR", workdir.join("logs"))
        .arg(command)
        .arg(SAMPLE_ID);
    for argument in arguments {
        process.arg(argument);
    }
    process.assert().success();
    Ok(serde_json::from_slice(&fs::read(
        workdir.join("results").join(document),
    )?)?)
}

#[test]
fn calls_a_sample_like_consensus_then_call_without_side_effects()
-> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let forward = directory.path().join("read-forward.ab1");
    let reverse = directory.path().join("read-reverse.ab1");
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config_with_notation(&config)?;
    write_abif(&forward, &alternate())?;
    write_abif(&reverse, &reverse_complement(&alternate()))?;
    let traces = [forward.clone(), reverse.clone()];

    let analysis =
        sample_analysis::analyze_sanger(&traces, &reference, &config, NonZeroUsize::MIN)?;
    let result = analysis.call(SAMPLE_ID)?;

    assert_eq!(analysis.reference.name, "synthetic");
    assert_eq!(analysis.profile.id, "synthetic-linear");
    assert!(analysis.rejected.is_empty());
    let notation = result
        .notation
        .clone()
        .ok_or("the profile declares notation")?;
    assert_eq!(notation.len(), 1);
    let call = &notation[0];
    assert_eq!((call.position, call.ordinal), (15, 0));
    assert_eq!((call.reference, call.alternate), (Some('G'), Some('A')));
    assert_eq!(call.text, "15A");
    assert!(!result.callable_segments.is_empty());
    assert!(!directory.path().join("results").exists());
    assert!(!directory.path().join("logs").exists());

    // The CLI pair, `call` given the consensus document, publishes the same
    // calls and identities.
    let workdir = directory.path().join("cli");
    fs::create_dir(&workdir)?;
    let consensus = cli(
        &workdir,
        &config,
        &[&forward, &reverse, Path::new("--reference"), &reference],
        "consensus",
        &format!("{SAMPLE_ID}.consensus.json"),
    )?;
    assert_eq!(consensus["rejected_reads"], serde_json::json!([]));
    let document: PathBuf = workdir
        .join("results")
        .join(format!("{SAMPLE_ID}.consensus.json"));
    let variants = cli(
        &workdir,
        &config,
        &[&document, Path::new("--reference"), &reference],
        "call",
        &format!("{SAMPLE_ID}.variants.json"),
    )?;
    assert_eq!(
        variants["provenance"]["plugins"]
            .as_array()
            .map(|plugins| plugins
                .iter()
                .map(|plugin| plugin["id"].clone())
                .collect::<Vec<_>>()),
        Some(vec![
            serde_json::json!("decided"),
            serde_json::json!("normalization"),
            serde_json::json!("nomenclature"),
        ])
    );
    let calls = variants["notation"]["calls"]
        .as_array()
        .ok_or("notation calls must be an array")?
        .iter()
        .map(|call| call["call"].as_str().map(str::to_owned))
        .collect::<Option<Vec<_>>>();
    let texts = notation.iter().map(|call| call.text.clone()).collect();
    assert_eq!(calls, Some(texts));
    let segments = variants["reads"]
        .as_array()
        .ok_or("reads must be an array")?
        .iter()
        .flat_map(|read| {
            read["alignment"]["callable_reference_segments"]
                .as_array()
                .cloned()
                .unwrap_or_default()
        })
        .map(|segment| (segment["start"].as_u64(), segment["end"].as_u64()))
        .collect::<Vec<_>>();
    let expected = result
        .callable_segments
        .iter()
        .map(|segment| {
            (
                u64::try_from(segment.start_0based).ok(),
                u64::try_from(segment.end_0based_exclusive).ok(),
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(segments, expected);
    assert_eq!(
        variants["provenance"]["configuration_sha256"],
        analysis.evidence.configuration_sha256.as_str()
    );
    assert_eq!(
        variants["provenance"]["profile"]["sha256"],
        analysis.profile.sha256.as_str()
    );

    // The consensus FASTA still takes the sequence path: the record is placed
    // again, and here it gives the same call.
    let sequences = directory.path().join("sequences");
    fs::create_dir(&sequences)?;
    let fasta: PathBuf = workdir
        .join("results")
        .join(format!("{SAMPLE_ID}.consensus.fasta"));
    let placed = cli(
        &sequences,
        &config,
        &[&fasta, Path::new("--reference"), &reference],
        "call",
        &format!("{SAMPLE_ID}.variants.json"),
    )?;
    assert_eq!(placed["provenance"]["plugins"][0]["id"], "sequence");
    assert_eq!(placed["notation"]["calls"], variants["notation"]["calls"]);
    Ok(())
}

#[test]
fn validates_the_sample_identifier_and_follows_the_profile()
-> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let trace = directory.path().join("read.ab1");
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;
    write_abif(&trace, QUERY)?;

    let analysis = sample_analysis::analyze_sanger(
        std::slice::from_ref(&trace),
        &reference,
        &config,
        NonZeroUsize::MIN,
    )?;
    assert!(analysis.call("../x").is_err());
    // Without notation in the profile there are callable segments but no
    // notation calls.
    let call = analysis.call(SAMPLE_ID)?;
    assert!(call.notation.is_none());
    assert!(!call.callable_segments.is_empty());
    assert!(!directory.path().join("results").exists());
    Ok(())
}
