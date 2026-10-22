//! CLI contract for the review of the stretches a call does not analyse
//! (SRS-CONS-012, review rules 36): every such stretch is one `not_analysed`
//! check item and one `not_analysed` zone, so an allele omitted there is
//! never passed silently.
//!
//! The reads are analysed from their recorded evidence
//! (`dna.read_evidence/v1`), with the modality's unreliable spans set where
//! each case needs them.

pub mod support;

use std::fs;
use std::path::{Path, PathBuf};

use assert_cmd::Command;
use serde_json::{Value, json};
use tempfile::tempdir;

use support::{write_abif_with_peak_heights, write_config_with_profile, write_reference};

const QUERY: &str = "ACGTCAGTACGATCGTACCTGAGTACGA";
const SAMPLE_ID: &str = "sample-1";
/// The reportable region: the reference positions `QUERY` covers.
const FIRST: usize = 5;
const LAST: usize = 32;

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

/// `QUERY` with the base at reference `position` changed to one that
/// equals neither it nor a neighbour, so the change makes no run.
fn with_change(position: usize) -> String {
    let mut bases = QUERY.as_bytes().to_vec();
    let index = position - FIRST;
    let near = |base: u8| {
        bases[index] == base
            || index.checked_sub(1).and_then(|before| bases.get(before)) == Some(&base)
            || bases.get(index + 1) == Some(&base)
    };
    let changed = b"ACGT".iter().copied().find(|&base| !near(base));
    bases[index] = changed.unwrap_or(b'N');
    String::from_utf8_lossy(&bases).into_owned()
}

/// The reference, and a configuration whose profile reports `FIRST..=LAST`.
fn inputs(directory: &Path) -> Result<(PathBuf, PathBuf), Box<dyn std::error::Error>> {
    let reference = directory.join("reference.fa");
    let config = directory.join("dna.toml");
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    let profile = directory.join("dna.profile.toml");
    fs::write(
        &profile,
        format!(
            "schema_version=2\nid='synthetic-region'\n[reference]\ntopology='linear'\n\
             [variant_calling]\nregions=[[{FIRST}, {LAST}]]\n"
        ),
    )?;
    write_config_with_profile(&config, Path::new("dna.profile.toml"))?;
    Ok((reference, config))
}

/// Writes `reads` (sequence, forward) as traces, records each one's
/// evidence with `unreliable(call_count, forward)` as its unreliable spans
/// (0-based half-open call intervals), and returns the evidence files.
fn evidence(
    directory: &Path,
    config: &Path,
    reads: &[(&str, bool)],
    unreliable: &dyn Fn(usize, bool) -> Vec<(usize, usize)>,
) -> Result<Vec<PathBuf>, Box<dyn std::error::Error>> {
    let mut files = Vec::new();
    for (index, &(sequence, forward)) in reads.iter().enumerate() {
        let name = format!("read-{index}");
        let trace = directory.join(format!("{name}.ab1"));
        let read = if forward {
            sequence.to_owned()
        } else {
            reverse_complement(sequence)
        };
        // Distinct heights give every read its own content.
        let height = 800 + 10 * i16::try_from(index)?;
        write_abif_with_peak_heights(&trace, &read, vec![height; read.len()])?;
        Command::new(dna_binary())
            .current_dir(directory)
            .env("DNA_CONFIG", config)
            .env("DNA_LOG_DIR", directory.join("logs"))
            .arg("evidence")
            .arg(&trace)
            .assert()
            .success();
        let file = directory
            .join("results")
            .join(format!("{name}.evidence.json"));
        let mut value: Value = serde_json::from_slice(&fs::read(&file)?)?;
        let calls = value["read"]["call_count"]
            .as_u64()
            .ok_or("evidence lacks its call count")?;
        value["read"]["unreliable"] = Value::Array(
            unreliable(usize::try_from(calls)?, forward)
                .into_iter()
                .map(|(start, end)| {
                    json!({"calls": {"start": start, "end": end}, "reason": "poor_start_signal"})
                })
                .collect(),
        );
        fs::write(&file, serde_json::to_vec(&value)?)?;
        files.push(file);
    }
    Ok(files)
}

/// `dna review --evidence` over `files`; returns the review document.
fn review(
    directory: &Path,
    reference: &Path,
    config: &Path,
    files: &[PathBuf],
) -> Result<Value, Box<dyn std::error::Error>> {
    let mut command = Command::new(dna_binary());
    command
        .current_dir(directory)
        .env("DNA_CONFIG", config)
        .env("DNA_LOG_DIR", directory.join("logs"))
        .arg("review")
        .arg(SAMPLE_ID)
        .arg("--evidence");
    for file in files {
        command.arg(file);
    }
    command.arg("--reference").arg(reference).assert().success();
    Ok(serde_json::from_slice(&fs::read(
        directory
            .join("results")
            .join(format!("{SAMPLE_ID}.review.json")),
    )?)?)
}

/// The `not_analysed` items: first, last, kind, and the reads giving
/// evidence.
fn stretch_items(value: &Value) -> Vec<(u64, u64, String, usize)> {
    value["items"]
        .as_array()
        .into_iter()
        .flatten()
        .filter(|item| {
            item["labels"]
                .as_array()
                .into_iter()
                .flatten()
                .any(|label| label["label"] == "not_analysed")
        })
        .map(|item| {
            (
                item["first"]["position"].as_u64().unwrap_or(0),
                item["last"]["position"].as_u64().unwrap_or(0),
                item["kind"].as_str().unwrap_or("").to_owned(),
                item["evidence"].as_array().map_or(0, Vec::len),
            )
        })
        .collect()
}

/// The zones: kind, first, last, and checks.
fn zones(value: &Value) -> Vec<(String, u64, u64, u64)> {
    value["zones"]
        .as_array()
        .into_iter()
        .flatten()
        .map(|zone| {
            (
                zone["kind"].as_str().unwrap_or("").to_owned(),
                zone["first"].as_u64().unwrap_or(0),
                zone["last"].as_u64().unwrap_or(0),
                zone["checks"].as_u64().unwrap_or(0),
            )
        })
        .collect()
}

/// Three reads, two forward and one reverse, each with `unreliable` spans.
fn three_reads(
    sequence: &str,
    unreliable: &dyn Fn(usize, bool) -> Vec<(usize, usize)>,
) -> Result<Value, Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let directory = directory.path();
    let (reference, config) = inputs(directory)?;
    let files = evidence(
        directory,
        &config,
        &[(sequence, true), (sequence, true), (sequence, false)],
        unreliable,
    )?;
    review(directory, &reference, &config, &files)
}

/// The spans of the first `n` reference positions on every read: a forward
/// read's first calls, a reverse read's last.
fn at_start(n: usize) -> impl Fn(usize, bool) -> Vec<(usize, usize)> {
    move |calls, forward| {
        if forward {
            vec![(0, n)]
        } else {
            vec![(calls - n, calls)]
        }
    }
}

/// The spans of the last `n` reference positions on every read.
fn at_end(n: usize) -> impl Fn(usize, bool) -> Vec<(usize, usize)> {
    move |calls, forward| {
        if forward {
            vec![(calls - n, calls)]
        } else {
            vec![(0, n)]
        }
    }
}

// SRS-CONS-012 SRS-API-014
/// A difference inside the poor start of every read: the region's start up
/// to it is not analysed, and is one check item with every covering read's
/// evidence and one zone of its own.
#[test]
fn a_stretch_at_a_region_start_is_a_check_and_a_zone() -> Result<(), Box<dyn std::error::Error>> {
    let value = three_reads(&with_change(7), &at_start(5))?;
    assert_eq!(value["not_analysed"][0]["first"], 5);
    assert_eq!(value["not_analysed"][0]["last"], 7);
    assert_eq!(
        stretch_items(&value),
        [(5, 7, "check".to_owned(), 3)],
        "{value}"
    );
    let zones = zones(&value);
    assert_eq!(zones[0], ("not_analysed".to_owned(), 5, 7, 1));
    assert!(zones[1..].iter().all(|zone| zone.0 != "not_analysed"));
    assert_eq!(value["summary"]["clear"], false);
    assert_eq!(value["summary"]["labels"]["not_analysed"], 1);
    assert_eq!(value["summary"]["positions"]["not_analysed"], 3);
    assert_eq!(value["summary"]["positions"]["check"], 0);
    Ok(())
}

// SRS-CONS-012
/// The same at a region's end.
#[test]
fn a_stretch_at_a_region_end_is_a_check_and_a_zone() -> Result<(), Box<dyn std::error::Error>> {
    let value = three_reads(&with_change(30), &at_end(5))?;
    assert_eq!(
        stretch_items(&value),
        [(30, 32, "check".to_owned(), 3)],
        "{value}"
    );
    let zones = zones(&value);
    assert_eq!(zones.last(), Some(&("not_analysed".to_owned(), 30, 32, 1)));
    Ok(())
}

// SRS-CONS-012
/// A region whose every position only poor spans read, with a difference at
/// its last position, is not analysed whole: one item, one zone, and no
/// other zone in it.
#[test]
fn a_region_not_analysed_whole_is_one_check_and_one_zone() -> Result<(), Box<dyn std::error::Error>>
{
    let value = three_reads(&with_change(LAST), &|calls, _| vec![(0, calls)])?;
    assert_eq!(
        stretch_items(&value),
        [(5, 32, "check".to_owned(), 3)],
        "{value}"
    );
    assert_eq!(zones(&value), [("not_analysed".to_owned(), 5, 32, 1)]);
    assert_eq!(value["callable"], json!([]));
    Ok(())
}

// SRS-CONS-012
/// A poor start that reads the reference throughout leaves the region whole:
/// no stretch, no item, and the sample is clear.
#[test]
fn a_poor_start_reading_the_reference_is_no_item() -> Result<(), Box<dyn std::error::Error>> {
    let value = three_reads(QUERY, &at_start(5))?;
    assert_eq!(value["not_analysed"], json!([]));
    assert!(stretch_items(&value).is_empty());
    assert!(zones(&value).iter().all(|zone| zone.0 != "not_analysed"));
    assert_eq!(value["summary"]["clear"], true);
    Ok(())
}
