//! CLI contract for a run context (PROP-0014): `dna run-context` over a
//! batch's read list, and `dna consensus --run-context` reading it.

pub mod support;

use std::fmt::Write as _;
use std::fs;
use std::path::{Path, PathBuf};

use assert_cmd::Command;
use serde_json::Value;
use tempfile::tempdir;

use support::{write_abif, write_config, write_reference};

const QUERY: &str = "ACGTCAGTACGATCGTACCTGAGTACGA";

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

/// A batch of two samples, each with a forward and a reverse read.
struct Batch {
    directory: tempfile::TempDir,
    reference: PathBuf,
    config: PathBuf,
    reads: Vec<(&'static str, PathBuf)>,
}

impl Batch {
    fn new() -> Result<Self, Box<dyn std::error::Error>> {
        let directory = tempdir()?;
        let reference = directory.path().join("reference.fa");
        let config = directory.path().join("dna.toml");
        write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
        write_config(&config, "linear")?;
        let mut reads = Vec::new();
        for sample in ["S1", "S2"] {
            for (strand, sequence) in [("F", QUERY.to_owned()), ("R", reverse_complement(QUERY))] {
                let path = directory.path().join(format!("{sample}-{strand}.ab1"));
                // A strand's reads differ between samples, so their contents
                // do too.
                let sequence = if sample == "S2" {
                    format!("{sequence}A")
                } else {
                    sequence
                };
                write_abif(&path, &sequence)?;
                reads.push((sample, path));
            }
        }
        Ok(Self {
            directory,
            reference,
            config,
            reads,
        })
    }

    fn workdir(&self) -> &Path {
        self.directory.path()
    }

    fn command(&self) -> Command {
        let mut command = Command::new(dna_binary());
        command
            .current_dir(self.workdir())
            .env("DNA_CONFIG", &self.config)
            .env("DNA_LOG_DIR", self.workdir().join("logs"));
        command
    }

    /// Writes a read list of `samples`, in the given order, and runs
    /// `dna run-context` over it.
    fn run_context(
        &self,
        run_id: &str,
        samples: &[&str],
    ) -> Result<assert_cmd::assert::Assert, Box<dyn std::error::Error>> {
        let list = self.workdir().join(format!("{run_id}.tsv"));
        let mut text = String::new();
        for sample in samples {
            for (name, path) in self.reads.iter().filter(|(name, _)| name == sample) {
                writeln!(text, "{name}\t{}", path.display())?;
            }
        }
        fs::write(&list, text)?;
        Ok(self
            .command()
            .args(["run-context", run_id, "--traces"])
            .arg(&list)
            .arg("--reference")
            .arg(&self.reference)
            .assert())
    }

    fn consensus(&self, sample: &str, context: &Path) -> assert_cmd::assert::Assert {
        let mut command = self.command();
        command
            .args(["consensus", sample, "--run-context"])
            .arg(context);
        for (_, path) in self.reads.iter().filter(|(name, _)| *name == sample) {
            command.arg(path);
        }
        command.arg("--reference").arg(&self.reference).assert()
    }

    fn document(&self, name: &str) -> Result<Value, Box<dyn std::error::Error>> {
        Ok(serde_json::from_slice(&fs::read(
            self.workdir().join("results").join(name),
        )?)?)
    }
}

/// The run context holds every read of the batch whatever the list's order,
/// and a sample of the batch is analysed with its own contribution left out.
// SRS-CONS-014 SRS-API-021
#[test]
fn writes_a_deterministic_run_context_and_leaves_the_sample_out()
-> Result<(), Box<dyn std::error::Error>> {
    let batch = Batch::new()?;
    batch.run_context("RUN1", &["S1", "S2"])?.success();
    let first = fs::read(batch.workdir().join("results/RUN1.run_context.json"))?;
    batch.run_context("RUN2", &["S2", "S1"])?.success();
    let value = batch.document("RUN1.run_context.json")?;
    let mut second = batch.document("RUN2.run_context.json")?;
    second["run_id"] = "RUN1".into();
    assert_eq!(serde_json::from_slice::<Value>(&first)?, second);
    assert_eq!(value["schema_version"], "dna.run_context/v1");
    assert_eq!(value["samples"], 2);
    assert_eq!(value["reads"].as_array().map(Vec::len), Some(4));
    assert_eq!(value["length"]["method"], 1);
    assert_eq!(value["length"]["runs"][0]["run"], 1);

    assert_eq!(value["provenance"]["skipped"], serde_json::json!([]));

    let context = batch.workdir().join("results/RUN1.run_context.json");
    batch.consensus("S1", &context).success();
    let used = &batch.document("S1.consensus.json")?["provenance"]["run_context"];
    assert_eq!(used["run_id"], "RUN1");
    assert_eq!(used["left_out"], true);
    // No anchored window, no pairs: τ stays the configured one.
    assert_eq!(used["length_scale"], serde_json::json!([null, null]));
    Ok(())
}

/// A sample outside the batch uses the context whole; a context that holds
/// only some of a sample's reads, or was made against another reference,
/// with another τ, or over other windows, is refused.
// SRS-CONS-014
#[test]
fn refuses_a_context_it_cannot_leave_the_sample_out_of() -> Result<(), Box<dyn std::error::Error>> {
    let batch = Batch::new()?;
    batch.run_context("RUN1", &["S1"])?.success();
    let context = batch.workdir().join("results/RUN1.run_context.json");
    batch.consensus("S2", &context).success();
    assert_eq!(
        batch.document("S2.consensus.json")?["provenance"]["run_context"]["left_out"],
        false
    );

    let mut partial = batch.document("RUN1.run_context.json")?;
    partial["reads"]
        .as_array_mut()
        .ok_or("reads must be an array")?
        .pop();
    let partial_path = batch.workdir().join("partial.json");
    fs::write(&partial_path, serde_json::to_vec(&partial)?)?;
    batch
        .consensus("S1", &partial_path)
        .failure()
        .stderr(predicates::str::contains("only partly"));

    let mut other = batch.document("RUN1.run_context.json")?;
    other["provenance"]["reference"]["sha256"] = "0".repeat(64).into();
    let other_path = batch.workdir().join("other.json");
    fs::write(&other_path, serde_json::to_vec(&other)?)?;
    batch
        .consensus("S1", &other_path)
        .failure()
        .stderr(predicates::str::contains("made against reference"));

    for (field, value) in [
        ("noise_forward", serde_json::json!([0.2, 0.111])),
        ("windows", serde_json::json!([[303, 315]])),
    ] {
        let mut other = batch.document("RUN1.run_context.json")?;
        other["length"][field] = value;
        let path = batch.workdir().join(format!("{field}.json"));
        fs::write(&path, serde_json::to_vec(&other)?)?;
        batch
            .consensus("S1", &path)
            .failure()
            .stderr(predicates::str::contains("noise, bin edges, or windows"));
    }
    Ok(())
}

/// A listed file every sample command fails on is left out and recorded;
/// the sample's other reads still count, so the sample, analysed without
/// it, has its own pairs left out.
// SRS-API-021
#[test]
fn leaves_out_an_empty_trace_and_records_it() -> Result<(), Box<dyn std::error::Error>> {
    let mut batch = Batch::new()?;
    let empty = batch.workdir().join("S1-empty.ab1");
    fs::write(&empty, b"")?;
    batch.reads.push(("S1", empty));
    batch.run_context("RUN1", &["S1", "S2"])?.success();
    let value = batch.document("RUN1.run_context.json")?;
    assert_eq!(
        value["provenance"]["skipped"],
        serde_json::json!([{"sample_id": "S1", "name": "S1-empty.ab1", "reason": "empty"}])
    );
    assert_eq!(value["reads"].as_array().map(Vec::len), Some(4));
    batch.reads.pop();
    let context = batch.workdir().join("results/RUN1.run_context.json");
    batch.consensus("S1", &context).success();
    assert_eq!(
        batch.document("S1.consensus.json")?["provenance"]["run_context"]["left_out"],
        true
    );
    Ok(())
}
