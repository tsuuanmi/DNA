//! CLI contract for the sample review list (PROP-0004, SRS-API-014).

pub mod support;

use std::fs;
use std::path::{Path, PathBuf};

use assert_cmd::Command;
use serde_json::Value;
use tempfile::tempdir;

use support::{
    write_abif, write_abif_with_incoherent_doubles, write_abif_with_secondary_signal, write_config,
    write_reference,
};

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

/// The synthetic reference and a configuration whose profile declares
/// notation, so the review list names the sample's calls.
fn inputs(directory: &Path) -> Result<(PathBuf, PathBuf), Box<dyn std::error::Error>> {
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
    Ok((reference, config))
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
        .arg("review")
        .arg(SAMPLE_ID);
    for trace in traces {
        command.arg(trace);
    }
    command.arg("--reference").arg(reference).assert()
}

/// `dna review` judged against the run signal `signal`.
fn run_with_signal(
    traces: &[&PathBuf],
    reference: &Path,
    config: &Path,
    workdir: &Path,
    signal: &Path,
) -> assert_cmd::assert::Assert {
    let mut command = Command::new(dna_binary());
    command
        .current_dir(workdir)
        .env("DNA_CONFIG", config)
        .arg("review")
        .arg(SAMPLE_ID);
    for trace in traces {
        command.arg(trace);
    }
    command
        .arg("--reference")
        .arg(reference)
        .arg("--signal")
        .arg(signal)
        .assert()
}

/// `dna signal <run-id>` over `traces`, in `workdir`; returns its output.
fn signal(run_id: &str, traces: &[&PathBuf], config: &Path, workdir: &Path) -> PathBuf {
    let mut command = Command::new(dna_binary());
    command
        .current_dir(workdir)
        .env("DNA_CONFIG", config)
        .arg("signal")
        .arg(run_id);
    for trace in traces {
        command.arg(trace);
    }
    command.assert().success();
    workdir
        .join("results")
        .join(format!("{run_id}.signal.json"))
}

fn output(workdir: &Path) -> PathBuf {
    workdir
        .join("results")
        .join(format!("{SAMPLE_ID}.review.json"))
}

fn document(workdir: &Path) -> Result<Value, Box<dyn std::error::Error>> {
    Ok(serde_json::from_slice(&fs::read(output(workdir))?)?)
}

fn items(value: &Value) -> Result<&Vec<Value>, Box<dyn std::error::Error>> {
    Ok(value["items"].as_array().ok_or("items must be an array")?)
}

fn label_names(item: &Value) -> Vec<&str> {
    item["labels"]
        .as_array()
        .map(|labels| {
            labels
                .iter()
                .filter_map(|label| label["label"].as_str())
                .collect()
        })
        .unwrap_or_default()
}

const RANKS: [(&str, u64); 15] = [
    ("strands_disagree", 0),
    ("read_disagrees", 0),
    ("possible_heteroplasmy", 1),
    ("estimated_length", 2),
    ("undecided", 2),
    ("not_called", 2),
    ("long_run", 2),
    ("edge_deletion", 2),
    ("one_read", 3),
    ("one_strand", 3),
    ("noisy_signal", 4),
    ("not_covered", 4),
    ("phase_shift", 4),
    ("dephased_end", 4),
    ("read_edge", 4),
];

/// The document's structural invariants: items ordered by priority then
/// position, kinds that follow the calls, and a summary that counts them.
fn assert_consistent(value: &Value) -> Result<(), Box<dyn std::error::Error>> {
    assert_eq!(value["schema_version"], "dna.review/v1");
    assert_eq!(value["provenance"]["rules"], 10);
    // An admitted read's masked stretches are ordered, disjoint, and well formed.
    for read in value["reads"].as_array().ok_or("reads must be an array")? {
        if read["admitted"] != true {
            assert!(read.get("masked").is_none());
            continue;
        }
        let mut end = 0;
        for span in read["masked"].as_array().ok_or("masked must be an array")? {
            let first = span["first"].as_u64().ok_or("masked first")?;
            let last = span["last"].as_u64().ok_or("masked last")?;
            assert!(end < first && first <= last, "{span}");
            end = last;
        }
    }
    let items = items(value)?;
    let mut previous = None;
    let mut counts = std::collections::BTreeMap::<&str, u64>::new();
    let mut checks = 0;
    for item in items {
        let names = label_names(item);
        let priority = names
            .iter()
            .filter_map(|name| RANKS.iter().find(|(label, _)| label == name))
            .map(|(_, rank)| *rank)
            .min()
            .unwrap_or(5);
        let position = item["first"]["position"].as_u64().ok_or("first position")?;
        let order = (priority, position);
        assert!(
            previous.is_none_or(|previous| previous <= order),
            "{order:?}"
        );
        previous = Some(order);
        let calls = item["calls"].as_array().ok_or("calls must be an array")?;
        let kind = if calls.is_empty() { "check" } else { "variant" };
        assert_eq!(item["kind"], kind);
        if kind == "check" {
            checks += 1;
            assert!(!names.is_empty());
        }
        for name in names.iter().collect::<std::collections::BTreeSet<_>>() {
            *counts.entry(name).or_default() += 1;
        }
    }
    let summary = &value["summary"];
    assert_eq!(summary["items"], items.len());
    assert_eq!(summary["checks"], checks);
    assert_eq!(summary["variants"], items.len() - checks);
    let labels = summary["labels"]
        .as_object()
        .ok_or("summary labels must be an object")?;
    assert_eq!(
        labels
            .keys()
            .map(String::as_str)
            .collect::<std::collections::BTreeSet<_>>(),
        RANKS.iter().map(|(name, _)| *name).collect()
    );
    for (name, count) in labels {
        assert_eq!(
            count.as_u64(),
            Some(counts.get(name.as_str()).copied().unwrap_or(0))
        );
    }
    let held = items
        .iter()
        .filter(|item| item["kind"] == "check")
        .flat_map(|item| {
            let first = item["first"]["position"].as_u64().unwrap_or(1);
            let last = item["last"]["position"].as_u64().unwrap_or(0);
            first..=last
        })
        .collect::<std::collections::BTreeSet<_>>();
    assert_eq!(summary["positions"]["check"], held.len());
    assert!(summary["positions"]["reportable"].as_u64() >= Some(held.len() as u64));
    assert_eq!(summary["clear"], value["status"] == "called" && checks == 0);
    Ok(())
}

/// Two reads that agree on a substitution give a variant item with each
/// read's evidence, and the reference ends no read covers are `not_covered`.
#[test]
fn lists_the_items_of_a_called_sample() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let (reference, config) = inputs(directory.path())?;
    let forward = directory.path().join("read-forward.ab1");
    let reverse = directory.path().join("read-reverse.ab1");
    write_abif(&forward, &alternate())?;
    write_abif(&reverse, &reverse_complement(&alternate()))?;

    run(&[&forward, &reverse], &reference, &config, directory.path()).success();

    let value = document(directory.path())?;
    assert_consistent(&value)?;
    assert_eq!(value["status"], "called");
    // The summary lists every label, in priority order.
    let text = fs::read_to_string(output(directory.path()))?;
    let summary = text.rfind("\"summary\"").ok_or("summary must be present")?;
    let offsets = RANKS
        .iter()
        .map(|(name, _)| text[summary..].find(&format!("\"{name}\":")))
        .collect::<Option<Vec<_>>>()
        .ok_or("every label must be counted")?;
    assert!(offsets.is_sorted());
    assert_eq!(
        value["provenance"]["plugins"]
            .as_array()
            .ok_or("plugins must be an array")?
            .iter()
            .map(|plugin| plugin["id"].as_str().unwrap_or_default())
            .collect::<Vec<_>>(),
        [
            "sanger",
            "core",
            "consensus",
            "decided",
            "normalization",
            "nomenclature",
            "review"
        ]
    );
    assert_eq!(
        value["reads"],
        serde_json::json!([
            {"name": "read-forward", "sha256": value["reads"][0]["sha256"], "strand": "forward", "admitted": true, "masked": []},
            {"name": "read-reverse", "sha256": value["reads"][1]["sha256"], "strand": "reverse", "admitted": true, "masked": []},
        ])
    );
    assert_eq!(value["calls"][0]["call"], "15A");
    assert_eq!(value["calls"][0]["position"], 15);
    let items = items(&value)?;
    let variant = items
        .iter()
        .find(|item| item["kind"] == "variant")
        .ok_or("the substitution must be an item")?;
    assert_eq!(variant["calls"], serde_json::json!(["15A"]));
    assert_eq!(variant["decision"], "called");
    assert_eq!(variant["sequence"], "AA");
    assert_eq!(
        variant["first"],
        serde_json::json!({"position": 15, "ordinal": 0})
    );
    assert_eq!(
        variant["last"],
        serde_json::json!({"position": 16, "ordinal": 0})
    );
    let evidence = variant["evidence"]
        .as_array()
        .ok_or("evidence must be an array")?;
    assert_eq!(evidence.len(), 2);
    let position = &evidence[0]["positions"][0];
    assert_eq!(position["position"], 15);
    assert_eq!(position["base"], "A");
    assert_eq!(position["state"], "alternate");
    assert!(position["scan"].is_u64() && position["quality"].is_u64());
    let uncovered = items
        .iter()
        .filter(|item| label_names(item) == ["not_covered"])
        .map(|item| (&item["first"]["position"], &item["last"]["position"]))
        .collect::<Vec<_>>();
    assert_eq!(
        uncovered,
        [(&1.into(), &4.into()), (&33.into(), &36.into())]
    );
    Ok(())
}

/// A read against a read on the other strand is a `strands_disagree` check,
/// listed first.
#[test]
fn reads_that_disagree_across_strands_come_first() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let (reference, config) = inputs(directory.path())?;
    let forward = directory.path().join("read-forward.ab1");
    let reverse = directory.path().join("read-reverse.ab1");
    write_abif(&forward, &alternate())?;
    write_abif(&reverse, &reverse_complement(QUERY))?;

    run(&[&forward, &reverse], &reference, &config, directory.path()).success();

    let value = document(directory.path())?;
    assert_consistent(&value)?;
    let first = &items(&value)?[0];
    assert_eq!(first["kind"], "check");
    assert_eq!(label_names(first)[0], "strands_disagree");
    assert_eq!(
        first["labels"][0]["reads"],
        serde_json::json!([
            {"read": "read-reverse", "strand": "reverse", "role": "supporting"},
            {"read": "read-forward", "strand": "forward", "role": "opposing"},
        ])
    );
    assert_eq!(value["summary"]["clear"], false);
    Ok(())
}

/// A second base on both strands at one position, beyond the reads'
/// neighbours, is a `possible_heteroplasmy` check with each read's share.
#[test]
fn a_second_base_on_both_strands_is_possible_heteroplasmy() -> Result<(), Box<dyn std::error::Error>>
{
    let directory = tempdir()?;
    let (reference, config) = inputs(directory.path())?;
    let forward = directory.path().join("read-forward.ab1");
    let reverse = directory.path().join("read-reverse.ab1");
    // Reference position 15 is G (call 10 of the forward read); its
    // neighbours are A and A, so a T there is no shadow.
    write_abif_with_secondary_signal(&forward, QUERY, 10, b'T', 400)?;
    let reverse_query = reverse_complement(QUERY);
    let call = QUERY.len() - 1 - 10;
    write_abif_with_secondary_signal(&reverse, &reverse_query, call, b'A', 400)?;

    run(&[&forward, &reverse], &reference, &config, directory.path()).success();

    let value = document(directory.path())?;
    assert_consistent(&value)?;
    let first = &items(&value)?[0];
    assert_eq!(first["first"]["position"], 15);
    assert_eq!(label_names(first), ["possible_heteroplasmy"]);
    let label = &first["labels"][0];
    assert_eq!(label["base"], "T");
    let reads = label["reads"].as_array().ok_or("reads must be an array")?;
    assert_eq!(reads.len(), 2);
    // Without a run signal the configured floor alone is the threshold.
    assert!(reads.iter().all(|read| {
        read["threshold"].as_f64() == Some(0.21)
            && read["share"].as_f64() > read["threshold"].as_f64()
            && read.get("calls").is_none()
    }));
    assert!(value["provenance"].get("signal").is_none());
    assert_eq!(first["sequence"], "G");
    assert_eq!(first["evidence"][0]["positions"][0]["second"]["base"], "T");
    Ok(())
}

/// A sample whose every read is rejected still gets its document, with the
/// rejected reads and their causes and no items.
#[test]
fn a_sample_without_an_admitted_read_lists_its_rejected_reads()
-> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let (reference, config) = inputs(directory.path())?;
    let first = directory.path().join("read-first.ab1");
    let second = directory.path().join("read-second.ab1");
    write_abif_with_incoherent_doubles(&first, QUERY, 0..QUERY.len(), 500)?;
    write_abif_with_incoherent_doubles(&second, QUERY, 0..QUERY.len(), 400)?;

    run(&[&first, &second], &reference, &config, directory.path()).success();

    let value = document(directory.path())?;
    assert_consistent(&value)?;
    assert_eq!(value["status"], "rejected");
    assert_eq!(value["items"], serde_json::json!([]));
    assert!(value.get("calls").is_none() && value.get("callable").is_none());
    let reads = value["reads"].as_array().ok_or("reads must be an array")?;
    assert_eq!(reads.len(), 2);
    for read in reads {
        assert_eq!(read["admitted"], false);
        assert_eq!(read["rejection"]["reason"], "callable_calls_below_minimum");
    }
    Ok(())
}

/// The same inputs give byte-identical documents, and an existing document
/// is never overwritten.
#[test]
fn reruns_are_identical_and_never_overwrite() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let (reference, config) = inputs(directory.path())?;
    let forward = directory.path().join("read-forward.ab1");
    let reverse = directory.path().join("read-reverse.ab1");
    write_abif(&forward, &alternate())?;
    write_abif(&reverse, &reverse_complement(QUERY))?;
    let mut documents = Vec::new();
    for name in ["one", "two"] {
        let workdir = directory.path().join(name);
        fs::create_dir(&workdir)?;
        run(&[&forward, &reverse], &reference, &config, &workdir).success();
        documents.push(fs::read(output(&workdir))?);
    }
    assert_eq!(documents[0], documents[1]);

    let workdir = directory.path().join("one");
    run(&[&forward, &reverse], &reference, &config, &workdir)
        .failure()
        .stderr(predicates::str::contains("target already exists"));
    assert_eq!(fs::read(output(&workdir))?, documents[0]);
    Ok(())
}

/// A run signal of the sample's run is recorded in provenance; with fewer
/// reads than the minimum it gives no run level, so the floor decides
/// (SRS-API-014, SRS-CALL-020).
#[test]
fn a_run_signal_of_the_samples_run_is_recorded() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let (reference, config) = inputs(directory.path())?;
    let forward = directory.path().join("read-forward.ab1");
    let reverse = directory.path().join("read-reverse.ab1");
    let mate = directory.path().join("read-mate.ab1");
    write_abif_with_secondary_signal(&forward, QUERY, 10, b'T', 400)?;
    let call = QUERY.len() - 1 - 10;
    write_abif_with_secondary_signal(&reverse, &reverse_complement(QUERY), call, b'A', 400)?;
    write_abif(&mate, QUERY)?;
    let run_signal = signal(
        "RUN1",
        &[&forward, &reverse, &mate],
        &config,
        directory.path(),
    );

    run_with_signal(
        &[&forward, &reverse],
        &reference,
        &config,
        directory.path(),
        &run_signal,
    )
    .success();

    let value = document(directory.path())?;
    assert_consistent(&value)?;
    let recorded = &value["provenance"]["signal"];
    assert_eq!(recorded["run_id"], "RUN1");
    assert_eq!(recorded["selection"], 2);
    assert_eq!(recorded["sha256"].as_str().map(str::len), Some(64));
    let first = &items(&value)?[0];
    assert_eq!(label_names(first), ["possible_heteroplasmy"]);
    let reads = first["labels"][0]["reads"]
        .as_array()
        .ok_or("reads must be an array")?;
    assert!(
        reads.iter().all(|read| {
            read["threshold"].as_f64() == Some(0.21) && read.get("calls").is_none()
        })
    );
    Ok(())
}

/// A run signal that holds none of the sample's traces is another run's:
/// the command fails without output.
#[test]
fn refuses_another_runs_signal() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let (reference, config) = inputs(directory.path())?;
    let forward = directory.path().join("read-forward.ab1");
    let reverse = directory.path().join("read-reverse.ab1");
    let other = directory.path().join("other.ab1");
    write_abif(&forward, &alternate())?;
    write_abif(&reverse, &reverse_complement(&alternate()))?;
    write_abif(&other, QUERY)?;
    let run_signal = signal("RUN2", &[&other], &config, directory.path());

    run_with_signal(
        &[&forward, &reverse],
        &reference,
        &config,
        directory.path(),
        &run_signal,
    )
    .failure();

    assert!(!output(directory.path()).exists());
    Ok(())
}
