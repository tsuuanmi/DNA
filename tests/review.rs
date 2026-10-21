//! CLI contract for the sample review list (PROP-0004, SRS-API-014).

pub mod support;

use std::fs;
use std::path::{Path, PathBuf};

use assert_cmd::Command;
use serde_json::Value;
use tempfile::tempdir;

use support::{
    write_abif, write_abif_with_incoherent_doubles, write_abif_with_peak_heights,
    write_abif_with_secondary_signal, write_config, write_reference,
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

const RANKS: [(&str, u64); 25] = [
    ("strands_disagree", 0),
    ("read_disagrees", 0),
    ("same_strand_dissent", 0),
    ("possible_heteroplasmy", 1),
    ("reference_second", 1),
    ("estimated_length", 2),
    ("length_unconfirmed", 2),
    ("undecided", 2),
    ("not_called", 2),
    ("long_run", 2),
    ("edge_deletion", 2),
    ("run_deletion", 2),
    ("run_indel", 2),
    ("one_strand_insertion", 2),
    ("beside_not_analysed", 2),
    ("one_read", 3),
    ("one_strand", 3),
    ("noisy_signal", 4),
    ("not_covered", 4),
    ("phase_shift", 4),
    ("dephased_end", 4),
    ("read_edge", 4),
    ("poor_signal", 4),
    ("length_not_recorded", 4),
    ("tract_length", 4),
];

/// The document's structural invariants: items ordered by priority then
/// position, kinds that follow the calls, and a summary that counts them.
fn assert_consistent(value: &Value) -> Result<(), Box<dyn std::error::Error>> {
    assert_eq!(value["schema_version"], "dna.review/v1");
    assert_eq!(value["provenance"]["rules"], 32);
    // An admitted read's masked stretches are ordered, disjoint, and well formed.
    for read in value["reads"].as_array().ok_or("reads must be an array")? {
        if read["admitted"] != true {
            assert!(read.get("masked").is_none());
            assert!(read.get("mixtures").is_none());
            assert!(read.get("unreliable").is_none());
            continue;
        }
        // Each phase mixture leads with its dominant sequence, offset 0.
        for mixture in read["mixtures"]
            .as_array()
            .ok_or("mixtures must be an array")?
        {
            let components = mixture["components"].as_array().ok_or("components")?;
            assert_eq!(
                components.first().map(|first| &first["offset"]),
                Some(&Value::from(0))
            );
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
    // Zones tile the reportable positions of a called sample in order, and
    // each counts the check items touching it.
    if value["status"] == "called" {
        let zones = value["zones"].as_array().ok_or("zones must be an array")?;
        let (mut end, mut covered, mut flagged) = (0, 0, 0);
        for zone in zones {
            let first = zone["first"].as_u64().ok_or("zone first")?;
            let last = zone["last"].as_u64().ok_or("zone last")?;
            assert!(end < first && first <= last, "{zone}");
            assert_eq!(
                zone["kind"] == "repeat",
                zone.get("name").is_some(),
                "{zone}"
            );
            let touching = items
                .iter()
                .filter(|item| {
                    item["kind"] == "check"
                        && item["first"]["position"].as_u64() <= Some(last)
                        && Some(first) <= item["last"]["position"].as_u64()
                })
                .count();
            assert_eq!(
                zone["checks"].as_u64(),
                u64::try_from(touching).ok(),
                "{zone}"
            );
            flagged += u64::from(touching > 0);
            covered += last - first + 1;
            end = last;
        }
        let summary = &value["summary"];
        assert_eq!(summary["positions"]["reportable"].as_u64(), Some(covered));
        assert_eq!(
            summary["zones"]["total"].as_u64(),
            u64::try_from(zones.len()).ok()
        );
        assert_eq!(summary["zones"]["check"].as_u64(), Some(flagged));
        // Reference-only positions are increasing, inside the zones, and counted.
        let confirmed = value["reference_only"]
            .as_array()
            .ok_or("reference_only must be an array")?;
        let mut previous = 0;
        for position in confirmed {
            let position = position.as_u64().ok_or("reference_only position")?;
            assert!(previous < position);
            assert!(zones.iter().any(|zone| {
                zone["first"].as_u64() <= Some(position) && Some(position) <= zone["last"].as_u64()
            }));
            previous = position;
        }
        assert_eq!(
            summary["positions"]["reference_only"].as_u64(),
            u64::try_from(confirmed.len()).ok()
        );
    } else {
        assert!(value.get("zones").is_none());
        assert!(value.get("reference_only").is_none());
    }
    // A second peak's cause is one of the modality's labels.
    for item in items {
        for read in item["evidence"].as_array().ok_or("evidence")? {
            for position in read["positions"].as_array().ok_or("positions")? {
                if let Some(cause) = position["second"].get("cause") {
                    assert!(
                        [
                            "noise",
                            "neighbour_tail",
                            "off_centre",
                            "dye_blob",
                            "phase_shift",
                            "context",
                            "second_template",
                            "mixed_stretch"
                        ]
                        .contains(&cause.as_str().ok_or("cause")?),
                        "{cause}"
                    );
                }
            }
        }
    }
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
        // An item without calls is a check; one with calls is a variant
        // unless a label leaves it unassured.
        let calls = item["calls"].as_array().ok_or("calls must be an array")?;
        let kind = item["kind"].as_str().ok_or("kind")?;
        assert!(
            kind == "check" || (kind == "variant" && !calls.is_empty()),
            "{item}"
        );
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
// SRS-IN-010 SRS-OUT-001 SRS-API-014
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
            {"name": "read-forward", "sha256": value["reads"][0]["sha256"], "strand": "forward", "admitted": true, "masked": [], "unreliable": [], "mixtures": []},
            {"name": "read-reverse", "sha256": value["reads"][1]["sha256"], "strand": "reverse", "admitted": true, "masked": [], "unreliable": [], "mixtures": []},
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
// SRS-API-014
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

/// A trace whose calls are `sequence`, with every peak `height` high, so that
/// reads with the same calls have different content.
fn write_read(path: &Path, sequence: &str, height: i16) -> Result<(), Box<dyn std::error::Error>> {
    write_abif_with_peak_heights(path, sequence, vec![height; sequence.len()])
}

/// Two forward reads against a third forward read: the majority decides,
/// and the clean dissent of the supporting strand is a
/// `same_strand_dissent` check that lists every read (ADR-0090).
// SRS-API-014
#[test]
fn a_clean_read_against_its_own_strand_is_a_check() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let (reference, config) = inputs(directory.path())?;
    let first = directory.path().join("read-a.ab1");
    let second = directory.path().join("read-b.ab1");
    let third = directory.path().join("read-c.ab1");
    write_read(&first, &alternate(), 1000)?;
    write_read(&second, &alternate(), 1001)?;
    write_read(&third, QUERY, 1002)?;

    run(
        &[&first, &second, &third],
        &reference,
        &config,
        directory.path(),
    )
    .success();

    let value = document(directory.path())?;
    assert_consistent(&value)?;
    let first = &items(&value)?[0];
    assert_eq!(first["kind"], "check");
    assert_eq!(first["sequence"], "AA");
    assert_eq!(label_names(first)[0], "same_strand_dissent");
    let mut sides = first["labels"][0]["reads"]
        .as_array()
        .ok_or("reads must be an array")?
        .iter()
        .map(|side| {
            (
                side["read"].as_str().unwrap_or_default(),
                side["strand"].as_str().unwrap_or_default(),
                side["role"].as_str().unwrap_or_default(),
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(sides.pop(), Some(("read-c", "forward", "opposing")));
    sides.sort_unstable();
    assert_eq!(
        sides,
        [
            ("read-a", "forward", "supporting"),
            ("read-b", "forward", "supporting")
        ]
    );
    assert!(!label_names(first).contains(&"strands_disagree"));
    assert_eq!(value["summary"]["labels"]["same_strand_dissent"], 1);
    assert_eq!(value["summary"]["clear"], false);
    Ok(())
}

/// A synthetic query that holds a reference run of five `A` at reference
/// positions 17 to 21.
const RUN_QUERY: &str = "ACGTCAGTACGTAAAAAGTACCTGAGTACGA";

/// Writes `TTTT{RUN_QUERY}CCCC` over the synthetic reference.
fn write_run_reference(reference: &Path) -> Result<(), Box<dyn std::error::Error>> {
    write_reference(reference, &format!("TTTT{RUN_QUERY}CCCC"))
}

/// Runs `dna review` over one trace per entry of `strands` that reads
/// `sequence` on that strand (`true` forward), each with its own peak
/// height, over `TTTT{RUN_QUERY}CCCC`, and returns its document.
fn reads_on_strands(
    directory: &Path,
    sequence: &str,
    strands: &[bool],
) -> Result<Value, Box<dyn std::error::Error>> {
    let (reference, config) = inputs(directory)?;
    write_run_reference(&reference)?;
    let mut traces = Vec::new();
    for (index, &forward) in strands.iter().enumerate() {
        let trace = directory.join(format!("read-{index}.ab1"));
        let read = if forward {
            sequence.to_owned()
        } else {
            reverse_complement(sequence)
        };
        write_read(&trace, &read, 1000 + i16::try_from(index)?)?;
        traces.push(trace);
    }
    run(
        &traces.iter().collect::<Vec<_>>(),
        &reference,
        &config,
        directory,
    )
    .success();
    let value = document(directory)?;
    assert_consistent(&value)?;
    Ok(value)
}

/// The item of `value` that carries `label`, and that label.
fn labelled<'a>(
    value: &'a Value,
    label: &str,
) -> Result<(&'a Value, &'a Value), Box<dyn std::error::Error>> {
    let item = items(value)?
        .iter()
        .find(|item| label_names(item).contains(&label))
        .ok_or_else(|| format!("an item must carry {label}: {value}"))?;
    let found = item["labels"]
        .as_array()
        .ok_or("labels must be an array")?
        .iter()
        .find(|found| found["label"] == label)
        .ok_or(label.to_owned())?;
    Ok((item, found))
}

/// Two or three reads of one strand that delete three bases of a reference
/// run are `run_deletion`, `run_indel`, and `edge_deletion` checks that list
/// every supporting read: reads of one strand share its slips, so their
/// number does not assure the deletion (ADR-0090). A read of the other strand
/// that shows it lifts the gates.
// SRS-API-014
#[test]
fn deletions_resting_on_one_strand_are_checks_however_many_reads()
-> Result<(), Box<dyn std::error::Error>> {
    let deleted = RUN_QUERY.replacen("AAAAA", "AA", 1);
    for strands in [&[true, true][..], &[true, true, true], &[false, false]] {
        let directory = tempdir()?;
        let value = reads_on_strands(directory.path(), &deleted, strands)?;
        let strand = if strands[0] { "forward" } else { "reverse" };
        let names = (0..strands.len())
            .map(|index| format!("read-{index}"))
            .collect::<Vec<_>>();
        for label in ["run_deletion", "run_indel", "edge_deletion"] {
            let (item, found) = labelled(&value, label)?;
            assert_eq!(item["kind"], "check", "{label}: {item}");
            assert_eq!(found["strand"], strand, "{label}: {found}");
            let mut reads = serde_json::from_value::<Vec<String>>(found["reads"].clone())?;
            reads.sort_unstable();
            assert_eq!(reads, names, "{label}: {found}");
            assert!(label_names(item).contains(&"one_strand"));
            assert!(!label_names(item).contains(&"one_read"));
        }
        assert_eq!(value["summary"]["clear"], false);
    }
    let directory = tempdir()?;
    let value = reads_on_strands(directory.path(), &deleted, &[true, true, false])?;
    for label in ["run_deletion", "run_indel", "edge_deletion", "one_strand"] {
        assert!(
            labelled(&value, label).is_err(),
            "{label} with both strands: {value}"
        );
    }
    Ok(())
}

/// Two or three reads of one strand that lengthen a reference run by three
/// bases are a `long_run` check: the length rests on one strand, however
/// many of its reads show it in phase (ADR-0090).
// SRS-API-014
#[test]
fn a_long_run_resting_on_one_strand_is_a_check_however_many_reads()
-> Result<(), Box<dyn std::error::Error>> {
    let longer = RUN_QUERY.replacen("AAAAA", "AAAAAAAA", 1);
    for strands in [&[true, true][..], &[true, true, true]] {
        let directory = tempdir()?;
        let value = reads_on_strands(directory.path(), &longer, strands)?;
        let (item, found) = labelled(&value, "long_run")?;
        assert_eq!(item["kind"], "check", "{item}");
        // Every read shows the run in phase: before, only one read made the
        // length a check.
        assert_eq!(
            found["runs"],
            serde_json::json!([{"base": "A", "length": 8, "reference_length": 5, "length_evidence": "in_phase"}])
        );
        let (_, indel) = labelled(&value, "run_indel")?;
        assert_eq!(indel["reads"].as_array().map(Vec::len), Some(strands.len()));
    }
    Ok(())
}

/// Two or three reads of one strand that read another base than the
/// reference while each shows the reference base as a second peak are a
/// `reference_second` check at every such read: reruns of one strand repeat
/// its dye blobs (ADR-0090). A read of the other strand that shows the base
/// lifts the gate.
// SRS-API-014
#[test]
fn a_reference_second_peak_on_one_strand_is_a_check_at_every_read()
-> Result<(), Box<dyn std::error::Error>> {
    for (count, other) in [(2, false), (3, false), (2, true)] {
        let directory = tempdir()?;
        let (reference, config) = inputs(directory.path())?;
        let mut traces = Vec::new();
        for index in 0..count {
            let trace = directory.path().join(format!("read-{index}.ab1"));
            // Reference position 15 (call 10) reads A over a G second peak.
            write_abif_with_secondary_signal(&trace, &alternate(), 10, b'G', 300 + index)?;
            traces.push(trace);
        }
        if other {
            let trace = directory.path().join("read-reverse.ab1");
            write_abif(&trace, &reverse_complement(&alternate()))?;
            traces.push(trace);
        }
        run(
            &traces.iter().collect::<Vec<_>>(),
            &reference,
            &config,
            directory.path(),
        )
        .success();
        let value = document(directory.path())?;
        assert_consistent(&value)?;
        let variant = items(&value)?
            .iter()
            .find(|item| item["calls"] == serde_json::json!(["15A"]))
            .ok_or("the substitution must be an item")?;
        let mut flagged = variant["labels"]
            .as_array()
            .ok_or("labels must be an array")?
            .iter()
            .filter(|label| label["label"] == "reference_second")
            .map(|label| {
                assert_eq!(label["position"], 15);
                assert_eq!(label["strand"], "forward");
                label["read"].as_str().unwrap_or_default()
            })
            .collect::<Vec<_>>();
        flagged.sort_unstable();
        if other {
            assert!(flagged.is_empty(), "{variant}");
            assert_eq!(variant["kind"], "variant");
        } else {
            assert_eq!(
                flagged,
                (0..count)
                    .map(|index| format!("read-{index}"))
                    .collect::<Vec<_>>(),
                "{variant}"
            );
            assert_eq!(variant["kind"], "check");
            assert!(!label_names(variant).contains(&"one_read"));
        }
    }
    Ok(())
}

/// Six traces of one file name from six folders, three per strand, are six
/// reads: each is known by its content, and its name gains the prefix of its
/// content SHA-256 that tells it from the others. Two copies of one trace are
/// refused.
// SRS-API-014
#[test]
fn traces_of_one_name_are_told_apart_by_content() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let (reference, config) = inputs(directory.path())?;
    let mut traces = Vec::new();
    for (index, folder) in ["a", "b", "c", "d", "e", "f"].into_iter().enumerate() {
        let folder = directory.path().join(folder);
        fs::create_dir(&folder)?;
        let trace = folder.join("read.ab1");
        let sequence = if index < 3 {
            alternate()
        } else {
            reverse_complement(&alternate())
        };
        write_read(&trace, &sequence, 1000 + i16::try_from(index)?)?;
        traces.push(trace);
    }

    let workdir = directory.path().join("run");
    fs::create_dir(&workdir)?;
    run(
        &traces.iter().collect::<Vec<_>>(),
        &reference,
        &config,
        &workdir,
    )
    .success();

    let value = document(&workdir)?;
    assert_consistent(&value)?;
    let reads = value["reads"].as_array().ok_or("reads must be an array")?;
    assert_eq!(reads.len(), 6);
    for read in reads {
        let sha256 = read["sha256"].as_str().ok_or("sha256")?;
        assert_eq!(read["name"], format!("read@{}", &sha256[..8]));
    }
    let strands = reads
        .iter()
        .filter(|read| read["strand"] == "forward")
        .count();
    assert_eq!(strands, 3);
    let variant = items(&value)?
        .iter()
        .find(|item| item["kind"] == "variant")
        .ok_or("the substitution must be an item")?;
    assert_eq!(variant["calls"], serde_json::json!(["15A"]));
    assert_eq!(variant["evidence"].as_array().map(Vec::len), Some(6));

    let copies = directory.path().join("copies");
    fs::create_dir(&copies)?;
    run(&[&traces[0], &traces[0]], &reference, &config, &copies)
        .failure()
        .stderr(predicates::str::contains("duplicate input trace content"));
    assert!(!output(&copies).exists());
    Ok(())
}

/// A second base on both strands at one position, beyond the reads'
/// neighbours, is a `possible_heteroplasmy` check with each read's share.
// SRS-API-014
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
    // Every call of these short reads is withheld for its relative quality,
    // so only reference confirmations (ADR-0080) decide the positions around
    // 15. Both reads show the second base there unexplained, so 15 is not
    // confirmed: it stays undecided, a check that is not callable.
    assert_eq!(label_names(first), ["possible_heteroplasmy", "undecided"]);
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
    // Undecided, the position has no sequence.
    assert!(first["sequence"].is_null());
    assert_eq!(first["evidence"][0]["positions"][0]["second"]["base"], "T");
    Ok(())
}

/// A sample whose every read is rejected still gets its document, with the
/// rejected reads and their causes and no items.
// SRS-API-014
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
// SRS-API-014 SRS-OUT-006 SRS-OUT-007
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
// SRS-API-014
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
    // As in `a_second_base_on_both_strands_is_possible_heteroplasmy`: the
    // position stays undecided between reference confirmations.
    assert_eq!(label_names(first), ["possible_heteroplasmy", "undecided"]);
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
// SRS-API-014
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
