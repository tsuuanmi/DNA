//! CLI contract for multi-read sample evidence aggregation.

pub mod support;

use std::fs;
use std::path::{Path, PathBuf};

use assert_cmd::Command;
use predicates::prelude::*;
use serde_json::Value;
use sha2::{Digest, Sha256};
use tempfile::tempdir;

use support::{
    human_mtdna_profile, write_abif, write_abif_with_incoherent_doubles,
    write_abif_with_secondary_signal, write_config, write_config_with_profile, write_reference,
};

const QUERY: &str = "ACGTCAGTACGATCGTACCTGAGTACGA";
/// A read unrelated to `TTTT{QUERY}CCCC`, so no placement reaches the identity floor.
const UNPLACEABLE: &str = "AGCATGAGTCCATGCTAGCATGACTGCA";
const SAMPLE_ID: &str = "sample-1";

fn dna_binary() -> String {
    std::env::var("CARGO_BIN_EXE_dna")
        .unwrap_or_else(|error| panic!("Cargo must provide the DNA test binary: {error}"))
}

#[test]
fn writes_deterministic_compact_sample_evidence_v9() -> Result<(), Box<dyn std::error::Error>> {
    let first = tempdir()?;
    let second = tempdir()?;

    for (directory, reverse_order) in [(first.path(), false), (second.path(), true)] {
        let reference = directory.join("reference.fa");
        let config = directory.join("dna.toml");
        let forward = directory.join("read-forward.ab1");
        let reverse = directory.join("read-reverse.ab1");

        let mut alternate = QUERY.as_bytes().to_vec();
        alternate[10] = b'A';
        let alternate = String::from_utf8(alternate)?;
        let reverse_alternate = reverse_complement(&alternate);

        write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
        write_config(&config, "linear")?;
        write_abif(&forward, QUERY)?;
        write_abif(&reverse, &reverse_alternate)?;

        let traces = if reverse_order {
            [&reverse, &forward]
        } else {
            [&forward, &reverse]
        };
        run(&traces, &reference, &config, directory)
            .success()
            .stdout(predicate::str::is_empty())
            .stderr(predicate::str::is_empty());

        let log_dir = directory.join("logs");
        let sample_log = log_dir.join(format!("{SAMPLE_ID}.log"));
        let log = fs::read_to_string(&sample_log)?;
        assert!(log.contains("event=sample_read_started"));
        assert!(log.contains("event=basecalling_completed"));
        assert!(log.contains("event=alignment_completed"));
        assert!(log.contains("event=variant_calling_completed"));
        assert!(log.contains("event=sample_read_completed"));
        assert!(!log_dir.join("read-forward.log").exists());
        assert!(!log_dir.join("read-reverse.log").exists());
    }

    let first_bytes = fs::read(sample_output_path(first.path()))?;
    let second_bytes = fs::read(sample_output_path(second.path()))?;
    assert_eq!(first_bytes, second_bytes);

    let value: Value = serde_json::from_slice(&first_bytes)?;
    assert_eq!(value["schema_version"], "dna.sample_evidence/v10");
    assert_eq!(value["sample_id"], SAMPLE_ID);
    assert_object_keys(
        &value,
        &[
            "schema_version",
            "sample_id",
            "provenance",
            "reads",
            "rejected_reads",
            "coverage",
            "overlaps",
            "locus_differences",
            "variants",
        ],
    )?;
    assert!(value.get("loci").is_none());
    assert_eq!(value["rejected_reads"], serde_json::json!([]));
    assert_eq!(
        value["provenance"]["plugins"],
        serde_json::json!([
            {"id": "sanger", "family": "modality", "version": 1},
            {"id": "core", "family": "core", "version": 1},
        ])
    );

    let reads = value["reads"].as_array().ok_or("reads must be an array")?;
    assert_eq!(reads.len(), 2);
    let forward_read = read_by_name(reads, "read-forward")?;
    let reverse_read = read_by_name(reads, "read-reverse")?;
    assert_eq!(forward_read["alignment"]["orientation"], "forward");
    assert_eq!(reverse_read["alignment"]["orientation"], "reverse");
    assert_eq!(forward_read["integrity"]["ploc_count"], QUERY.len());
    assert_eq!(reverse_read["integrity"]["ploc_count"], QUERY.len());
    assert_eq!(forward_read["integrity"]["clipped_channel_samples"], 0);
    assert_eq!(reverse_read["integrity"]["clipped_channel_samples"], 0);
    for read in [forward_read, reverse_read] {
        assert_object_keys(
            read,
            &["name", "sha256", "integrity", "callability", "alignment"],
        )?;
        assert_eq!(read["callability"]["masked_calls"], 0);
        assert_eq!(read["callability"]["callable_span"]["end"], QUERY.len());
    }
    assert!(log_text(first.path())?.contains("masked_calls_total=0 callable_calls_total=56"));

    let coverage = value["coverage"]
        .as_array()
        .ok_or("coverage must be an array")?;
    assert_eq!(coverage.len(), 1);
    assert_eq!(coverage[0]["reference"]["start"], 4);
    assert_eq!(coverage[0]["reference"]["end"], 4 + QUERY.len());
    assert_eq!(coverage[0]["read_depth"], 2);
    assert_eq!(coverage[0]["forward_depth"], 1);
    assert_eq!(coverage[0]["reverse_depth"], 1);

    let overlaps = value["overlaps"]
        .as_array()
        .ok_or("overlaps must be an array")?;
    assert_eq!(overlaps.len(), 1);
    let overlap = &overlaps[0];
    let pair = [
        overlap["left"]
            .as_str()
            .ok_or("overlap left must be a string")?,
        overlap["right"]
            .as_str()
            .ok_or("overlap right must be a string")?,
    ];
    assert!(pair.contains(&"read-forward"));
    assert!(pair.contains(&"read-reverse"));
    assert_eq!(overlap["shared_positions"], QUERY.len());
    assert_eq!(overlap["comparable_bases"], QUERY.len());
    assert_eq!(overlap["agreements"], QUERY.len() - 1);
    assert_eq!(overlap["conflicts"], 1);
    let agreement = overlap["agreement"]
        .as_f64()
        .ok_or("overlap agreement must be numeric")?;
    let expected_agreement = (QUERY.len() - 1) as f64 / QUERY.len() as f64;
    assert!((agreement - expected_agreement).abs() < 1e-12);
    assert_eq!(overlap["eligible"], true);
    assert_eq!(overlap["exclusion_reasons"], serde_json::json!([]));

    let differences = value["locus_differences"]
        .as_array()
        .ok_or("locus_differences must be an array")?;
    assert_eq!(differences.len(), 1);
    let difference = &differences[0];
    assert_eq!(difference["position"], 15);
    assert_eq!(difference["reference"], "G");
    let topology = &difference["support_topology"];
    assert_eq!(topology["reads"], 2);
    assert_eq!(topology["forward_reads"], 1);
    assert_eq!(topology["reverse_reads"], 1);
    assert_eq!(topology["reference_reads"], 1);
    assert_eq!(topology["alternate_reads"], 1);
    assert_eq!(topology["unresolved_reads"], 0);
    assert_eq!(topology["deletion_reads"], 0);
    assert_eq!(topology["profile_reads"], 2);
    assert_eq!(topology["profile_forward_reads"], 1);
    assert_eq!(topology["profile_reverse_reads"], 1);
    let observations = difference["observations"]
        .as_array()
        .ok_or("observations must be an array")?;
    assert_eq!(observations.len(), 2);
    let forward_observation = observation_for_read(observations, "read-forward")?;
    let reverse_observation = observation_for_read(observations, "read-reverse")?;
    assert_eq!(forward_observation["state"], "reference");
    assert_eq!(forward_observation["base"], "G");
    assert_eq!(forward_observation["in_noisy_region"], false);
    assert_profile(forward_observation)?;
    assert!(
        forward_observation["profile"]["G"]
            .as_f64()
            .ok_or("forward G profile must be numeric")?
            > 0.9
    );
    assert_eq!(reverse_observation["state"], "alternate");
    assert_eq!(reverse_observation["base"], "A");
    assert_eq!(reverse_observation["in_noisy_region"], false);
    assert_profile(reverse_observation)?;
    assert!(
        reverse_observation["profile"]["A"]
            .as_f64()
            .ok_or("reverse A profile must be numeric")?
            > 0.9
    );

    let variants = value["variants"]
        .as_array()
        .ok_or("variants must be an array")?;
    assert_eq!(variants.len(), 1);
    let variant = &variants[0];
    assert_eq!(variant["position"], 15);
    assert_eq!(variant["reference"], "G");
    assert_eq!(variant["alternate"], "A");
    assert_eq!(variant["kind"], "SNV");
    assert_eq!(variant["support_topology"]["reads"], 1);
    assert_eq!(variant["support_topology"]["eligible_reads"], 1);
    assert_eq!(variant["support_topology"]["forward_reads"], 0);
    assert_eq!(variant["support_topology"]["reverse_reads"], 1);
    assert_eq!(variant["support_topology"]["eligible_forward_reads"], 0);
    assert_eq!(variant["support_topology"]["eligible_reverse_reads"], 1);
    // The forward read observes the same position callably without the variant.
    assert_eq!(
        variant["opposition"],
        serde_json::json!({"reads": ["read-forward"], "forward_reads": 1, "reverse_reads": 0})
    );
    let support = variant["support"]
        .as_array()
        .ok_or("variant support must be an array")?;
    assert_eq!(support.len(), 1);
    assert_eq!(support[0]["read"], "read-reverse");
    assert_eq!(support[0]["eligible"], true);
    assert_eq!(support[0]["exclusion_reasons"], serde_json::json!([]));
    let call = &support[0]["calls"][0];
    assert_eq!(call["role"], "supporting");
    assert_eq!(call["base"], "A");
    assert_eq!(call["peaks"]["A"], 1000);
    assert!(call["quality"].is_number());
    assert!(call.get("index").is_none());
    assert!(call.get("position").is_none());
    assert!(call.get("ploc").is_none());

    let text = std::str::from_utf8(&first_bytes)?;
    for repeated in ["read_name", "read_sha256", "\"loci\""] {
        assert!(!text.contains(repeated));
    }

    Ok(())
}

#[test]
fn preserves_mixed_snv_as_ineligible_sample_evidence() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let trace = directory.path().join("mixed-read.ab1");
    let mut query = QUERY.to_owned();
    query.replace_range(10..11, "T");

    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;
    write_abif_with_secondary_signal(&trace, &query, 10, b'C', 400)?;

    let mut command = Command::new(dna_binary());
    command
        .current_dir(directory.path())
        .env("DNA_CONFIG", &config)
        .arg("sample")
        .arg(SAMPLE_ID)
        .arg(&trace)
        .arg("--reference")
        .arg(&reference)
        .assert()
        .success();

    let value: Value = serde_json::from_slice(&fs::read(sample_output_path(directory.path()))?)?;
    assert_eq!(value["schema_version"], "dna.sample_evidence/v10");
    assert_eq!(value["overlaps"], serde_json::json!([]));
    let coverage = value["coverage"]
        .as_array()
        .ok_or("coverage must be an array")?;
    assert_eq!(coverage.len(), 1);
    assert_eq!(coverage[0]["read_depth"], 1);
    let variants = value["variants"]
        .as_array()
        .ok_or("variants must be an array")?;
    assert_eq!(variants.len(), 1);
    assert_eq!(variants[0]["support_topology"]["reads"], 1);
    assert_eq!(variants[0]["support_topology"]["eligible_reads"], 0);
    assert_eq!(variants[0]["support_topology"]["forward_reads"], 1);
    assert_eq!(variants[0]["support_topology"]["reverse_reads"], 0);
    assert_eq!(variants[0]["support_topology"]["eligible_forward_reads"], 0);
    assert_eq!(variants[0]["support_topology"]["eligible_reverse_reads"], 0);
    let support = variants[0]["support"]
        .as_array()
        .ok_or("support must be an array")?;
    assert_eq!(support.len(), 1);
    assert_eq!(support[0]["read"], "mixed-read");
    assert_eq!(support[0]["eligible"], false);
    let reasons = support[0]["exclusion_reasons"]
        .as_array()
        .ok_or("exclusion_reasons must be an array")?;
    assert!(
        reasons
            .iter()
            .any(|reason| reason == "mixed_supporting_dna")
    );
    let call = &support[0]["calls"][0];
    assert_eq!(call["base"], "T");
    assert_eq!(call["peaks"]["T"], 1000);
    assert_eq!(call["peaks"]["C"], 400);
    Ok(())
}

/// Under the human-mtDNA profile, each read's calls are right-aligned and given
/// the HVS-II representation, then published as per-base notation with
/// supporting reads, and the profile identity is recorded.
#[test]
fn publishes_mtdna_notation_against_the_rcrs() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let config = directory.path().join("dna.toml");
    let trace = directory.path().join("hv2-read.ab1");
    let rcrs =
        fs::read_to_string(Path::new(env!("CARGO_MANIFEST_DIR")).join("references/rCRS.fasta"))?
            .lines()
            .filter(|line| !line.starts_with('>'))
            .collect::<String>();
    let read = format!("{}C{}", &rcrs[270..304], &rcrs[304..350]);
    write_config_with_profile(&config, &human_mtdna_profile())?;
    write_abif(&trace, &read)?;

    let mut command = Command::new(dna_binary());
    command
        .current_dir(directory.path())
        .env("DNA_CONFIG", &config)
        .arg("sample")
        .arg(SAMPLE_ID)
        .arg(&trace)
        .arg("--reference")
        .arg(Path::new(env!("CARGO_MANIFEST_DIR")).join("references/rCRS.fasta"))
        .assert()
        .success();

    let value: Value = serde_json::from_slice(&fs::read(sample_output_path(directory.path()))?)?;
    assert_eq!(value["schema_version"], "dna.sample_evidence/v10");
    assert_eq!(
        value["notation"],
        serde_json::json!({
            "style": "per_base_decimal",
            "calls": [{"call": "309.1C", "reads": ["hv2-read"]}],
        })
    );
    assert_eq!(value["provenance"]["profile"]["id"], "human-mtdna-rcrs");
    assert_eq!(
        value["provenance"]["plugins"],
        serde_json::json!([
            {"id": "sanger", "family": "modality", "version": 1},
            {"id": "core", "family": "core", "version": 1},
            {"id": "normalization", "family": "post_calling", "version": 1},
            {"id": "nomenclature", "family": "post_calling", "version": 1},
        ])
    );
    assert_eq!(
        value["provenance"]["profile"]["sha256"],
        format!("{:x}", Sha256::digest(fs::read(human_mtdna_profile())?))
    );
    Ok(())
}

/// Regression: a sample read whose indel flank is an unresolved tied call still
/// yields sample evidence; the unresolved flank is omitted from public calls.
#[test]
fn omits_unresolved_indel_flank_from_sample_evidence() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let trace = directory.path().join("tied-flank.ab1");
    write_reference(
        &reference,
        &format!("TTTT{}A{}CCCC", &QUERY[..14], &QUERY[14..]),
    )?;
    write_config(&config, "linear")?;
    write_abif_with_secondary_signal(&trace, QUERY, 13, b'A', 1000)?;

    let mut command = Command::new(dna_binary());
    command
        .current_dir(directory.path())
        .env("DNA_CONFIG", &config)
        .arg("sample")
        .arg(SAMPLE_ID)
        .arg(&trace)
        .arg("--reference")
        .arg(&reference)
        .assert()
        .success();

    let value: Value = serde_json::from_slice(&fs::read(sample_output_path(directory.path()))?)?;
    let variants = value["variants"]
        .as_array()
        .ok_or("variants must be an array")?;
    assert_eq!(variants.len(), 1);
    assert_eq!(variants[0]["kind"], "DEL");
    let calls = variants[0]["support"][0]["calls"]
        .as_array()
        .ok_or("support calls must be an array")?;
    assert_eq!(calls.len(), 1, "the unresolved flank is omitted");
    assert_eq!(calls[0]["role"], "flanking");
    assert!(matches!(
        calls[0]["base"].as_str(),
        Some("A" | "C" | "G" | "T")
    ));
    Ok(())
}

/// Writes the forward reference read and a reverse read carrying `15A` against
/// `TTTT{QUERY}CCCC`, returning the two trace paths.
/// Read callability (ADR-0067): a read with too few callable calls is recorded
/// as rejected and the remaining reads are still aggregated.
#[test]
fn records_a_read_with_too_few_callable_calls_as_rejected() -> Result<(), Box<dyn std::error::Error>>
{
    let directory = tempdir()?;
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let forward = directory.path().join("read-forward.ab1");
    let reverse = directory.path().join("read-reverse.ab1");
    let mixed = directory.path().join("read-mixed.ab1");
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;
    write_abif(&forward, QUERY)?;
    write_abif(&reverse, &reverse_complement(QUERY))?;
    write_abif_with_incoherent_doubles(&mixed, QUERY, 0..QUERY.len(), 500)?;

    run(
        &[&forward, &mixed, &reverse],
        &reference,
        &config,
        directory.path(),
    )
    .success();

    let value: Value = serde_json::from_slice(&fs::read(sample_output_path(directory.path()))?)?;
    let reads = value["reads"].as_array().ok_or("reads must be an array")?;
    assert_eq!(reads.len(), 2);
    assert!(reads.iter().all(|read| read["name"] != "read-mixed"));
    let rejected = value["rejected_reads"]
        .as_array()
        .ok_or("rejected_reads must be an array")?;
    assert_eq!(rejected.len(), 1);
    assert_object_keys(
        &rejected[0],
        &["name", "sha256", "integrity", "callability", "rejection"],
    )?;
    assert_eq!(rejected[0]["name"], "read-mixed");
    assert_eq!(
        rejected[0]["rejection"],
        serde_json::json!({
            "reason": "callable_calls_below_minimum",
            "callable_calls": 0,
            "minimum_callable_calls": 20
        })
    );
    assert_eq!(rejected[0]["callability"]["masked_calls"], QUERY.len());
    assert_eq!(value["coverage"][0]["read_depth"], 2);
    let log = log_text(directory.path())?;
    assert!(log.contains("event=sample_read_rejected read_index=1"));
    assert!(log.contains(
        "reason=callable_calls_below_minimum callable_calls=0 minimum_callable_calls=20"
    ));
    assert!(log.contains("rejected_reads=1"));
    Ok(())
}

/// A read the core cannot place is recorded as rejected and the remaining
/// reads are still aggregated (SRS-SAMPLE-029).
#[test]
fn records_an_unplaceable_read_as_rejected() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let [forward, reverse] = write_two_reads(directory.path(), &reference)?;
    let foreign = directory.path().join("read-foreign.ab1");
    write_config(&config, "linear")?;
    write_abif(&foreign, UNPLACEABLE)?;

    run(
        &[&forward, &foreign, &reverse],
        &reference,
        &config,
        directory.path(),
    )
    .success();

    let value: Value = serde_json::from_slice(&fs::read(sample_output_path(directory.path()))?)?;
    assert_eq!(value["reads"].as_array().map(Vec::len), Some(2));
    let rejected = value["rejected_reads"]
        .as_array()
        .ok_or("rejected_reads must be an array")?;
    assert_eq!(rejected.len(), 1);
    assert_eq!(rejected[0]["name"], "read-foreign");
    assert_eq!(rejected[0]["rejection"]["reason"], "identity_below_minimum");
    assert!(
        rejected[0]["rejection"]["identity"]
            .as_f64()
            .is_some_and(|identity| identity < 0.8)
    );
    assert_eq!(rejected[0]["rejection"]["minimum_identity"], 0.8);
    let log = log_text(directory.path())?;
    assert!(log.contains("event=sample_read_rejected read_index=1"));
    assert!(log.contains("reason=identity_below_minimum"));
    Ok(())
}

/// A sample whose every read is rejected fails typed and publishes nothing.
#[test]
fn fails_when_every_read_is_rejected() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let first = directory.path().join("read-first.ab1");
    let second = directory.path().join("read-second.ab1");
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;
    write_abif_with_incoherent_doubles(&first, QUERY, 0..QUERY.len(), 500)?;
    write_abif_with_incoherent_doubles(&second, QUERY, 0..QUERY.len(), 400)?;

    run(&[&first, &second], &reference, &config, directory.path())
        .failure()
        .stderr(predicate::str::contains(
            "all 2 reads were rejected and none could be analyzed",
        ));
    assert!(!sample_output_path(directory.path()).exists());
    Ok(())
}

fn write_two_reads(
    directory: &Path,
    reference: &Path,
) -> Result<[PathBuf; 2], Box<dyn std::error::Error>> {
    let mut alternate = QUERY.as_bytes().to_vec();
    alternate[10] = b'A';
    let forward = directory.join("read-forward.ab1");
    let reverse = directory.join("read-reverse.ab1");
    write_reference(reference, &format!("TTTT{QUERY}CCCC"))?;
    write_abif(&forward, QUERY)?;
    write_abif(
        &reverse,
        &reverse_complement(&String::from_utf8(alternate)?),
    )?;
    Ok([forward, reverse])
}

/// A profile that declares notation drives the generic notation path for any
/// target, not only human mtDNA.
#[test]
fn publishes_notation_for_a_non_mtdna_profile() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let [forward, reverse] = write_two_reads(directory.path(), &reference)?;
    write_config(&config, "linear")?;
    let profile = config.with_extension("profile.toml");
    let mut text = fs::read_to_string(&profile)?;
    text.push_str(
        "[normalization]\nindel_placement='right'\n[notation]\nstyle='per_base_decimal'\n",
    );
    fs::write(&profile, text)?;

    run(&[&forward, &reverse], &reference, &config, directory.path()).success();
    let value: Value = serde_json::from_slice(&fs::read(sample_output_path(directory.path()))?)?;
    assert_eq!(value["provenance"]["profile"]["id"], "synthetic-linear");
    assert_eq!(
        value["notation"],
        serde_json::json!({
            "style": "per_base_decimal",
            "calls": [{"call": "15A", "reads": ["read-reverse"]}],
        })
    );
    Ok(())
}

#[test]
fn sample_fails_closed_when_the_profile_names_another_reference()
-> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let [forward, reverse] = write_two_reads(directory.path(), &reference)?;
    write_config_with_profile(&config, &human_mtdna_profile())?;

    run(&[&forward, &reverse], &reference, &config, directory.path())
        .failure()
        .stderr(predicate::str::contains(
            "invalid target profile: reference sequence does not match profile human-mtdna-rcrs",
        ));
    assert!(!sample_output_path(directory.path()).exists());
    Ok(())
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
        .arg("sample")
        .arg(SAMPLE_ID);
    for trace in traces {
        command.arg(trace);
    }
    command.arg("--reference").arg(reference).assert()
}

fn log_text(workdir: &Path) -> Result<String, Box<dyn std::error::Error>> {
    Ok(fs::read_to_string(
        workdir.join("logs").join(format!("{SAMPLE_ID}.log")),
    )?)
}

fn sample_output_path(workdir: &Path) -> PathBuf {
    workdir
        .join("results")
        .join(format!("{SAMPLE_ID}.sample.json"))
}

fn read_by_name<'a>(
    reads: &'a [Value],
    name: &str,
) -> Result<&'a Value, Box<dyn std::error::Error>> {
    reads
        .iter()
        .find(|read| read["name"] == name)
        .ok_or_else(|| format!("missing read {name}").into())
}

fn observation_for_read<'a>(
    observations: &'a [Value],
    read: &str,
) -> Result<&'a Value, Box<dyn std::error::Error>> {
    observations
        .iter()
        .find(|observation| observation["read"] == read)
        .ok_or_else(|| format!("missing observation for read {read}").into())
}

fn assert_profile(observation: &Value) -> Result<(), Box<dyn std::error::Error>> {
    let profile = observation["profile"]
        .as_object()
        .ok_or("profile must be an object")?;
    assert_eq!(
        profile
            .keys()
            .map(String::as_str)
            .collect::<std::collections::BTreeSet<_>>(),
        ["A", "C", "G", "T"].into_iter().collect()
    );
    let total = ["A", "C", "G", "T"]
        .into_iter()
        .map(|base| {
            profile[base]
                .as_f64()
                .ok_or("profile weight must be numeric")
        })
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .sum::<f64>();
    assert!((total - 1.0).abs() < 1e-12);
    Ok(())
}

fn reverse_complement(sequence: &str) -> String {
    sequence
        .chars()
        .rev()
        .map(|base| match base {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            other => panic!("unsupported synthetic base {other}"),
        })
        .collect()
}

fn assert_object_keys(value: &Value, expected: &[&str]) -> Result<(), Box<dyn std::error::Error>> {
    let object = value.as_object().ok_or("expected object")?;
    let keys = object
        .keys()
        .map(String::as_str)
        .collect::<std::collections::BTreeSet<_>>();
    let expected = expected.iter().copied().collect();
    assert_eq!(keys, expected);
    Ok(())
}
