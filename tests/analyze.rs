//! CLI contract for single-read reference analysis: JSON output, logs, and failures.

pub mod support;

use std::collections::BTreeSet;
use std::fs;
use std::path::Path;

use assert_cmd::Command;
use predicates::prelude::*;
use serde_json::Value;
use tempfile::tempdir;

use support::{
    analysis_output_path, write_abif, write_abif_with_amplitude_decay,
    write_abif_with_background_noise, write_abif_with_channel_order,
    write_abif_with_incoherent_doubles, write_abif_with_peak_heights, write_abif_with_ploc,
    write_abif_with_secondary_signal, write_abif_with_secondary_signals,
    write_abif_with_shadow_ladder, write_abif_with_short_pbas, write_abif_with_unused_p2ba,
    write_abif_with_vendor, write_config, write_reference,
};

const QUERY: &str = "ACGTCAGTACGATCGTACCTGAGTACGA";

fn dna_binary() -> String {
    std::env::var("CARGO_BIN_EXE_dna")
        .unwrap_or_else(|error| panic!("Cargo must provide the DNA test binary: {error}"))
}
const REPRESENTATIVE_TRACE_STEM: &str = "E01_20260504_2544522_HV1F_13";

#[test]
fn writes_deterministic_compact_json() -> Result<(), Box<dyn std::error::Error>> {
    let first = tempdir()?;
    let second = tempdir()?;
    for directory in [first.path(), second.path()] {
        let trace = directory.join(format!("{REPRESENTATIVE_TRACE_STEM}.ab1"));
        let reference = directory.join("reference.fa");
        let config = directory.join("dna.toml");
        write_abif(&trace, QUERY)?;
        write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
        write_config(&config, "linear")?;
        run(&trace, &reference, &config, directory)
            .success()
            .stdout(predicate::str::is_empty())
            .stderr(predicate::str::is_empty());
    }

    let first_trace = first
        .path()
        .join(format!("{REPRESENTATIVE_TRACE_STEM}.ab1"));
    let second_trace = second
        .path()
        .join(format!("{REPRESENTATIVE_TRACE_STEM}.ab1"));
    let first_output = analysis_output_path(first.path(), &first_trace);
    assert_eq!(
        first_output,
        first
            .path()
            .join("results")
            .join(format!("{REPRESENTATIVE_TRACE_STEM}.json"))
    );
    let first_bytes = fs::read(first_output)?;
    let second_bytes = fs::read(analysis_output_path(second.path(), &second_trace))?;
    assert_eq!(first_bytes, second_bytes);
    let value: Value = serde_json::from_slice(&first_bytes)?;
    assert_eq!(value["schema_version"], "dna.analysis/v9");
    assert_object_keys(
        &value,
        &[
            "schema_version",
            "provenance",
            "read",
            "signal_quality",
            "alignment",
            "variants",
            "warnings",
        ],
    );
    assert_object_keys(
        &value["provenance"],
        &[
            "input",
            "reference",
            "configuration_sha256",
            "profile",
            "plugins",
        ],
    );
    assert_eq!(
        value["provenance"]["plugins"],
        serde_json::json!([
            {"id": "sanger", "family": "modality", "version": 1},
            {"id": "core", "family": "core", "version": 1},
        ])
    );
    assert_eq!(value["provenance"]["profile"]["id"], "synthetic-linear");
    assert_object_keys(&value["provenance"]["profile"], &["id", "sha256"]);
    assert_object_keys(&value["read"], &["call_count", "trim", "callability"]);
    assert_object_keys(
        &value["read"]["callability"],
        &["callable_span", "segments", "masked_calls"],
    );
    assert_object_keys(
        &value["read"]["callability"]["segments"][0],
        &["calls", "state", "after_repeat"],
    );
    assert_eq!(value["read"]["callability"]["masked_calls"], 0);
    assert_eq!(
        value["read"]["callability"]["segments"][0]["state"],
        "in_phase"
    );
    assert_object_keys(
        &value["alignment"],
        &[
            "orientation",
            "callable_bases",
            "identity",
            "unresolved_bases",
            "masked_bases",
            "gap_opens",
            "reference_segments",
            "callable_reference_segments",
            "wraps_origin",
        ],
    );
    assert_object_keys(
        &value["warnings"],
        &[
            "unresolved_primary_calls",
            "multi_channel_unresolved_calls",
            "ploc_vendor_length_mismatches",
            "clipped_channel_samples",
            "excluded_variant_candidates",
        ],
    );
    assert_eq!(value["read"]["call_count"], 28);
    assert_eq!(value["alignment"]["orientation"], "forward");
    assert!(value["signal_quality"]["noisy_regions"].is_array());
    assert_eq!(value["signal_quality"]["integrity"]["ploc_count"], 28);
    assert_eq!(
        value["signal_quality"]["integrity"]["clipped_channel_samples"],
        0
    );
    assert!(value.get("meta").is_none());
    assert!(value.get("sequence").is_none());
    assert!(value.pointer("/signal_quality/windows").is_none());
    let text = std::str::from_utf8(&first_bytes)?;
    for obsolete in ["evidence", "position_1based", "_0based", "_exclusive"] {
        assert!(!text.contains(obsolete));
    }
    assert!(
        !first
            .path()
            .join(format!("results/{REPRESENTATIVE_TRACE_STEM}.vcf"))
            .exists()
    );
    let log = fs::read_to_string(
        first
            .path()
            .join(format!("logs/{REPRESENTATIVE_TRACE_STEM}.log")),
    )?;
    let mut search_start = 0;
    for event in [
        "event=analysis_started",
        "event=inputs_loaded",
        "event=basecalling_completed",
        "event=signal_processing_completed",
        "event=callability_completed",
        "event=quality_control_completed",
        "event=read_evidence_completed",
        "event=alignment_completed",
        "event=variant_calling_completed",
        "event=result_ready_for_publication",
    ] {
        let offset = log[search_start..]
            .find(event)
            .ok_or_else(|| format!("missing ordered log event {event}"))?;
        search_start += offset + event.len();
    }
    assert!(log.lines().all(|line| line.contains("run_id=")));
    assert!(log.contains("calls=28 canonical_primary=28 unresolved_primary=0"));
    assert!(log.contains("retained=28"));
    assert!(log.contains(
        "calls=28 repeats=0 segments=1 in_phase_segments=1 dephased_segments=0 mixed_segments=0 weak_segments=0 irregular_segments=0 masked_calls=0 callable=0..28 callable_fraction=1.0000 segment_map=0..28:in_phase window_calls=16"
    ));
    assert!(log.contains(
        "windows=19 noisy_windows=0 noisy_regions=0 noisy_calls=0 window_size_bases=10 minimum_noisy_windows=2"
    ));
    assert!(log.contains("orientation=Forward"));
    assert!(log.contains("reported=0 snv=0 insertion=0 deletion=0 excluded=0"));
    assert!(log.contains("minimum_peak_height=150"));
    assert!(!log.contains(QUERY));
    assert!(!log.contains("[[1, 50000]]"));
    assert!(!log.contains("\"schema_version\""));
    assert!(!log.contains("gapped_query"));
    Ok(())
}

#[test]
fn fails_closed_when_the_profile_names_another_reference() -> Result<(), Box<dyn std::error::Error>>
{
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    write_abif(&trace, QUERY)?;
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;
    let profile = config.with_extension("profile.toml");
    let pinned = fs::read_to_string(&profile)?.replace(
        "topology='linear'",
        &format!("topology='linear'\nsequence_sha256='{}'", "0".repeat(64)),
    );
    fs::write(&profile, pinned)?;

    run(&trace, &reference, &config, directory.path())
        .failure()
        .stderr(predicate::str::contains(
            "invalid target profile: reference sequence does not match profile synthetic-linear",
        ));
    assert!(!analysis_output_path(directory.path(), &trace).exists());
    let log = fs::read_to_string(directory.path().join("logs/trace.log"))?;
    assert!(log.contains("event=analysis_failed stage=input_loading"));
    Ok(())
}

#[test]
fn requires_the_configured_profile_file() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    write_abif(&trace, QUERY)?;
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;
    fs::remove_file(config.with_extension("profile.toml"))?;

    run(&trace, &reference, &config, directory.path())
        .failure()
        .stderr(predicate::str::contains("failed to read profile file"));
    assert!(!analysis_output_path(directory.path(), &trace).exists());
    Ok(())
}

#[test]
fn reorders_noncanonical_fwo_channels() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    write_abif_with_channel_order(&trace, QUERY, *b"TGCA")?;
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;

    run(&trace, &reference, &config, directory.path()).success();
    let value = read_result(directory.path(), &trace)?;
    assert_eq!(value["read"]["call_count"], QUERY.len());
    assert_eq!(value["alignment"]["identity"], 1.0);
    assert_eq!(value["variants"].as_array().map(Vec::len), Some(0));
    Ok(())
}

#[test]
fn ignores_unused_p2ba_content() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    write_abif_with_unused_p2ba(&trace, QUERY, vec![b'!'])?;
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;

    run(&trace, &reference, &config, directory.path()).success();
    Ok(())
}

#[test]
fn rejects_non_increasing_ploc_without_output() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let mut ploc: Vec<usize> = (0..QUERY.len()).map(|index| 2 + 4 * index).collect();
    ploc[5] = ploc[4];
    write_abif_with_ploc(&trace, QUERY, ploc)?;
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;

    run(&trace, &reference, &config, directory.path())
        .failure()
        .stderr(predicate::str::contains("strictly increasing"));
    assert!(!analysis_output_path(directory.path(), &trace).exists());
    Ok(())
}

#[test]
fn rejects_out_of_range_ploc_without_output() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let mut ploc: Vec<usize> = (0..QUERY.len()).map(|index| 2 + 4 * index).collect();
    *ploc.last_mut().ok_or("missing synthetic PLOC")? = 30_000;
    write_abif_with_ploc(&trace, QUERY, ploc)?;
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;

    run(&trace, &reference, &config, directory.path())
        .failure()
        .stderr(predicate::str::contains("outside channel samples"));
    assert!(!analysis_output_path(directory.path(), &trace).exists());
    Ok(())
}

#[test]
fn preserves_vendor_length_mismatch_as_integrity_evidence() -> Result<(), Box<dyn std::error::Error>>
{
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    write_abif_with_short_pbas(&trace, QUERY)?;
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;

    run(&trace, &reference, &config, directory.path()).success();
    let value = read_result(directory.path(), &trace)?;
    assert_eq!(value["read"]["call_count"], QUERY.len());
    assert_eq!(
        value["signal_quality"]["integrity"]["vendor_primary_count"],
        QUERY.len() - 1
    );
    assert_eq!(
        value["signal_quality"]["integrity"]["vendor_quality_count"],
        QUERY.len()
    );
    assert_eq!(value["warnings"]["ploc_vendor_length_mismatches"], 1);
    Ok(())
}

#[test]
fn reports_exact_signal_clipping_without_reclassifying_the_read()
-> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let mut heights = vec![1000_i16; QUERY.len()];
    heights[10] = i16::MAX;
    write_abif_with_peak_heights(&trace, QUERY, heights)?;
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;

    run(&trace, &reference, &config, directory.path()).success();
    let value = read_result(directory.path(), &trace)?;
    assert_eq!(
        value["signal_quality"]["integrity"]["clipped_channel_samples"],
        1
    );
    assert_eq!(value["warnings"]["clipped_channel_samples"], 1);
    assert!(
        value["signal_quality"]["integrity"]["maximum_to_median_event_signal_ratio"]
            .as_f64()
            .is_some_and(|ratio| ratio > 30.0)
    );
    assert_eq!(value["alignment"]["orientation"], "forward");
    Ok(())
}

#[test]
fn processes_only_valid_ploc_loci_when_vendor_series_are_longer()
-> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let mut ploc: Vec<usize> = (0..QUERY.len()).map(|index| 2 + 4 * index).collect();
    ploc.pop();
    write_abif_with_ploc(&trace, QUERY, ploc)?;
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;

    run(&trace, &reference, &config, directory.path()).success();
    let value = read_result(directory.path(), &trace)?;
    assert_eq!(value["read"]["call_count"], QUERY.len() - 1);
    assert_eq!(
        value["signal_quality"]["integrity"]["ploc_count"],
        QUERY.len() - 1
    );
    assert_eq!(
        value["signal_quality"]["integrity"]["vendor_primary_count"],
        QUERY.len()
    );
    assert_eq!(
        value["signal_quality"]["integrity"]["vendor_quality_count"],
        QUERY.len()
    );
    assert_eq!(value["warnings"]["ploc_vendor_length_mismatches"], 2);
    Ok(())
}

#[test]
fn reports_snv_with_peaks_and_quality() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    write_abif(&trace, QUERY)?;
    let mut reference_query = QUERY.as_bytes().to_vec();
    reference_query[10] = if reference_query[10] == b'A' {
        b'C'
    } else {
        b'A'
    };
    write_reference(
        &reference,
        &format!("TTTT{}CCCC", String::from_utf8(reference_query)?),
    )?;
    write_config(&config, "linear")?;

    run(&trace, &reference, &config, directory.path()).success();
    let value = read_result(directory.path(), &trace)?;
    let variants = value["variants"]
        .as_array()
        .ok_or("variants is not an array")?;
    assert_eq!(variants.len(), 1);
    let variant = &variants[0];
    assert_object_keys(
        variant,
        &["position", "reference", "alternate", "kind", "calls"],
    );
    assert_eq!(variant["kind"], "SNV");
    assert_eq!(variant["position"], 15);
    let call = &variant["calls"][0];
    assert_call_evidence(call);
    assert_eq!(call["role"], "supporting");
    assert_eq!(call["base"], variant["alternate"]);
    assert!(call["quality"].is_number());
    assert!(
        !call
            .as_object()
            .ok_or("call is not an object")?
            .contains_key("index")
    );
    assert!(
        !call
            .as_object()
            .ok_or("call is not an object")?
            .contains_key("ploc")
    );
    Ok(())
}

#[test]
fn excludes_mixed_supporting_snv_without_erasing_the_observation()
-> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let mut query = QUERY.to_owned();
    query.replace_range(10..11, "T");

    write_abif_with_secondary_signal(&trace, &query, 10, b'C', 400)?;
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;

    run(&trace, &reference, &config, directory.path()).success();
    let value = read_result(directory.path(), &trace)?;
    assert_eq!(value["variants"].as_array().map(Vec::len), Some(0));
    assert_eq!(value["warnings"]["excluded_variant_candidates"], 1);

    let log = fs::read_to_string(directory.path().join("logs/trace.log"))?;
    assert!(log.contains("event=variant_removed kind=SNV"));
    assert!(log.contains("mixed_supporting_dna"));
    Ok(())
}

#[test]
fn annotates_noisy_region_without_filtering_supported_snv() -> Result<(), Box<dyn std::error::Error>>
{
    let low_threshold = tempdir()?;
    let high_threshold = tempdir()?;
    let mut query = QUERY.to_owned();
    query.replace_range(10..11, "T");
    let mut results = Vec::new();

    for (directory, threshold) in [(&low_threshold, 0.1), (&high_threshold, 10_000.0)] {
        let trace = directory.path().join("trace.ab1");
        let reference = directory.path().join("reference.fa");
        let config = directory.path().join("dna.toml");
        write_abif_with_background_noise(&trace, &query, 7..14, 100)?;
        write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
        write_config(&config, "linear")?;
        let config_text = fs::read_to_string(&config)?;
        fs::write(
            &config,
            config_text.replace(
                "minimum_primary_snr=3.0",
                &format!("minimum_primary_snr={threshold}"),
            ),
        )?;
        run(&trace, &reference, &config, directory.path()).success();
        results.push(read_result(directory.path(), &trace)?);
    }

    assert_eq!(results[0]["read"], results[1]["read"]);
    assert_eq!(results[0]["alignment"], results[1]["alignment"]);
    assert_eq!(results[0]["variants"], results[1]["variants"]);
    assert!(
        results[0]["signal_quality"]["noisy_regions"]
            .as_array()
            .is_some_and(Vec::is_empty)
    );

    let variants = results[1]["variants"]
        .as_array()
        .ok_or("variants is not an array")?;
    assert_eq!(variants.len(), 1);
    assert_eq!(variants[0]["kind"], "SNV");
    assert_eq!(variants[0]["calls"][0]["base"], variants[0]["alternate"]);
    assert!(variants[0]["calls"][0]["quality"].is_number());
    let regions = results[1]["signal_quality"]["noisy_regions"]
        .as_array()
        .ok_or("noisy_regions is not an array")?;
    assert!(regions.iter().any(|region| {
        region["calls"]["start"]
            .as_u64()
            .is_some_and(|start| start <= 10)
            && region["calls"]["end"].as_u64().is_some_and(|end| end > 10)
    }));
    assert_eq!(results[1]["warnings"]["excluded_variant_candidates"], 0);
    Ok(())
}

#[test]
fn reports_extracted_snv_at_peak_floor() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let mut query = QUERY.to_owned();
    query.replace_range(10..11, "T");
    let mut peaks = vec![1000; query.len()];
    peaks[10] = 150;
    write_abif_with_peak_heights(&trace, &query, peaks)?;
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;

    run(&trace, &reference, &config, directory.path()).success();
    let value = read_result(directory.path(), &trace)?;
    assert_eq!(value["variants"].as_array().map(Vec::len), Some(1));
    assert_eq!(value["variants"][0]["kind"], "SNV");
    assert_eq!(value["warnings"]["excluded_variant_candidates"], 0);
    Ok(())
}

#[test]
fn filters_extracted_snv_below_peak_floor() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let mut query = QUERY.to_owned();
    query.replace_range(10..11, "T");
    let mut peaks = vec![1000; query.len()];
    peaks[10] = 149;
    write_abif_with_peak_heights(&trace, &query, peaks)?;
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;

    run(&trace, &reference, &config, directory.path()).success();
    let value = read_result(directory.path(), &trace)?;
    assert_eq!(value["variants"].as_array().map(Vec::len), Some(0));
    assert_eq!(value["warnings"]["excluded_variant_candidates"], 1);
    let log = fs::read_to_string(directory.path().join("logs/trace.log"))?;
    assert!(log.contains("event=variant_calling_completed"));
    assert!(log.contains("reported=0 snv=0 insertion=0 deletion=0 excluded=1"));
    let removed = log
        .lines()
        .find(|line| line.contains("event=variant_removed"))
        .ok_or("missing removed-variant log record")?;
    assert!(removed.contains("kind=SNV contig=\"synthetic\" position=15"));
    assert!(removed.contains("reasons=peak_below_minimum"));
    assert!(!removed.contains("reference="));
    assert!(!removed.contains("alternate="));
    assert!(log.contains(" | WARN     | "));
    assert!(log.contains("event=warning_summary"));
    assert!(log.contains("excluded_variant_candidates=1"));
    Ok(())
}

#[test]
fn filters_by_normalized_anchor_region() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let mut query = QUERY.to_owned();
    query.replace_range(10..11, "T");
    write_abif(&trace, &query)?;
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;
    let profile = config.with_extension("profile.toml");
    let restricted =
        fs::read_to_string(&profile)?.replace("regions=[[1, 50000]]", "regions=[[16, 16]]");
    fs::write(&profile, restricted)?;

    run(&trace, &reference, &config, directory.path()).success();
    let value = read_result(directory.path(), &trace)?;
    assert_eq!(value["variants"].as_array().map(Vec::len), Some(0));
    assert_eq!(value["warnings"]["excluded_variant_candidates"], 1);
    let log = fs::read_to_string(directory.path().join("logs/trace.log"))?;
    assert!(log.contains(
        "event=variant_removed kind=SNV contig=\"synthetic\" position=15 reasons=outside_target_region"
    ));
    Ok(())
}

#[test]
fn maps_reverse_snv_to_original_call_and_ploc() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let mut biological_query = QUERY.to_owned();
    biological_query.replace_range(10..11, "T");
    write_abif(&trace, &reverse_complement(&biological_query))?;
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;

    run(&trace, &reference, &config, directory.path()).success();
    let value = read_result(directory.path(), &trace)?;
    assert_eq!(value["alignment"]["orientation"], "reverse");
    let variant = &value["variants"][0];
    assert_eq!(variant["position"], 15);
    assert_eq!(variant["reference"], "G");
    assert_eq!(variant["alternate"], "T");
    let call = &variant["calls"][0];
    assert_call_evidence(call);
    assert_eq!(call["base"], "T");
    assert_eq!(call["peaks"]["T"], 1000);
    Ok(())
}

#[test]
fn reports_insertion_support_and_flanks() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let mut peaks = vec![1000; QUERY.len()];
    peaks[11] = 1;
    peaks[13] = 1;
    write_abif_with_peak_heights(&trace, QUERY, peaks)?;
    let reference_query = format!("{}{}", &QUERY[..12], &QUERY[13..]);
    write_reference(&reference, &format!("TTTT{reference_query}CCCC"))?;
    write_config(&config, "linear")?;

    run(&trace, &reference, &config, directory.path()).success();
    let value = read_result(directory.path(), &trace)?;
    let variant = &value["variants"][0];
    assert_eq!(variant["kind"], "INS");
    assert_eq!(variant["position"], 16);
    assert_eq!(variant["reference"], "A");
    assert_eq!(variant["alternate"], "AT");
    let calls = variant["calls"].as_array().ok_or("calls is not an array")?;
    assert_eq!(calls.len(), 3);
    assert_eq!(calls[0]["role"], "supporting");
    assert_eq!(calls[0]["base"], "T");
    assert_eq!(calls[1]["role"], "flanking");
    assert_eq!(calls[2]["role"], "flanking");
    for call in calls {
        assert_call_evidence(call);
    }
    Ok(())
}

#[test]
fn filters_multibase_insertion_when_any_inserted_peak_is_low()
-> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let mut peaks = vec![1000; QUERY.len()];
    peaks[13] = 149;
    write_abif_with_peak_heights(&trace, QUERY, peaks)?;
    let reference_query = format!("{}{}", &QUERY[..12], &QUERY[14..]);
    write_reference(&reference, &format!("TTTT{reference_query}CCCC"))?;
    write_config(&config, "linear")?;

    run(&trace, &reference, &config, directory.path()).success();
    let value = read_result(directory.path(), &trace)?;
    assert_eq!(value["variants"].as_array().map(Vec::len), Some(0));
    assert_eq!(value["warnings"]["excluded_variant_candidates"], 1);
    Ok(())
}

#[test]
fn reports_deletion_with_flanks_only() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let mut peaks = vec![1000; QUERY.len()];
    peaks[10] = 1;
    peaks[11] = 1;
    write_abif_with_peak_heights(&trace, QUERY, peaks)?;
    let reference_query = format!("{}A{}", &QUERY[..12], &QUERY[12..]);
    write_reference(&reference, &format!("TTTT{reference_query}CCCC"))?;
    write_config(&config, "linear")?;

    run(&trace, &reference, &config, directory.path()).success();
    let value = read_result(directory.path(), &trace)?;
    let variant = &value["variants"][0];
    assert_eq!(variant["kind"], "DEL");
    assert_eq!(variant["position"], 16);
    assert_eq!(variant["reference"], "AA");
    assert_eq!(variant["alternate"], "A");
    let calls = variant["calls"].as_array().ok_or("calls is not an array")?;
    assert_eq!(calls.len(), 2);
    assert!(calls.iter().all(|item| item["role"] == "flanking"));
    for call in calls {
        assert_call_evidence(call);
    }
    Ok(())
}

#[test]
fn represents_circular_origin_wrap() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let rotated = format!("{}{}", &QUERY[18..], &QUERY[..18]);
    write_abif(&trace, &rotated)?;
    write_reference(&reference, QUERY)?;
    write_config(&config, "circular")?;

    run(&trace, &reference, &config, directory.path()).success();
    let value = read_result(directory.path(), &trace)?;
    assert_eq!(value["alignment"]["wraps_origin"], true);
    assert_eq!(
        value["alignment"]["reference_segments"]
            .as_array()
            .map(Vec::len),
        Some(2)
    );
    assert_eq!(value["alignment"]["reference_segments"][0]["start"], 18);
    assert_eq!(value["alignment"]["reference_segments"][0]["end"], 28);
    assert_eq!(value["alignment"]["reference_segments"][1]["start"], 0);
    assert_eq!(value["alignment"]["reference_segments"][1]["end"], 18);
    assert!(value["warnings"].get("reference_origin_wrap").is_none());
    Ok(())
}

#[test]
fn maps_circular_origin_snv_to_original_call_and_ploc() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let mut rotated = format!("{}{}", &QUERY[18..], &QUERY[..18]);
    rotated.replace_range(12..13, "A");
    write_abif(&trace, &rotated)?;
    write_reference(&reference, QUERY)?;
    write_config(&config, "circular")?;

    run(&trace, &reference, &config, directory.path()).success();
    let value = read_result(directory.path(), &trace)?;
    assert_eq!(value["alignment"]["wraps_origin"], true);
    let variant = &value["variants"][0];
    assert_eq!(variant["position"], 3);
    assert_eq!(variant["reference"], "G");
    assert_eq!(variant["alternate"], "A");
    let call = &variant["calls"][0];
    assert_call_evidence(call);
    assert_eq!(call["base"], "A");
    Ok(())
}

#[test]
fn accepts_iupac_vendor_calls_and_char_pcon() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let mut vendor = QUERY.to_owned();
    vendor.replace_range(5..6, "K");
    write_abif_with_vendor(&trace, QUERY, &vendor, 2)?;
    let mut reference_query = QUERY.to_owned();
    reference_query.replace_range(5..6, if &QUERY[5..6] == "A" { "C" } else { "A" });
    write_reference(&reference, &format!("TTTT{reference_query}CCCC"))?;
    write_config(&config, "linear")?;

    run(&trace, &reference, &config, directory.path()).success();
    let value = read_result(directory.path(), &trace)?;
    let call = &value["variants"][0]["calls"][0];
    assert_call_evidence(call);
    assert!(call["quality"].is_number());
    assert!(call.get("relative_quality").is_none());
    Ok(())
}

#[test]
fn malformed_abif_leaves_no_output() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    fs::write(&trace, b"not an ABIF file")?;
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;

    run(&trace, &reference, &config, directory.path())
        .failure()
        .stderr(predicate::str::contains("invalid ABIF input"));
    assert!(!analysis_output_path(directory.path(), &trace).exists());
    let log = fs::read_to_string(directory.path().join("logs/trace.log"))?;
    assert!(log.contains(" | ERROR    | "));
    assert!(log.contains("event=analysis_failed stage=input_loading"));
    assert!(log.contains("invalid ABIF input"));
    Ok(())
}

#[test]
fn refuses_to_overwrite_completed_output() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    write_abif(&trace, QUERY)?;
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;
    let output = analysis_output_path(directory.path(), &trace);
    fs::create_dir_all(output.parent().ok_or("output has no parent")?)?;
    fs::write(&output, b"owned")?;

    run(&trace, &reference, &config, directory.path())
        .failure()
        .stderr(predicate::str::contains("target already exists"));
    assert_eq!(fs::read(output)?, b"owned");
    let log = fs::read_to_string(directory.path().join("logs/trace.log"))?;
    assert!(log.contains("event=analysis_started"));
    assert!(log.contains("event=analysis_failed stage=input_loading"));
    assert!(log.contains("target already exists"));
    assert!(log.contains(" | ERROR    | "));
    Ok(())
}

/// Read callability (ADR-0067): an SNV right after a long homopolymer is
/// reported while the signal stays in phase, and is no candidate at all once
/// a shadow ladder masks the calls after the run.
#[test]
fn reports_an_snv_after_a_homopolymer_only_while_the_signal_is_in_phase()
-> Result<(), Box<dyn std::error::Error>> {
    let read = format!(
        "ACGTAGTCAGTACG{}TAGCTAGCATGCATGACTGACTAGCATGCA",
        "C".repeat(9)
    );
    let mut reference_read = read.clone();
    reference_read.replace_range(26..27, "T");
    for (shadows, reported) in [(&[][..], 1), (&[(-1, 0.4)][..], 0)] {
        let directory = tempdir()?;
        let trace = directory.path().join("trace.ab1");
        let reference = directory.path().join("reference.fa");
        let config = directory.path().join("dna.toml");
        write_abif_with_shadow_ladder(&trace, &read, 23, shadows)?;
        write_reference(&reference, &format!("TTTT{reference_read}CCCC"))?;
        write_config(&config, "linear")?;

        run(&trace, &reference, &config, directory.path()).success();

        let value = read_result(directory.path(), &trace)?;
        assert_eq!(
            value["variants"].as_array().map(Vec::len),
            Some(reported),
            "shadows {shadows:?}"
        );
        // The dephased tail is trimmed: the trim interval is the callable span.
        let trim_end = if reported == 1 { read.len() } else { 23 };
        assert_eq!(
            value["read"]["trim"]["end"], trim_end,
            "shadows {shadows:?}"
        );
        assert_eq!(value["alignment"]["masked_bases"], 0, "shadows {shadows:?}");
        let log = fs::read_to_string(directory.path().join("logs/trace.log"))?;
        assert!(
            !log.contains("event=variant_removed"),
            "shadows {shadows:?}"
        );
    }
    Ok(())
}

/// Read callability (ADR-0067): an internal mixed-signal stretch stays inside
/// the trim interval, aligns as unresolved, and supports no variant, while an
/// SNV in the callable signal after it is reported.
#[test]
fn aligns_an_internal_masked_stretch_as_unresolved() -> Result<(), Box<dyn std::error::Error>> {
    let read = "ACGTCAGTACGATCGTACCTGAGTACGATCGATCGTAGCTGACTAGCTAGCATGACGTCAGTCATGCATCGATGCTAGCTAGTCGATCGA";
    let mut reference_read = read.to_owned();
    reference_read.replace_range(40..41, "C");
    reference_read.replace_range(70..71, "A");
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    write_abif_with_incoherent_doubles(&trace, read, 32..52, 500)?;
    write_reference(&reference, &format!("TTTT{reference_read}CCCC"))?;
    write_config(&config, "linear")?;

    run(&trace, &reference, &config, directory.path()).success();
    let value = read_result(directory.path(), &trace)?;
    let segments = value["read"]["callability"]["segments"]
        .as_array()
        .ok_or("segments must be an array")?;
    let masked = segments
        .iter()
        .find(|segment| segment["state"] != "in_phase")
        .ok_or("an internal masked segment is expected")?;
    assert_eq!(masked["state"], "mixed");
    let masked_calls = value["read"]["callability"]["masked_calls"]
        .as_u64()
        .ok_or("masked_calls must be a count")?;
    assert_eq!(value["read"]["trim"]["start"], 0);
    assert_eq!(value["read"]["trim"]["end"], read.len());
    assert_eq!(
        value["alignment"]["masked_bases"].as_u64(),
        Some(masked_calls)
    );
    let start = masked["calls"]["start"].as_u64().ok_or("segment start")?;
    let end = masked["calls"]["end"].as_u64().ok_or("segment end")?;
    assert_eq!(
        value["alignment"]["callable_reference_segments"],
        serde_json::json!([
            {"start": 4, "end": 4 + start},
            {"start": 4 + end, "end": 4 + read.len()}
        ])
    );
    let variants = value["variants"]
        .as_array()
        .ok_or("variants must be an array")?;
    assert_eq!(variants.len(), 1, "{variants:?}");
    assert_eq!(variants[0]["position"], 4 + 70 + 1);
    Ok(())
}

/// Read callability (ADR-0067): a slippage shadow ladder behind a long
/// homopolymer is published as a dephased segment attributed to the repeat,
/// while the in-phase prefix through the run's last call stays callable.
#[test]
fn publishes_a_dephased_segment_after_a_long_homopolymer() -> Result<(), Box<dyn std::error::Error>>
{
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let read = format!(
        "ACGTAGTCAGTACG{}TAGCTAGCATGCATGACTGACTAGCATGCA",
        "C".repeat(9)
    );
    write_abif_with_shadow_ladder(&trace, &read, 23, &[(-1, 0.4)])?;
    write_reference(&reference, &format!("TTTT{read}CCCC"))?;
    write_config(&config, "linear")?;

    run(&trace, &reference, &config, directory.path()).success();
    let value = read_result(directory.path(), &trace)?;
    let callability = &value["read"]["callability"];
    assert_eq!(
        callability["segments"],
        serde_json::json!([
            {"calls": {"start": 0, "end": 23}, "state": "in_phase", "after_repeat": false},
            {"calls": {"start": 23, "end": read.len()}, "state": "dephased", "after_repeat": true, "shadow_offsets": [-1]}
        ])
    );
    assert_eq!(callability["callable_span"]["start"], 0);
    assert_eq!(callability["callable_span"]["end"], 23);
    assert_eq!(callability["masked_calls"], read.len() - 23);
    let log = fs::read_to_string(directory.path().join("logs/trace.log"))?;
    assert!(log.contains("event=callability_completed"));
    assert!(log.contains(&format!(
        "repeats=1 segments=2 in_phase_segments=1 dephased_segments=1 mixed_segments=0 weak_segments=0 irregular_segments=0 masked_calls={} callable=0..23",
        read.len() - 23
    )));
    let segment_map = log
        .split_whitespace()
        .find(|field| field.starts_with("segment_map="))
        .ok_or("callability_completed must carry a segment map")?;
    assert!(
        segment_map.starts_with("segment_map=0..23:in_phase,23..")
            && segment_map.contains(":dephased(-1)[")
            && segment_map.ends_with("]+repeat"),
        "{segment_map}"
    );
    Ok(())
}

/// Two-sided slippage (shorter and longer length populations) behind a long
/// homopolymer is dephased with shadows on both sides, not mixed.
#[test]
fn publishes_two_sided_shadows_after_a_long_homopolymer() -> Result<(), Box<dyn std::error::Error>>
{
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let read = format!(
        "ACGTAGTCAGTACG{}TAGCTAGCATGCATGACTGACTAGCATGCA",
        "C".repeat(9)
    );
    write_abif_with_shadow_ladder(&trace, &read, 23, &[(-1, 0.4), (1, 0.4)])?;
    write_reference(&reference, &format!("TTTT{read}CCCC"))?;
    write_config(&config, "linear")?;

    run(&trace, &reference, &config, directory.path()).success();
    let value = read_result(directory.path(), &trace)?;
    let segments = value["read"]["callability"]["segments"]
        .as_array()
        .ok_or("segments must be an array")?;
    let last = segments.last().ok_or("segments must not be empty")?;
    assert_eq!(last["state"], "dephased");
    assert_eq!(last["after_repeat"], true);
    assert_eq!(last["shadow_offsets"], serde_json::json!([-1, 1]));
    let log = fs::read_to_string(directory.path().join("logs/trace.log"))?;
    assert!(log.contains(":dephased(-1,+1)["), "{log}");
    Ok(())
}

/// Double peaks that no neighbour offset explains are a mixed segment, and a
/// collapsed tail is weak; neither is attributed to a repeat.
#[test]
fn publishes_mixed_and_weak_segments_without_a_repeat() -> Result<(), Box<dyn std::error::Error>> {
    let read = "ACGTCAGTACGATCGTACCTGAGTACGATCGATCGTAGCTGACTAGCTAGCATGAC";
    let mixed = tempdir()?;
    let trace = mixed.path().join("trace.ab1");
    let reference = mixed.path().join("reference.fa");
    let config = mixed.path().join("dna.toml");
    write_abif_with_incoherent_doubles(&trace, read, 30..read.len(), 500)?;
    write_reference(&reference, &format!("TTTT{read}CCCC"))?;
    write_config(&config, "linear")?;
    run(&trace, &reference, &config, mixed.path()).success();
    let value = read_result(mixed.path(), &trace)?;
    let segments = value["read"]["callability"]["segments"]
        .as_array()
        .ok_or("segments must be an array")?;
    assert_eq!(segments.len(), 2);
    assert_eq!(segments[1]["state"], "mixed");
    assert_eq!(segments[1]["after_repeat"], false);
    assert!(segments[1].get("shadow_offsets").is_none());
    assert_eq!(segments[1]["calls"]["start"], 30);

    let weak = tempdir()?;
    let trace = weak.path().join("trace.ab1");
    let reference = weak.path().join("reference.fa");
    let config = weak.path().join("dna.toml");
    write_abif_with_amplitude_decay(&trace, read, 40, 0.001)?;
    write_reference(&reference, &format!("TTTT{read}CCCC"))?;
    write_config(&config, "linear")?;
    run(&trace, &reference, &config, weak.path()).success();
    let value = read_result(weak.path(), &trace)?;
    let segments = value["read"]["callability"]["segments"]
        .as_array()
        .ok_or("segments must be an array")?;
    assert_eq!(segments.len(), 2);
    assert_eq!(segments[1]["state"], "weak");
    assert_eq!(segments[1]["calls"]["start"], 40);
    assert_eq!(
        value["read"]["callability"]["masked_calls"],
        read.len() - 40
    );
    Ok(())
}

/// Regression: an indel flank whose two strongest channels tie is an unresolved
/// `N` with no primary event; it is omitted from public calls instead of
/// aborting the analysis.
#[test]
fn omits_unresolved_indel_flank_without_failing() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    write_abif_with_secondary_signal(&trace, QUERY, 13, b'A', 1000)?;
    write_reference(
        &reference,
        &format!("TTTT{}A{}CCCC", &QUERY[..14], &QUERY[14..]),
    )?;
    write_config(&config, "linear")?;

    run(&trace, &reference, &config, directory.path()).success();

    let value = read_result(directory.path(), &trace)?;
    let variants = value["variants"]
        .as_array()
        .ok_or("variants must be an array")?;
    assert_eq!(variants.len(), 1);
    assert_eq!(variants[0]["kind"], "DEL");
    let calls = variants[0]["calls"]
        .as_array()
        .ok_or("calls must be an array")?;
    assert_eq!(calls.len(), 1, "the unresolved flank is omitted");
    assert_eq!(calls[0]["role"], "flanking");
    assert!(matches!(
        calls[0]["base"].as_str(),
        Some("A" | "C" | "G" | "T")
    ));
    assert_call_evidence(&calls[0]);
    Ok(())
}

/// Regression: a flank whose four channels all qualify is called `N` despite a
/// primary event; the result schema allows only A/C/G/T call bases, so the
/// unresolved flank is omitted rather than published.
#[test]
fn omits_mixed_signal_indel_flank_called_n() -> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    write_abif_with_secondary_signals(
        &trace,
        QUERY,
        &[(13, b'A', 900), (13, b'G', 900), (13, b'T', 900)],
    )?;
    write_reference(
        &reference,
        &format!("TTTT{}A{}CCCC", &QUERY[..14], &QUERY[14..]),
    )?;
    write_config(&config, "linear")?;

    run(&trace, &reference, &config, directory.path()).success();

    let value = read_result(directory.path(), &trace)?;
    let variants = value["variants"]
        .as_array()
        .ok_or("variants must be an array")?;
    assert_eq!(variants.len(), 1);
    let calls = variants[0]["calls"]
        .as_array()
        .ok_or("calls must be an array")?;
    assert_eq!(calls.len(), 1, "the unresolved flank is omitted");
    assert!(
        calls
            .iter()
            .all(|call| matches!(call["base"].as_str(), Some("A" | "C" | "G" | "T")))
    );
    Ok(())
}

/// Every operation-log write fails (ENOSPC): the run fails fast on the first
/// record with the log error and publishes no result.
#[cfg(target_os = "linux")]
#[test]
fn unwritable_operation_log_fails_fast_without_publishing_a_result()
-> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    write_abif(&trace, QUERY)?;
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;
    let logs = directory.path().join("logs");
    fs::create_dir(&logs)?;
    std::os::unix::fs::symlink("/dev/full", logs.join("trace.log"))?;

    run(&trace, &reference, &config, directory.path())
        .failure()
        .stdout(predicate::str::is_empty())
        .stderr(predicate::str::starts_with(
            "error: failed to access log path logs/trace.log: No space left on device",
        ))
        .stderr(predicate::str::contains("additionally").not());
    assert!(!analysis_output_path(directory.path(), &trace).exists());
    Ok(())
}

fn assert_object_keys(value: &Value, expected: &[&str]) {
    let actual = value
        .as_object()
        .map(|object| object.keys().map(String::as_str).collect::<BTreeSet<_>>())
        .unwrap_or_default();
    let expected = expected.iter().copied().collect::<BTreeSet<_>>();
    assert_eq!(actual, expected);
}

fn read_result(workdir: &Path, trace: &Path) -> Result<Value, Box<dyn std::error::Error>> {
    Ok(serde_json::from_slice(&fs::read(analysis_output_path(
        workdir, trace,
    ))?)?)
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
            _ => 'N',
        })
        .collect()
}

fn assert_call_evidence(call: &Value) {
    assert_object_keys(call, &["role", "base", "peaks", "quality"]);
    assert_object_keys(&call["peaks"], &["A", "C", "G", "T"]);
    assert!(call["quality"].is_number());
}

fn run(
    trace: &Path,
    reference: &Path,
    config: &Path,
    workdir: &Path,
) -> assert_cmd::assert::Assert {
    let mut command = Command::new(dna_binary());
    command
        .current_dir(workdir)
        .env("DNA_CONFIG", config)
        .arg("analyze")
        .arg(trace)
        .arg("--reference")
        .arg(reference)
        .assert()
}
