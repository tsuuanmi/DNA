//! Multi-read sample evidence orchestration and publication.

use std::time::Instant;

use crate::cli::SampleArgs;
use crate::error::Result;
use crate::logger::Logger;
use crate::pipeline::input;
use crate::report::{self, CompletedSampleEvidence};
use crate::sample as sample_science;

use super::{sample_metrics, sample_reads};

/// Runs one sample-evidence operation with one sample-level append-only log.
pub(crate) fn run(args: &SampleArgs) -> Result<()> {
    input::validate_sample_id(&args.sample_id)?;
    let mut logger = Logger::open(&args.sample_id)?;
    let started = Instant::now();
    logger.info(
        module_path!(),
        line!(),
        format_args!(
            "event=sample_started version={} sample_id={:?} traces={} reference_path={:?}",
            env!("CARGO_PKG_VERSION"),
            args.sample_id,
            args.traces.len(),
            args.reference.display().to_string()
        ),
    )?;

    let mut stage = "input_loading";
    match run_logged(args, &mut logger, &mut stage, started) {
        Ok(()) => Ok(()),
        Err(error) => Err(super::record_failure(
            &mut logger,
            "sample_failed",
            stage,
            started,
            error,
        )),
    }
}

fn run_logged(
    args: &SampleArgs,
    logger: &mut Logger,
    stage: &mut &'static str,
    started: Instant,
) -> Result<()> {
    *stage = "input_loading";
    let stage_started = Instant::now();
    let inputs = input::load_sample(args)?;
    let output = input::sample_output(&args.sample_id)?;
    logger.info(
        module_path!(),
        line!(),
        format_args!(
            concat!(
                "event=sample_inputs_loaded elapsed_ms={} sample_id={:?} traces={} ",
                "reference_name={:?} reference_sha256={} topology={:?} reference_bases={} ",
                "config_path={:?} config_sha256={} output_path={:?}"
            ),
            stage_started.elapsed().as_millis(),
            args.sample_id,
            inputs.traces.len(),
            inputs.reference.name,
            inputs.reference.sequence_sha256,
            inputs.reference.topology,
            inputs.reference.len(),
            inputs.config.source_path.display().to_string(),
            inputs.config.source_sha256,
            output.display().to_string()
        ),
    )?;

    let completed_reads = sample_reads::build(
        &inputs.traces,
        &inputs.reference,
        &inputs.config,
        logger,
        stage,
    )?;
    let reads = completed_reads.reads;
    let warning_total = completed_reads.warning_total;

    *stage = "sample_aggregation";
    let stage_started = Instant::now();
    let evidence = sample_science::aggregate(&reads, &inputs.config.sample_reconciliation)?;
    let metrics = sample_metrics::summarize(&evidence);
    logger.info(
        module_path!(),
        line!(),
        format_args!(
            concat!(
                "event=sample_aggregation_completed elapsed_ms={} reads={} coverage_segments={} ",
                "overlaps={} eligible_overlaps={} locus_differences={} profiled_locus_observations={} ",
                "profiled_locus_forward_reads={} profiled_locus_reverse_reads={} noisy_locus_observations={} ",
                "eligible_nucleotide_locus_observations={} missing_profile_locus_observations={} ",
                "deletion_event_locus_observations={} nucleotide_support_loci={} ",
                "bidirectional_nucleotide_support_loci={} unweighted_nucleotide_profile_mass={:.6} ",
                "profile_geometry_loci={} within_profile_impurity_sum={:.6} ",
                "between_profile_dispersion_sum={:.6} total_profile_heterogeneity_sum={:.6} ",
                "forward_profile_geometry_loci={} reverse_profile_geometry_loci={} ",
                "directional_profile_distance_loci={} directional_profile_distance_sum={:.6} ",
                "locus_positive_corrected_channels={} locus_positive_snr_channels={} ",
                "locus_forward_reads={} locus_reverse_reads={} locus_reference_reads={} ",
                "locus_alternate_reads={} locus_unresolved_reads={} locus_deletion_reads={} variants={} ",
                "profiled_variant_calls={} noisy_variant_calls={} variant_positive_corrected_channels={} ",
                "variant_positive_snr_channels={}"
            ),
            stage_started.elapsed().as_millis(),
            evidence.reads.len(),
            evidence.coverage.len(),
            evidence.overlaps.len(),
            evidence
                .overlaps
                .iter()
                .filter(|overlap| overlap.eligible)
                .count(),
            evidence.locus_differences.len(),
            metrics.profiled_locus_observations,
            metrics.profiled_locus_forward_reads,
            metrics.profiled_locus_reverse_reads,
            metrics.noisy_locus_observations,
            metrics.eligible_nucleotide_locus_observations,
            metrics.missing_profile_locus_observations,
            metrics.deletion_event_locus_observations,
            metrics.nucleotide_support_loci,
            metrics.bidirectional_nucleotide_support_loci,
            metrics.unweighted_nucleotide_profile_mass,
            metrics.profile_geometry_loci,
            metrics.within_profile_impurity_sum,
            metrics.between_profile_dispersion_sum,
            metrics.total_profile_heterogeneity_sum,
            metrics.forward_profile_geometry_loci,
            metrics.reverse_profile_geometry_loci,
            metrics.directional_profile_distance_loci,
            metrics.directional_profile_distance_sum,
            metrics.locus_positive_corrected_channels,
            metrics.locus_positive_snr_channels,
            metrics.locus_forward_reads,
            metrics.locus_reverse_reads,
            metrics.locus_reference_reads,
            metrics.locus_alternate_reads,
            metrics.locus_unresolved_reads,
            metrics.locus_deletion_reads,
            evidence.variants.len(),
            metrics.profiled_variant_calls,
            metrics.noisy_variant_calls,
            metrics.variant_positive_corrected_channels,
            metrics.variant_positive_snr_channels
        ),
    )?;

    *stage = "reporting";
    let stage_started = Instant::now();
    let result = report::build_sample(CompletedSampleEvidence {
        sample_id: args.sample_id.clone(),
        reference: inputs.reference,
        evidence,
    })?;
    let reads = result.reads.len();
    let coverage_segments = result.coverage.len();
    let overlaps = result.overlaps.len();
    let locus_differences = result.locus_differences.len();
    let variants = result.variants.len();
    let schema_version = result.schema_version;
    let bytes = report::serialize(&result)?;

    *stage = "result_publication";
    logger.info(
        module_path!(),
        line!(),
        format_args!(
            concat!(
                "event=sample_result_ready_for_publication elapsed_ms={} total_elapsed_ms={} ",
                "schema={} reads={} coverage_segments={} overlaps={} locus_differences={} variants={} ",
                "read_warnings={} output_path={:?} bytes={}"
            ),
            stage_started.elapsed().as_millis(),
            started.elapsed().as_millis(),
            schema_version,
            reads,
            coverage_segments,
            overlaps,
            locus_differences,
            variants,
            warning_total,
            output.display().to_string(),
            bytes.len()
        ),
    )?;
    logger.sync()?;
    report::publish(&output, &bytes)
}
