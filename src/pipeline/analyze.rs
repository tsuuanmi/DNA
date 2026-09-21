//! Single-read analysis orchestration, reporting, and publication.

use std::time::Instant;

use crate::cli::AnalyzeArgs;
use crate::error::Result;
use crate::logger::Logger;
use crate::pipeline::{input, observation};
use crate::report::{self, CompletedAnalysis};

/// Runs one complete AB1-to-JSON analysis with one per-trace append-only log.
pub(crate) fn run(args: &AnalyzeArgs) -> Result<()> {
    let trace_stem = input::trace_stem(&args.trace)?;
    let mut logger = Logger::open(trace_stem)?;
    let analysis_started = Instant::now();
    logger.info(
        module_path!(),
        line!(),
        format_args!(
            "event=analysis_started version={} trace_path={:?} reference_path={:?}",
            env!("CARGO_PKG_VERSION"),
            args.trace.display().to_string(),
            args.reference.display().to_string()
        ),
    )?;

    let mut stage = "input_loading";
    match run_logged(args, &mut logger, &mut stage, analysis_started) {
        Ok(()) => Ok(()),
        Err(error) => Err(super::record_failure(
            &mut logger,
            "analysis_failed",
            stage,
            analysis_started,
            error,
        )),
    }
}

fn run_logged(
    args: &AnalyzeArgs,
    logger: &mut Logger,
    stage: &mut &'static str,
    analysis_started: Instant,
) -> Result<()> {
    *stage = "input_loading";
    let stage_started = Instant::now();
    let inputs = input::load_analysis(args)?;
    logger.info(
        module_path!(),
        line!(),
        format_args!(
            concat!(
                "event=inputs_loaded elapsed_ms={} trace_name={:?} trace_sha256={} ",
                "samples={} call_loci={} vendor_primary={} vendor_quality={} ",
                "reference_name={:?} reference_sha256={} topology={:?} reference_bases={} ",
                "config_path={:?} config_sha256={} output_path={:?}"
            ),
            stage_started.elapsed().as_millis(),
            inputs.trace.source_name,
            inputs.trace.source_sha256,
            inputs.trace.sample_count(),
            inputs.trace.call_count(),
            inputs.trace.vendor.primary.is_some(),
            inputs.trace.vendor.qualities.is_some(),
            inputs.reference.name,
            inputs.reference.sequence_sha256,
            inputs.reference.topology,
            inputs.reference.len(),
            inputs.config.source_path.display().to_string(),
            inputs.config.source_sha256,
            inputs.output.display().to_string()
        ),
    )?;

    let completed = observation::build(
        &inputs.trace,
        &inputs.reference,
        &inputs.config,
        logger,
        stage,
    )?;

    *stage = "reporting";
    let stage_started = Instant::now();
    let output = inputs.output.clone();
    let warning_total = completed.warning_total;
    let result = report::build_analysis(CompletedAnalysis {
        reference: inputs.reference,
        read: completed.read,
    })?;
    let result_variants = result.variants.len();
    let schema_version = result.schema_version;
    let bytes = report::serialize(&result)?;

    *stage = "result_publication";
    logger.info(
        module_path!(),
        line!(),
        format_args!(
            concat!(
                "event=result_ready_for_publication elapsed_ms={} total_elapsed_ms={} ",
                "schema={} variants={} warnings={} output_path={:?} bytes={}"
            ),
            stage_started.elapsed().as_millis(),
            analysis_started.elapsed().as_millis(),
            schema_version,
            result_variants,
            warning_total,
            output.display().to_string(),
            bytes.len()
        ),
    )?;
    logger.sync()?;
    report::publish(&output, &bytes)
}
