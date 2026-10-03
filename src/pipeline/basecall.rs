//! Reference-free base re-calling, quality trimming, and result publication.

use std::path::Path;
use std::time::Instant;

use crate::error::Result;
use crate::input::sanger;
use crate::logger::Logger;
use crate::pipeline::path;
use crate::read_processing::{self, ProcessedRead};
use crate::report::{self, CompletedBasecall};

/// Runs one complete AB1-to-basecalls JSON operation.
pub(crate) fn run(trace: &Path, config_path: &Path) -> Result<()> {
    let trace_stem = path::trace_stem(trace)?;
    let mut logger = Logger::open(trace_stem)?;
    let started = Instant::now();
    logger.info(
        module_path!(),
        line!(),
        format_args!(
            "event=basecall_started version={} trace_path={:?}",
            env!("CARGO_PKG_VERSION"),
            trace.display().to_string()
        ),
    )?;

    let mut stage = "input_loading";
    match run_logged(trace, config_path, &mut logger, &mut stage, started) {
        Ok(()) => Ok(()),
        Err(error) => Err(super::record_failure(
            &mut logger,
            "basecall_failed",
            stage,
            started,
            error,
        )),
    }
}

fn run_logged(
    trace: &Path,
    config_path: &Path,
    logger: &mut Logger,
    stage: &mut &'static str,
    started: Instant,
) -> Result<()> {
    *stage = "input_loading";
    let stage_started = Instant::now();
    let prepared = sanger::prepare_basecall(trace, config_path)?;
    let output = path::basecall_output(trace)?;
    let inputs = sanger::load_basecall(prepared)?;
    logger.info(
        module_path!(),
        line!(),
        format_args!(
            concat!(
                "event=basecall_inputs_loaded elapsed_ms={} trace_name={:?} trace_sha256={} ",
                "samples={} call_loci={} vendor_primary={} vendor_quality={} ",
                "config_path={:?} config_sha256={} output_path={:?}"
            ),
            stage_started.elapsed().as_millis(),
            inputs.trace.source_name,
            inputs.trace.source_sha256,
            inputs.trace.sample_count(),
            inputs.trace.call_count(),
            inputs.trace.vendor.primary.is_some(),
            inputs.trace.vendor.qualities.is_some(),
            inputs.config.source_path.display().to_string(),
            inputs.config.source_sha256,
            output.display().to_string()
        ),
    )?;

    let ProcessedRead {
        calls,
        signal,
        quality,
        warnings,
    } = read_processing::process(&inputs.trace, &inputs.config, logger, stage)?;
    let warning_total = warnings.unresolved_primary_calls
        + warnings.multi_channel_unresolved_calls
        + warnings.vendor_disagreements
        + warnings.locus_vendor_length_mismatches
        + warnings.clipped_channel_samples;
    if warning_total > 0 {
        logger.warn(
            module_path!(),
            line!(),
            format_args!(
                concat!(
                    "event=basecall_warning_summary total={} unresolved_primary_calls={} ",
                    "multi_channel_unresolved_calls={} vendor_disagreements={} ",
                    "ploc_vendor_length_mismatches={} clipped_channel_samples={}"
                ),
                warning_total,
                warnings.unresolved_primary_calls,
                warnings.multi_channel_unresolved_calls,
                warnings.vendor_disagreements,
                warnings.locus_vendor_length_mismatches,
                warnings.clipped_channel_samples
            ),
        )?;
    }

    *stage = "reporting";
    let stage_started = Instant::now();
    let result = report::build_basecall(CompletedBasecall {
        config: inputs.config,
        trace: inputs.trace,
        calls,
        signal,
        quality,
    })?;
    let schema_version = result.schema_version;
    let call_count = result.read.call_count;
    let retained = result.read.retained.len();
    let bytes = report::serialize(&result)?;

    *stage = "result_publication";
    logger.info(
        module_path!(),
        line!(),
        format_args!(
            concat!(
                "event=basecall_ready_for_publication elapsed_ms={} total_elapsed_ms={} ",
                "schema={} calls={} retained={} warnings={} output_path={:?} bytes={}"
            ),
            stage_started.elapsed().as_millis(),
            started.elapsed().as_millis(),
            schema_version,
            call_count,
            retained,
            warning_total,
            output.display().to_string(),
            bytes.len()
        ),
    )?;
    logger.sync()?;
    report::publish(&output, &bytes)
}
