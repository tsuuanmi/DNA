//! Shared multi-trace read processing for sample operations.

use crate::config::Config;
use crate::error::Result;
use crate::logger::Logger;
use crate::model::read_observation::ReadObservation;
use crate::model::reference::Reference;
use crate::model::trace::Chromatogram;

use super::observation;

pub(crate) struct CompletedSampleReads {
    pub(crate) reads: Vec<ReadObservation>,
    pub(crate) warning_total: usize,
}

pub(crate) fn build(
    traces: &[Chromatogram],
    reference: &Reference,
    config: &Config,
    logger: &mut Logger,
    stage: &mut &'static str,
) -> Result<CompletedSampleReads> {
    let mut reads = Vec::with_capacity(traces.len());
    let mut warning_total = 0usize;

    for (index, trace) in traces.iter().enumerate() {
        logger.info(
            module_path!(),
            line!(),
            format_args!(
                "event=sample_read_started read_index={} trace_name={:?} trace_sha256={}",
                index, trace.source_name, trace.source_sha256
            ),
        )?;
        let completed = observation::build(trace, reference, config, logger, stage)?;
        warning_total += completed.warning_total;
        logger.info(
            module_path!(),
            line!(),
            format_args!(
                "event=sample_read_completed read_index={} trace_sha256={} orientation={:?} segments={} variants={}",
                index,
                completed.read.input_sha256,
                completed.read.alignment.orientation,
                completed.read.alignment.reference_segments.len(),
                completed.read.variants.reported.len()
            ),
        )?;
        reads.push(completed.read);
    }

    Ok(CompletedSampleReads {
        reads,
        warning_total,
    })
}
