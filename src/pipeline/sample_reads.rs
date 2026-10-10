//! Shared multi-trace read processing for sample operations.

use crate::config::Config;
use crate::model::read_observation::{ReadObservation, RejectedRead, RejectionCause};
use crate::variant_analysis;
use dna_core::read_call::PlacementRejection;
use dna_kernel::error::Result;
use dna_kernel::model::reference::Reference;
use dna_kernel::profile::Profile;
use dna_sanger::model::attachment::SangerRejection;
use dna_sanger::model::sanger::Chromatogram;
use dna_sanger::read_processing;

pub(crate) struct CompletedSampleReads {
    pub(crate) reads: Vec<ReadObservation>,
    pub(crate) rejected: Vec<RejectedRead>,
    pub(crate) warning_total: usize,
}

/// Processes every trace independently; a read with too few callable calls,
/// or one the core cannot place, is recorded as rejected and the remaining
/// reads continue.
pub(crate) fn build(
    traces: &[Chromatogram],
    reference: &Reference,
    config: &Config,
    profile: &Profile,
) -> Result<CompletedSampleReads> {
    let mut reads = Vec::with_capacity(traces.len());
    let mut rejected = Vec::new();
    let mut warning_total = 0usize;

    for (index, trace) in traces.iter().enumerate() {
        tracing::info!(
            event = "sample_read_started",
            read_index = index,
            trace_name = ?trace.source_name,
            trace_sha256 = %trace.source_sha256,
        );
        let prepared = read_processing::prepare(trace, &config.sanger)?;
        if let Some(rejection) = prepared.rejection(&config.sanger) {
            tracing::warn!(
                event = "sample_read_rejected",
                read_index = index,
                trace_sha256 = %trace.source_sha256,
                reason = "callable_calls_below_minimum",
                callable_calls = rejection.callable_calls,
                minimum_callable_calls = rejection.minimum_callable_calls,
            );
            rejected.push(RejectedRead {
                sanger: SangerRejection {
                    input_name: trace.source_name.clone(),
                    input_sha256: trace.source_sha256.clone(),
                    integrity: prepared.signal.integrity,
                    callability: prepared.callability,
                },
                cause: RejectionCause::Callability(rejection),
            });
            continue;
        }
        let processed = read_processing::finish(
            trace,
            prepared,
            &config.sanger,
            config.core.variant_calling.read_end_margin,
        )?;
        let sanger = SangerRejection {
            input_name: trace.source_name.clone(),
            input_sha256: trace.source_sha256.clone(),
            integrity: processed.signal.integrity.clone(),
            callability: processed.callability.clone(),
        };
        let completed = match variant_analysis::observation::observe(
            trace, processed, reference, config, profile,
        ) {
            Ok(completed) => completed,
            Err(error) => {
                let Some(placement) = PlacementRejection::of(&error) else {
                    return Err(error);
                };
                tracing::warn!(
                    event = "sample_read_rejected",
                    read_index = index,
                    trace_sha256 = %trace.source_sha256,
                    reason = placement.label(),
                    detail = %error,
                );
                rejected.push(RejectedRead {
                    sanger,
                    cause: RejectionCause::Placement(placement),
                });
                continue;
            }
        };
        warning_total += completed.warning_total;
        tracing::info!(
            event = "sample_read_completed",
            read_index = index,
            trace_sha256 = %completed.read.called.input_sha256,
            orientation = ?completed.read.called.alignment.orientation,
            segments = completed.read.called.alignment.reference_segments.len(),
            variants = completed.read.called.variants.reported.len(),
        );
        reads.push(completed.read);
    }

    Ok(CompletedSampleReads {
        reads,
        rejected,
        warning_total,
    })
}
