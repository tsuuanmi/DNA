//! Compact typed assembly and deterministic JSON serialization.

use crate::error::Result;
use crate::model::basecalls::BaseCalls;
use crate::model::read_observation::ReadObservation;
use crate::model::reference::Reference;
use crate::model::result::{
    AlignmentResult, AnalysisResult, InputResult, IntervalResult, ProvenanceResult, ReadResult,
    ReferenceResult, WarningSummaryResult,
};
use crate::report::{signal, variant};

/// Inputs consumed to build the immutable analysis document.
pub(crate) struct CompletedAnalysis {
    pub(crate) reference: Reference,
    pub(crate) read: ReadObservation,
}

/// Builds the compact v7 document without filesystem side effects.
pub(crate) fn build_analysis(completed: CompletedAnalysis) -> Result<AnalysisResult> {
    let CompletedAnalysis { reference, read } = completed;
    let ReadObservation {
        input_name: _,
        input_sha256,
        reference_sha256,
        configuration_sha256,
        calls,
        signal,
        quality,
        alignment,
        variants,
    } = read;
    if reference_sha256 != reference.sequence_sha256 {
        return Err(crate::error::Error::Report(
            "read observation reference identity does not match report reference".into(),
        ));
    }
    let warnings = warning_summary(&calls, &signal, variants.excluded_count());
    let variant_results =
        variant::project(variants.reported, &calls, &quality, alignment.orientation)?;
    let signal_quality = signal::project(signal);
    let reference_segments = alignment
        .reference_segments
        .into_iter()
        .map(|segment| IntervalResult {
            start: segment.start_0based,
            end: segment.end_0based_exclusive,
        })
        .collect();

    Ok(AnalysisResult {
        schema_version: "dna.analysis/v7",
        provenance: ProvenanceResult {
            input: InputResult {
                sha256: input_sha256,
            },
            reference: ReferenceResult {
                name: reference.name,
                topology: reference.topology,
                sha256: reference.sequence_sha256,
            },
            configuration_sha256,
        },
        read: ReadResult {
            call_count: calls.len(),
            trim: IntervalResult {
                start: quality.trim_start_0based,
                end: quality.trim_end_0based_exclusive,
            },
        },
        signal_quality,
        alignment: AlignmentResult {
            orientation: alignment.orientation,
            callable_bases: alignment.metrics.callable_columns,
            identity: alignment.metrics.callable_identity,
            unresolved_bases: alignment.metrics.unresolved_query_bases,
            gap_opens: alignment.metrics.gap_opens,
            reference_segments,
            wraps_origin: alignment.wraps_origin,
        },
        variants: variant_results,
        warnings,
    })
}

/// Serializes any typed result with a trailing newline for stable text files.
pub(crate) fn serialize<T: serde::Serialize>(result: &T) -> Result<Vec<u8>> {
    let mut bytes = serde_json::to_vec_pretty(result)?;
    bytes.push(b'\n');
    Ok(bytes)
}

fn warning_summary(
    calls: &BaseCalls,
    signal: &crate::model::signal::DNAAnalysis,
    excluded_variant_candidates: usize,
) -> WarningSummaryResult {
    let unresolved_primary_calls = calls
        .calls
        .iter()
        .filter(|call| call.primary == 'N')
        .count();
    let multi_channel_unresolved_calls = calls
        .calls
        .iter()
        .filter(|call| call.ambiguity == 'N' && call.qualifying_channels.len() > 2)
        .count();
    WarningSummaryResult {
        unresolved_primary_calls,
        multi_channel_unresolved_calls,
        ploc_vendor_length_mismatches: signal.integrity.vendor_length_mismatch_count(),
        clipped_channel_samples: signal.integrity.clipped_channel_samples,
        excluded_variant_candidates,
    }
}
