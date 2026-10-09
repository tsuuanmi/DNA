//! Compact typed assembly and deterministic JSON serialization.

use crate::model::read_observation::ReadObservation;
use crate::model::result::{
    AlignmentResult, AnalysisResult, InputResult, IntervalResult, PluginResult, ProfileResult,
    ProvenanceResult, ReadResult, ReferenceResult, WarningSummaryResult,
};
use crate::report::{callability, signal, variant};
use dna_core::model::called_read::CalledRead;
use dna_kernel::error::{Error, ReportError, Result};
use dna_kernel::model::reference::Reference;
use dna_kernel::plugin::PluginDescriptor;
use dna_kernel::profile::ProfileIdentity;
use dna_sanger::model::attachment::SangerAttachment;
use dna_sanger::model::basecalls::BaseCalls;

/// Inputs consumed to build the immutable analysis document.
pub(crate) struct CompletedAnalysis {
    pub(crate) reference: Reference,
    pub(crate) profile: ProfileIdentity,
    pub(crate) read: ReadObservation,
    /// Plugins of the workflow, in execution order.
    pub(crate) plugins: &'static [&'static PluginDescriptor],
}

/// Builds the compact v9 document without filesystem side effects.
pub(crate) fn build_analysis(completed: CompletedAnalysis) -> Result<AnalysisResult> {
    let CompletedAnalysis {
        reference,
        profile,
        read,
        plugins,
    } = completed;
    let ReadObservation {
        called:
            CalledRead {
                input_name: _,
                input_sha256,
                reference_sha256,
                configuration_sha256,
                evidence,
                alignment,
                variants,
            },
        sanger,
    } = read;
    if reference_sha256 != reference.sequence_sha256 {
        return Err(ReportError::Inconsistent(
            "read observation reference identity does not match report reference",
        )
        .into());
    }
    let warnings = warning_summary(&sanger.calls, &sanger.signal, variants.excluded_count());
    let variant_results =
        variant::project(variants.reported, &evidence, &sanger, alignment.orientation)?;
    let SangerAttachment {
        calls,
        signal,
        callability: read_callability,
        quality,
    } = sanger;
    let signal_quality = signal::project(signal);
    let reference_segments = alignment
        .reference_segments
        .into_iter()
        .map(|segment| IntervalResult {
            start: segment.start_0based,
            end: segment.end_0based_exclusive,
        })
        .collect();
    let callable_reference_segments = alignment
        .callable_segments
        .into_iter()
        .map(|segment| IntervalResult {
            start: segment.start_0based,
            end: segment.end_0based_exclusive,
        })
        .collect();

    Ok(AnalysisResult {
        schema_version: "dna.analysis/v9",
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
            profile: project_profile(profile),
            plugins: project_plugins(plugins),
        },
        read: ReadResult {
            call_count: calls.len(),
            trim: IntervalResult {
                start: quality.trim_start_0based,
                end: quality.trim_end_0based_exclusive,
            },
            callability: callability::project(&read_callability),
        },
        signal_quality,
        alignment: AlignmentResult {
            orientation: alignment.orientation,
            callable_bases: alignment.metrics.callable_columns,
            identity: alignment.metrics.callable_identity,
            unresolved_bases: alignment.metrics.unresolved_query_bases,
            masked_bases: alignment.metrics.masked_query_bases,
            gap_opens: alignment.metrics.gap_opens,
            reference_segments,
            callable_reference_segments,
            wraps_origin: alignment.wraps_origin,
        },
        variants: variant_results,
        warnings,
    })
}

/// Projects the workflow's plugin identities recorded in result provenance.
pub(crate) fn project_plugins(plugins: &[&PluginDescriptor]) -> Vec<PluginResult> {
    plugins
        .iter()
        .map(|plugin| PluginResult {
            id: plugin.id,
            family: plugin.family.label(),
            version: plugin.version,
        })
        .collect()
}

/// Projects the profile identity recorded in result provenance.
pub(crate) fn project_profile(profile: ProfileIdentity) -> ProfileResult {
    ProfileResult {
        id: profile.id,
        sha256: profile.sha256,
    }
}

/// Serializes any typed result with a trailing newline for stable text files.
pub(crate) fn serialize<T: serde::Serialize>(result: &T) -> Result<Vec<u8>> {
    let mut bytes =
        serde_json::to_vec_pretty(result).map_err(|error| Error::Serialize(Box::new(error)))?;
    bytes.push(b'\n');
    Ok(bytes)
}

fn warning_summary(
    calls: &BaseCalls,
    signal: &dna_sanger::model::signal::SignalAnalysis,
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
