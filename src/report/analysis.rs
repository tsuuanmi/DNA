//! The analysis document, `dna.analysis/v10`: one read's compact result,
//! assembled from the core's products and the read's Sanger evidence.

use crate::error::{ReportError, Result};
use crate::model::result::{
    AlignmentResult, AnalysisResult, InputResult, IntervalResult, ProvenanceResult, ReadResult,
    ReferenceResult, WarningSummaryResult,
};
use crate::model::sanger_read::SangerRead;
use crate::report::json::{project_plugins, project_profile};
use crate::report::{aligned_calls, variant};
use dna_call::model::called_read::CalledRead;
use dna_call::model::trust::CallTrust;
use dna_kernel::model::reference::Reference;
use dna_kernel::plugin::PluginDescriptor;
use dna_kernel::profile::ProfileIdentity;

/// Inputs consumed to build the immutable analysis document.
pub(crate) struct CompletedAnalysis {
    pub(crate) reference: Reference,
    pub(crate) profile: ProfileIdentity,
    /// The core's products for the read.
    pub(crate) called: CalledRead,
    /// The read's Sanger evidence.
    pub(crate) sanger: SangerRead,
    /// The core's trust in every call of the read, in call order.
    pub(crate) trust: Vec<CallTrust>,
    /// Plugins of the workflow, in execution order.
    pub(crate) plugins: &'static [&'static PluginDescriptor],
}

/// Builds the compact v10 document without filesystem side effects.
pub(crate) fn build_analysis(completed: CompletedAnalysis) -> Result<AnalysisResult> {
    let CompletedAnalysis {
        reference,
        profile,
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
        trust,
        plugins,
    } = completed;
    if reference_sha256 != reference.sequence_sha256 {
        return Err(ReportError::Inconsistent(
            "read observation reference identity does not match report reference",
        )
        .into());
    }
    if trust.len() != evidence.calls().len() {
        return Err(
            ReportError::Inconsistent("analysis trust and evidence differ in call count").into(),
        );
    }
    if sanger.calls.len() != evidence.calls().len() {
        return Err(ReportError::Inconsistent(
            "analysis Sanger evidence and core evidence differ in call count",
        )
        .into());
    }
    let calls_result = aligned_calls::project(alignment.columns(), &trust)?;
    let warnings = WarningSummaryResult {
        unresolved_primary_calls: sanger.warnings.unresolved_primary_calls,
        multi_channel_unresolved_calls: sanger.warnings.multi_channel_unresolved_calls,
        ploc_vendor_length_mismatches: sanger.warnings.ploc_vendor_length_mismatches,
        clipped_channel_samples: sanger.warnings.clipped_channel_samples,
        excluded_variant_candidates: variants.excluded_count(),
    };
    let variant_results = variant::project(
        variants.reported,
        &evidence,
        &sanger.calls,
        alignment.orientation(),
    )?;
    let SangerRead {
        trim,
        callability,
        signal_quality,
        warnings: _,
        calls,
    } = sanger;
    let reference_segments = alignment
        .reference_segments()
        .iter()
        .map(|segment| IntervalResult {
            start: segment.start_0based,
            end: segment.end_0based_exclusive,
        })
        .collect();
    let callable_reference_segments = alignment
        .callable_segments()
        .iter()
        .map(|segment| IntervalResult {
            start: segment.start_0based,
            end: segment.end_0based_exclusive,
        })
        .collect();

    Ok(AnalysisResult {
        schema_version: "dna.analysis/v10",
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
            trim,
            callability,
        },
        signal_quality,
        alignment: AlignmentResult {
            orientation: alignment.orientation(),
            callable_bases: alignment.metrics().callable_columns,
            identity: alignment.metrics().callable_identity,
            unresolved_bases: alignment.metrics().unresolved_query_bases,
            masked_bases: alignment.metrics().masked_query_bases,
            gap_opens: alignment.metrics().gap_opens,
            reference_segments,
            callable_reference_segments,
            wraps_origin: alignment.wraps_origin(),
        },
        calls: calls_result,
        variants: variant_results,
        warnings,
    })
}
