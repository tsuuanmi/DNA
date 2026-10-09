//! Projection of the core's called reads into `dna.variants/v1`.

use crate::model::result::{AlignmentResult, IntervalResult, ReferenceResult};
use crate::model::sample_result::SampleProvenanceResult;
use crate::model::variants_result::{ObservedVariantResult, ReadVariantsResult, VariantsResult};
use crate::report::json::{project_plugins, project_profile};
use crate::report::notation::{self, SampleNotation};
use dna_core::model::called_read::CalledRead;
use dna_kernel::error::{ReportError, Result};
use dna_kernel::model::reference::Reference;
use dna_kernel::plugin::PluginDescriptor;
use dna_kernel::profile::ProfileIdentity;

/// Inputs consumed to build one immutable variants document.
pub(crate) struct CompletedCall {
    pub(crate) sample_id: String,
    pub(crate) reference: Reference,
    pub(crate) profile: ProfileIdentity,
    /// Called reads in input order.
    pub(crate) reads: Vec<CalledRead>,
    /// Per-read represented calls, present only when the profile declares notation.
    pub(crate) notation: Option<SampleNotation>,
    /// Plugins of the workflow, in execution order.
    pub(crate) plugins: &'static [&'static PluginDescriptor],
}

/// Builds `dna.variants/v1` without filesystem side effects.
pub(crate) fn build(completed: CompletedCall) -> Result<VariantsResult> {
    let CompletedCall {
        sample_id,
        reference,
        profile,
        reads,
        notation,
        plugins,
    } = completed;
    let first = reads.first().ok_or(ReportError::Inconsistent(
        "a call document needs at least one read",
    ))?;
    let configuration_sha256 = first.configuration_sha256.clone();
    if reads.iter().any(|read| {
        read.reference_sha256 != reference.sequence_sha256
            || read.configuration_sha256 != configuration_sha256
    }) {
        return Err(ReportError::Inconsistent(
            "called reads disagree on reference or configuration identity",
        )
        .into());
    }
    let names = reads
        .iter()
        .map(|read| read.input_name.clone())
        .collect::<Vec<_>>();
    let notation = notation
        .map(|notation| {
            let identities = reads
                .iter()
                .map(|read| read.input_sha256.as_str())
                .collect::<Vec<_>>();
            notation::project(&reference, &identities, &names, notation)
        })
        .transpose()?;
    let reads = reads.into_iter().map(project_read).collect();
    Ok(VariantsResult {
        schema_version: "dna.variants/v1",
        sample_id,
        provenance: SampleProvenanceResult {
            reference: ReferenceResult {
                name: reference.name,
                topology: reference.topology,
                sha256: reference.sequence_sha256,
            },
            configuration_sha256,
            profile: project_profile(profile),
            plugins: project_plugins(plugins),
        },
        reads,
        notation,
    })
}

fn project_read(read: CalledRead) -> ReadVariantsResult {
    let intervals = |segments: Vec<dna_core::model::alignment::ReferenceSegment>| {
        segments
            .into_iter()
            .map(|segment| IntervalResult {
                start: segment.start_0based,
                end: segment.end_0based_exclusive,
            })
            .collect()
    };
    let alignment = read.alignment;
    let mut variants = read
        .variants
        .observed
        .into_iter()
        .map(|observed| ObservedVariantResult {
            eligible: observed.eligible(),
            position: observed.variant.position_1based,
            reference: observed.variant.reference,
            alternate: observed.variant.alternate,
            kind: observed.variant.kind,
            exclusion_reasons: observed.exclusion_reasons,
        })
        .collect::<Vec<_>>();
    variants.sort_by(|left, right| {
        (left.position, &left.reference, &left.alternate).cmp(&(
            right.position,
            &right.reference,
            &right.alternate,
        ))
    });
    ReadVariantsResult {
        name: read.input_name,
        sha256: read.input_sha256,
        alignment: AlignmentResult {
            orientation: alignment.orientation,
            callable_bases: alignment.metrics.callable_columns,
            identity: alignment.metrics.callable_identity,
            gap_opens: alignment.metrics.gap_opens,
            unresolved_bases: alignment.metrics.unresolved_query_bases,
            masked_bases: alignment.metrics.masked_query_bases,
            reference_segments: intervals(alignment.reference_segments),
            callable_reference_segments: intervals(alignment.callable_segments),
            wraps_origin: alignment.wraps_origin,
        },
        variants,
    }
}
