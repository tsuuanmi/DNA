//! The core's one-read path: placement and variant calling from modality evidence.

use dna_kernel::plugin::{Contract, PluginDescriptor, PluginFamily};
use std::time::Instant;

use crate::alignment::{self, AlignmentConfig, RawAlignmentConfig};
use crate::model::called_read::CalledRead;
use crate::model::variant::VariantKind;
use crate::sample::{RawSampleReconciliationConfig, SampleReconciliationConfig};
use crate::variant_calling::{self, RawVariantCallingConfig, VariantCallingConfig};
use dna_kernel::error::Result;
use dna_kernel::model::reference::Reference;
use dna_kernel::read_evidence::ReadEvidence;

/// Every configuration section the core plugin owns.
#[derive(Debug, Clone)]
pub struct CoreConfig {
    pub(crate) alignment: AlignmentConfig,
    /// The `[variant_calling]` section.
    pub variant_calling: VariantCallingConfig,
    /// The `[sample_reconciliation]` section.
    pub sample_reconciliation: SampleReconciliationConfig,
}

/// The core plugin's sections as written in the configuration.
pub struct RawCoreConfig {
    /// The `[alignment]` section as written.
    pub alignment: RawAlignmentConfig,
    /// The `[variant_calling]` section as written.
    pub variant_calling: RawVariantCallingConfig,
    /// The `[sample_reconciliation]` section as written.
    pub sample_reconciliation: RawSampleReconciliationConfig,
}

impl RawCoreConfig {
    /// Validates every section.
    ///
    /// # Errors
    ///
    /// Returns `ConfigError` naming the first section value that violates its
    /// rule.
    pub fn validate(self) -> Result<CoreConfig> {
        Ok(CoreConfig {
            alignment: self.alignment.validate()?,
            variant_calling: self.variant_calling.validate()?,
            sample_reconciliation: self.sample_reconciliation.validate()?,
        })
    }
}

/// Run-level inputs of the core that do not come from the read.
pub struct CoreRun<'a> {
    /// The core's configuration sections.
    pub config: &'a CoreConfig,
    /// SHA-256 of the configuration the run was made with.
    pub configuration_sha256: &'a str,
    /// Inclusive 1-based reportable regions of the target.
    pub regions: &'a [[usize; 2]],
}

/// Content identity and reviewer-facing name of one read.
pub struct ReadIdentity {
    /// Reviewer-facing name of the read's source.
    pub input_name: String,
    /// SHA-256 content identity of the read's source.
    pub input_sha256: String,
}

/// Places one read's evidence on the reference and calls its variants; the
/// modality that produced the evidence is unknown here (ADR-0069).
///
/// # Errors
///
/// Returns `AlignmentError` when the read cannot be placed and `VariantError`
/// when its differences cannot be called.
pub fn call_read(
    identity: ReadIdentity,
    evidence: ReadEvidence,
    reference: &Reference,
    run: &CoreRun<'_>,
) -> Result<CalledRead> {
    let config = run.config;
    let stage = tracing::info_span!("alignment").entered();
    let stage_started = Instant::now();
    let alignment = alignment::align_best(&evidence, reference, &config.alignment)?;
    tracing::info!(
        event = "alignment_completed",
        elapsed_ms = stage_started.elapsed().as_millis(),
        orientation = ?alignment.orientation,
        profile_score_units = alignment.score,
        exact_matches = alignment.metrics.exact_matches,
        mismatches = alignment.metrics.mismatches,
        gap_opens = alignment.metrics.gap_opens,
        callable_columns = alignment.metrics.callable_columns,
        callable_identity = %format_args!("{:.4}", alignment.metrics.callable_identity),
        unresolved_query_bases = alignment.metrics.unresolved_query_bases,
        masked_query_bases = alignment.metrics.masked_query_bases,
        segments = alignment.reference_segments.len(),
        segment_bounds = ?alignment
            .reference_segments
            .iter()
            .map(|segment| format!("{}..{}", segment.start_0based, segment.end_0based_exclusive))
            .collect::<Vec<_>>()
            .join(","),
        wraps_origin = alignment.wraps_origin,
    );

    drop(stage);
    let _stage = tracing::info_span!("variant_calling").entered();
    let stage_started = Instant::now();
    let variants = variant_calling::call(
        &alignment,
        reference,
        &evidence,
        &config.variant_calling,
        run.regions,
    )?;
    let snvs = variants
        .reported
        .iter()
        .filter(|variant| variant.kind == VariantKind::Snv)
        .count();
    let insertions = variants
        .reported
        .iter()
        .filter(|variant| variant.kind == VariantKind::Ins)
        .count();
    let deletions = variants
        .reported
        .iter()
        .filter(|variant| variant.kind == VariantKind::Del)
        .count();
    tracing::info!(
        event = "variant_calling_completed",
        elapsed_ms = stage_started.elapsed().as_millis(),
        reported = variants.reported.len(),
        snv = snvs,
        insertion = insertions,
        deletion = deletions,
        excluded = variants.excluded_count(),
        region_count = run.regions.len(),
        max_indel_length = config.variant_calling.max_indel_length,
    );
    for excluded in &variants.excluded {
        tracing::warn!(
            event = "variant_removed",
            kind = excluded.kind.label(),
            contig = ?excluded.contig,
            position = %excluded
                .position_1based
                .map_or_else(|| "unknown".to_owned(), |position| position.to_string()),
            reasons = %excluded
                .reasons
                .iter()
                .map(|reason| reason.label())
                .collect::<Vec<_>>()
                .join(","),
        );
    }

    Ok(CalledRead {
        input_name: identity.input_name,
        input_sha256: identity.input_sha256,
        reference_sha256: reference.sequence_sha256.clone(),
        configuration_sha256: run.configuration_sha256.to_owned(),
        evidence,
        alignment,
        variants,
    })
}

/// Core caller: evidence-profile alignment, per-read variant calling, and
/// sample aggregation.
pub const PLUGIN: PluginDescriptor = PluginDescriptor {
    id: "core",
    family: PluginFamily::Core,
    version: 1,
    provides: &[Contract::CalledVariants],
    requires: &[Contract::ReadEvidence],
    config_sections: &["alignment", "sample_reconciliation", "variant_calling"],
};
