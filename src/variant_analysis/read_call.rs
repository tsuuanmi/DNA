//! The core's one-read path: placement and variant calling from modality evidence.

use std::time::Instant;

use crate::alignment;
use crate::config::Config;
use crate::error::Result;
use crate::model::called_read::CalledRead;
use crate::model::reference::Reference;
use crate::model::variant::VariantKind;
use crate::profile::Profile;
use crate::read_evidence::ReadEvidence;
use crate::variant_calling;

/// Content identity and reviewer-facing name of one read.
pub(crate) struct ReadIdentity {
    pub(crate) input_name: String,
    pub(crate) input_sha256: String,
}

/// Places one read's evidence on the reference and calls its variants; the
/// modality that produced the evidence is unknown here (ADR-0069).
pub(crate) fn call_read(
    identity: ReadIdentity,
    evidence: ReadEvidence,
    reference: &Reference,
    config: &Config,
    profile: &Profile,
) -> Result<CalledRead> {
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
        &profile.regions,
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
        region_count = profile.regions.len(),
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
        configuration_sha256: config.source_sha256.clone(),
        evidence,
        alignment,
        variants,
    })
}
