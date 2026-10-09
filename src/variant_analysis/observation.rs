//! One authoritative read-level reference-guided scientific path.

use std::time::Instant;

use crate::alignment;
use crate::config::Config;
use crate::error::Result;
use crate::model::read_observation::ReadObservation;
use crate::model::reference::Reference;
use crate::model::sanger::Chromatogram;
use crate::model::variant::VariantKind;
use crate::profile::Profile;
use crate::read_processing::{self, ProcessedRead};
use crate::variant_calling;

/// Completed one-read observation plus operational warning total.
pub(crate) struct CompletedObservation {
    pub(crate) read: ReadObservation,
    pub(crate) warning_total: usize,
}

/// Runs the shared read, alignment, and variant stages for one trace.
pub(crate) fn build(
    trace: &Chromatogram,
    reference: &Reference,
    config: &Config,
    profile: &Profile,
) -> Result<CompletedObservation> {
    let processed = read_processing::process(trace, config)?;
    observe(trace, processed, reference, config, profile)
}

/// Runs the alignment and variant stages for one processed read.
pub(crate) fn observe(
    trace: &Chromatogram,
    processed: ProcessedRead,
    reference: &Reference,
    config: &Config,
    profile: &Profile,
) -> Result<CompletedObservation> {
    let ProcessedRead {
        calls,
        signal,
        callability,
        quality,
        warnings: read_warnings,
    } = processed;

    let stage = tracing::info_span!("alignment").entered();
    let stage_started = Instant::now();
    let evidence = read_processing::read_evidence(
        &calls,
        &signal,
        &callability,
        &quality,
        &config.sanger_evidence,
    )?;
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
        minimum_peak_height = config.sanger_evidence.minimum_peak_height,
        relative_quality_threshold = config.sanger_evidence.relative_quality_threshold,
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

    let excluded_variant_candidates = variants.excluded_count();
    let reference_origin_wrap = alignment.wraps_origin;
    let warning_total = read_warnings.unresolved_primary_calls
        + read_warnings.multi_channel_unresolved_calls
        + read_warnings.vendor_disagreements
        + read_warnings.locus_vendor_length_mismatches
        + read_warnings.clipped_channel_samples
        + excluded_variant_candidates
        + usize::from(reference_origin_wrap);
    if warning_total > 0 {
        tracing::warn!(
            event = "warning_summary",
            total = warning_total,
            unresolved_primary_calls = read_warnings.unresolved_primary_calls,
            multi_channel_unresolved_calls = read_warnings.multi_channel_unresolved_calls,
            vendor_disagreements = read_warnings.vendor_disagreements,
            ploc_vendor_length_mismatches = read_warnings.locus_vendor_length_mismatches,
            clipped_channel_samples = read_warnings.clipped_channel_samples,
            excluded_variant_candidates,
            reference_origin_wrap,
        );
    }

    Ok(CompletedObservation {
        read: ReadObservation {
            input_name: trace.source_name.clone(),
            input_sha256: trace.source_sha256.clone(),
            reference_sha256: reference.sequence_sha256.clone(),
            configuration_sha256: config.source_sha256.clone(),
            calls,
            signal,
            callability,
            quality,
            alignment,
            variants,
        },
        warning_total,
    })
}
