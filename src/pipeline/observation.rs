//! One authoritative read-level reference-guided scientific path.

use std::time::Instant;

use crate::alignment;
use crate::config::Config;
use crate::error::Result;
use crate::logger::Logger;
use crate::model::read_observation::ReadObservation;
use crate::model::reference::Reference;
use crate::model::trace::Chromatogram;
use crate::model::variant::VariantKind;
use crate::pipeline::read::{self, ProcessedRead};
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
    logger: &mut Logger,
    stage: &mut &'static str,
) -> Result<CompletedObservation> {
    let ProcessedRead {
        calls,
        signal,
        quality,
        warnings: read_warnings,
    } = read::process(trace, config, logger, stage)?;

    *stage = "alignment";
    let stage_started = Instant::now();
    let alignment = alignment::align_best(&quality, &signal, reference, &config.alignment)?;
    let reference_segments = alignment
        .reference_segments
        .iter()
        .map(|segment| format!("{}..{}", segment.start_0based, segment.end_0based_exclusive))
        .collect::<Vec<_>>()
        .join(",");
    logger.info(
        module_path!(),
        line!(),
        format_args!(
            concat!(
                "event=alignment_completed elapsed_ms={} orientation={:?} profile_score_units={} ",
                "exact_matches={} mismatches={} gap_opens={} callable_columns={} ",
                "callable_identity={:.4} unresolved_query_bases={} segments={} ",
                "segment_bounds={:?} wraps_origin={}"
            ),
            stage_started.elapsed().as_millis(),
            alignment.orientation,
            alignment.score,
            alignment.metrics.exact_matches,
            alignment.metrics.mismatches,
            alignment.metrics.gap_opens,
            alignment.metrics.callable_columns,
            alignment.metrics.callable_identity,
            alignment.metrics.unresolved_query_bases,
            alignment.reference_segments.len(),
            reference_segments,
            alignment.wraps_origin
        ),
    )?;

    *stage = "variant_calling";
    let stage_started = Instant::now();
    let variants = variant_calling::call(
        &alignment,
        reference,
        &calls,
        &quality,
        &config.variant_calling,
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
    logger.info(
        module_path!(),
        line!(),
        format_args!(
            concat!(
                "event=variant_calling_completed elapsed_ms={} reported={} snv={} insertion={} ",
                "deletion={} excluded={} region_count={} minimum_peak_height={} ",
                "relative_quality_threshold={} max_indel_length={}"
            ),
            stage_started.elapsed().as_millis(),
            variants.reported.len(),
            snvs,
            insertions,
            deletions,
            variants.excluded_count(),
            config.variant_calling.regions.len(),
            config.variant_calling.minimum_peak_height,
            config.variant_calling.relative_quality_threshold,
            config.variant_calling.max_indel_length
        ),
    )?;
    for excluded in &variants.excluded {
        let position = excluded
            .position_1based
            .map_or_else(|| "unknown".to_owned(), |position| position.to_string());
        let reasons = excluded
            .reasons
            .iter()
            .map(|reason| reason.label())
            .collect::<Vec<_>>()
            .join(",");
        logger.warn(
            module_path!(),
            line!(),
            format_args!(
                "event=variant_removed kind={} contig={:?} position={} reasons={}",
                excluded.kind.label(),
                excluded.contig,
                position,
                reasons
            ),
        )?;
    }

    let excluded_variant_candidates = variants.excluded_count();
    let reference_origin_wrap = alignment.wraps_origin;
    let warning_total = read_warnings.unresolved_primary_calls
        + read_warnings.multi_channel_unresolved_calls
        + read_warnings.vendor_disagreements
        + read_warnings.ploc_vendor_length_mismatches
        + read_warnings.clipped_channel_samples
        + excluded_variant_candidates
        + usize::from(reference_origin_wrap);
    if warning_total > 0 {
        logger.warn(
            module_path!(),
            line!(),
            format_args!(
                concat!(
                    "event=warning_summary total={} unresolved_primary_calls={} ",
                    "multi_channel_unresolved_calls={} vendor_disagreements={} ",
                    "ploc_vendor_length_mismatches={} clipped_channel_samples={} ",
                    "excluded_variant_candidates={} reference_origin_wrap={}"
                ),
                warning_total,
                read_warnings.unresolved_primary_calls,
                read_warnings.multi_channel_unresolved_calls,
                read_warnings.vendor_disagreements,
                read_warnings.ploc_vendor_length_mismatches,
                read_warnings.clipped_channel_samples,
                excluded_variant_candidates,
                reference_origin_wrap
            ),
        )?;
    }

    Ok(CompletedObservation {
        read: ReadObservation {
            input_name: trace.source_name.clone(),
            input_sha256: trace.source_sha256.clone(),
            reference_sha256: reference.sequence_sha256.clone(),
            configuration_sha256: config.source_sha256.clone(),
            calls,
            signal,
            quality,
            alignment,
            variants,
        },
        warning_total,
    })
}
