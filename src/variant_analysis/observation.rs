//! One authoritative read-level reference-guided scientific path.

use crate::config::Config;
use crate::error::Result;
use crate::model::attachment::SangerAttachment;
use crate::model::read_observation::ReadObservation;
use crate::model::reference::Reference;
use crate::model::sanger::Chromatogram;
use crate::profile::Profile;
use crate::read_evidence::VetoSet;
use crate::read_processing::{self, ProcessedRead};

use crate::read_call::{self, CoreRun, ReadIdentity};

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
    let processed = read_processing::process(
        trace,
        &config.sanger,
        config.core.variant_calling.read_end_margin,
    )?;
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

    let stage = tracing::info_span!("read_evidence").entered();
    let evidence = read_processing::read_evidence(
        &calls,
        &signal,
        &callability,
        &quality,
        &config.sanger.sanger_evidence,
    )?;
    tracing::info!(
        event = "read_evidence_completed",
        calls = evidence.calls().len(),
        informative = ?evidence.informative(),
        masked_calls = evidence.calls().iter().filter(|call| call.mask.is_some()).count(),
        vetoed_calls = evidence
            .calls()
            .iter()
            .filter(|call| call.vetoes != VetoSet::default())
            .count(),
        minimum_peak_height = config.sanger.sanger_evidence.minimum_peak_height,
        relative_quality_threshold = config.sanger.sanger_evidence.relative_quality_threshold,
    );
    drop(stage);
    let called = read_call::call_read(
        ReadIdentity {
            input_name: trace.source_name.clone(),
            input_sha256: trace.source_sha256.clone(),
        },
        evidence,
        reference,
        &CoreRun {
            config: &config.core,
            configuration_sha256: &config.source_sha256,
            regions: &profile.regions,
        },
    )?;

    let excluded_variant_candidates = called.variants.excluded_count();
    let reference_origin_wrap = called.alignment.wraps_origin;
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
            called,
            sanger: SangerAttachment {
                calls,
                signal,
                callability,
                quality,
            },
        },
        warning_total,
    })
}
