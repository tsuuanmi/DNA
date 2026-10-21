//! One authoritative read-level reference-guided scientific path.

use crate::config::Config;
use crate::error::{Error, Result};
use crate::model::read_observation::ReadObservation;
use dna_kernel::model::reference::Reference;
use dna_kernel::profile::Profile;
use dna_kernel::read_evidence::{Placement, ReadEvidence};
use dna_sanger::model::attachment::{SangerAttachment, SangerRejection};
use dna_sanger::model::sanger::Chromatogram;
use dna_sanger::read_processing::{self, ProcessedRead, ReadWarnings};

use dna_align::model::alignment::Alignment;
use dna_call::model::called_read::CalledRead;
use dna_call::read_call::{self, CoreRun, PlacementRejection, ReadIdentity};

/// Completed one-read observation plus operational warning total.
pub(crate) struct CompletedObservation {
    pub(crate) read: ReadObservation,
    pub(crate) warning_total: usize,
}

/// What the core made of one processed read.
pub(crate) enum Observation {
    /// The core placed the read.
    Placed(Box<CompletedObservation>),
    /// The core could not place the read (SRS-SAMPLE-029); its integrity and
    /// callability are kept for a rejection record.
    Unplaced {
        rejection: PlacementRejection,
        error: Error,
        sanger: Box<SangerRejection>,
    },
}

/// Runs the shared read, alignment, and variant stages for one trace.
pub(crate) fn build(
    trace: &Chromatogram,
    reference: &Reference,
    config: &Config,
    profile: &Profile,
) -> Result<CompletedObservation> {
    let processed = read_processing::process(trace, &config.sanger)?;
    match observe(trace, processed, reference, config, profile)? {
        Observation::Placed(completed) => Ok(*completed),
        Observation::Unplaced { error, .. } => Err(error),
    }
}

/// Runs the alignment and variant stages for one processed read; a read the
/// core cannot place is an [`Observation::Unplaced`] rather than an error.
pub(crate) fn observe(
    trace: &Chromatogram,
    processed: ProcessedRead,
    reference: &Reference,
    config: &Config,
    profile: &Profile,
) -> Result<Observation> {
    let stage = tracing::info_span!("read_evidence").entered();
    let evidence = read_processing::read_evidence(&processed, &config.sanger.sanger_evidence)?;
    record_evidence(&evidence, config);
    drop(stage);
    let ProcessedRead {
        calls,
        signal,
        callability,
        quality,
        second,
        mixtures: _,
        tracts: _,
        warnings: read_warnings,
    } = processed;
    let called = match place(
        ReadIdentity {
            input_name: trace.source_name.clone(),
            input_sha256: trace.source_sha256.clone(),
            declared: None,
        },
        evidence,
        None,
        reference,
        config,
        profile,
    )? {
        Called::Placed(called) => *called,
        Called::Unplaced { rejection, error } => {
            return Ok(Observation::Unplaced {
                rejection,
                error,
                sanger: Box::new(SangerRejection {
                    integrity: signal.integrity,
                    callability,
                }),
            });
        }
    };

    let warning_total = warning_total(&read_warnings, &called);

    Ok(Observation::Placed(Box::new(CompletedObservation {
        read: ReadObservation {
            called,
            sanger: SangerAttachment {
                calls,
                signal,
                callability,
                quality,
                second,
            },
        },
        warning_total,
    })))
}

/// Records what one read's evidence holds, as the core receives it.
pub(crate) fn record_evidence(evidence: &ReadEvidence, config: &Config) {
    tracing::info!(
        event = "read_evidence_completed",
        calls = evidence.calls().len(),
        informative = ?evidence.informative(),
        masked_calls = evidence
            .calls()
            .iter()
            .filter(|call| call.placement != Placement::Scored)
            .count(),
        vetoed_calls = evidence
            .calls()
            .iter()
            .filter(|call| {
                evidence
                    .reason_order()
                    .iter()
                    .any(|reason| call.identity.withheld.contains(*reason))
            })
            .count(),
        minimum_peak_height = config.sanger.sanger_evidence.minimum_peak_height,
        relative_quality_threshold = config.sanger.sanger_evidence.relative_quality_threshold,
    );
}

/// What the core made of one read's evidence.
pub(crate) enum Called {
    /// The core placed and called the read.
    Placed(Box<CalledRead>),
    /// The core could not place the read (SRS-SAMPLE-029).
    Unplaced {
        rejection: PlacementRejection,
        error: Error,
    },
}

/// Places and calls one read's evidence with the core, as [`call`] does; a
/// read the core cannot place is [`Called::Unplaced`] rather than an error.
pub(crate) fn place(
    identity: ReadIdentity,
    evidence: ReadEvidence,
    placed: Option<Alignment>,
    reference: &Reference,
    config: &Config,
    profile: &Profile,
) -> Result<Called> {
    match call(identity, evidence, placed, reference, config, profile) {
        Ok(called) => Ok(Called::Placed(Box::new(called))),
        Err(error) => match PlacementRejection::of(&error) {
            Some(rejection) => Ok(Called::Unplaced {
                rejection,
                error: error.into(),
            }),
            None => Err(error.into()),
        },
    }
}

/// Places and calls one read's evidence with the core: evidence from a
/// trace's Sanger products, or from a recorded read evidence document
/// (PROP-0011), with its recorded placement when a `dna.placed_read/v1`
/// document gives one.
pub(crate) fn call(
    identity: ReadIdentity,
    evidence: ReadEvidence,
    placed: Option<Alignment>,
    reference: &Reference,
    config: &Config,
    profile: &Profile,
) -> dna_kernel::error::Result<CalledRead> {
    read_call::call_read_placed(
        identity,
        evidence,
        placed,
        reference,
        &CoreRun {
            config: &config.core,
            configuration_sha256: &config.source_sha256,
            regions: &profile.regions,
        },
    )
}

/// The read's operational warning total, logged when there are any.
pub(crate) fn warning_total(read_warnings: &ReadWarnings, called: &CalledRead) -> usize {
    let excluded_variant_candidates = called.variants.excluded_count();
    let reference_origin_wrap = called.alignment.wraps_origin();
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
    warning_total
}
