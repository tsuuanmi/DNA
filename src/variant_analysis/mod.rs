//! Public Variant Analysis capability.
//!
//! This module exposes reusable, typed DNA analysis independent of CLI logging
//! and JSON publication. Source-specific adapters normalize their output into
//! the same canonical result contract.

pub(crate) mod observation;
pub(crate) mod read_call;

use std::path::Path;

use crate::error::Result;
use crate::input::sanger;
use crate::profile::ProfileIdentity;
use crate::variant::{CalledVariantSet, ReferenceIdentity, Variant};

/// Typed result of one reference-guided variant analysis.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct VariantAnalysisResult {
    /// SHA-256 identity of the analyzed source artifact.
    pub input_sha256: String,
    /// Identity of the reference sequence used for alignment and normalization.
    pub reference: ReferenceIdentity,
    /// SHA-256 identity of the validated scientific configuration.
    pub configuration_sha256: String,
    /// Identity of the target profile the configuration references.
    pub profile: ProfileIdentity,
    /// Reference intervals covered by the selected alignment.
    pub reference_segments: Vec<ReferenceSegment>,
    /// Normalized reportable primary-sequence differences.
    pub variants: Vec<Variant>,
}

/// Zero-based, half-open covered interval on the reference.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReferenceSegment {
    /// First covered reference index (0-based, inclusive).
    pub start_0based: usize,
    /// End of the covered interval (0-based, exclusive).
    pub end_0based_exclusive: usize,
}

impl VariantAnalysisResult {
    /// Projects this analysis result into the common called-variant boundary.
    #[must_use]
    pub fn called_variants(&self) -> CalledVariantSet {
        CalledVariantSet {
            reference: self.reference.clone(),
            variants: self.variants.clone(),
        }
    }
}

/// Runs the current Sanger AB1 adapter through the canonical Variant Analysis
/// capability without CLI logging or JSON publication side effects.
///
/// # Errors
///
/// Returns [`Error`](crate::error::Error) when an input or configuration file is
/// unreadable or invalid, or when a scientific stage cannot produce a uniquely
/// interpretable result (for example an unaligned or low-identity read).
pub fn analyze_sanger(
    trace: &Path,
    reference: &Path,
    config: &Path,
) -> Result<VariantAnalysisResult> {
    let prepared = sanger::prepare_analysis(trace, reference, config)?;
    let inputs = sanger::load_analysis(prepared)?;
    let reference_identity = ReferenceIdentity {
        name: inputs.reference.name.clone(),
        sha256: inputs.reference.sequence_sha256.clone(),
    };
    let completed = observation::build(
        &inputs.trace,
        &inputs.reference,
        &inputs.config,
        &inputs.profile,
    )?;

    let read = completed.read.called;
    let reference_segments = read
        .alignment
        .reference_segments
        .iter()
        .map(|segment| ReferenceSegment {
            start_0based: segment.start_0based,
            end_0based_exclusive: segment.end_0based_exclusive,
        })
        .collect();
    let variants = read.variants.reported.iter().map(Variant::from).collect();

    Ok(VariantAnalysisResult {
        input_sha256: read.input_sha256,
        reference: reference_identity,
        configuration_sha256: read.configuration_sha256,
        profile: inputs.profile.identity().clone(),
        reference_segments,
        variants,
    })
}
