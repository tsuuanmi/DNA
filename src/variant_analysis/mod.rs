//! Public Variant Analysis capability.
//!
//! This module exposes reusable, typed DNA analysis independent of CLI logging
//! and JSON publication. Source-specific adapters normalize their output into
//! the same canonical result contract.

pub(crate) mod observation;

use std::path::Path;

use crate::error::Result;
use crate::input::sanger;
use crate::model::variant as internal_variant;

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
    /// Reference intervals covered by the selected alignment.
    pub reference_segments: Vec<ReferenceSegment>,
    /// Normalized reportable primary-sequence differences.
    pub variants: Vec<Variant>,
}

/// Stable reference identity carried by analysis results.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReferenceIdentity {
    /// FASTA record identifier of the reference.
    pub name: String,
    /// SHA-256 of the normalized reference sequence.
    pub sha256: String,
}

/// Zero-based, half-open covered interval on the reference.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReferenceSegment {
    /// First covered reference index (0-based, inclusive).
    pub start_0based: usize,
    /// End of the covered interval (0-based, exclusive).
    pub end_0based_exclusive: usize,
}

/// Supported called-variant type.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[non_exhaustive]
pub enum VariantKind {
    /// Single-nucleotide substitution.
    Snv,
    /// Anchored insertion of one or more bases.
    Ins,
    /// Anchored deletion of one or more bases.
    Del,
}

/// One reportable evidence-backed primary-sequence difference.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Variant {
    /// Reference record identifier the variant is placed on.
    pub contig: String,
    /// 1-based position of the first reference allele base.
    pub position_1based: usize,
    /// Reference allele, including the anchor base for indels.
    pub reference: String,
    /// Alternate allele, including the anchor base for indels.
    pub alternate: String,
    /// Variant type.
    pub kind: VariantKind,
}

/// Cross-modality boundary for evidence-backed variants against one reference.
///
/// This type does not imply right/left alignment, nomenclature, VCF
/// normalization, or another representation policy.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CalledVariantSet {
    /// Reference every variant in the set is placed on.
    pub reference: ReferenceIdentity,
    /// Evidence-backed called variants.
    pub variants: Vec<Variant>,
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
    let completed = observation::build(&inputs.trace, &inputs.reference, &inputs.config)?;

    let read = completed.read;
    let reference_segments = read
        .alignment
        .reference_segments
        .iter()
        .map(|segment| ReferenceSegment {
            start_0based: segment.start_0based,
            end_0based_exclusive: segment.end_0based_exclusive,
        })
        .collect();
    let variants = read.variants.reported.iter().map(project_variant).collect();

    Ok(VariantAnalysisResult {
        input_sha256: read.input_sha256,
        reference: reference_identity,
        configuration_sha256: read.configuration_sha256,
        reference_segments,
        variants,
    })
}

/// Projects an internal called variant into the public variant boundary.
pub(crate) fn project_variant(variant: &internal_variant::Variant) -> Variant {
    Variant {
        contig: variant.contig.clone(),
        position_1based: variant.position_1based,
        reference: variant.reference.clone(),
        alternate: variant.alternate.clone(),
        kind: match variant.kind {
            internal_variant::VariantKind::Snv => VariantKind::Snv,
            internal_variant::VariantKind::Ins => VariantKind::Ins,
            internal_variant::VariantKind::Del => VariantKind::Del,
        },
    }
}
