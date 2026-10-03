//! Public Variant Analysis capability.
//!
//! This module exposes reusable, typed DNA analysis independent of CLI logging
//! and JSON publication. Source-specific adapters normalize their output into
//! the same canonical result contract.

pub(crate) mod observation;

use std::fmt;
use std::path::Path;

use crate::config::Config;
use crate::error::Result;
use crate::input::sanger;
use crate::logger::StageLog;
use crate::model::reference::Reference;
use crate::model::sanger::Chromatogram;
use crate::model::variant as internal_variant;

/// Typed result of one reference-guided variant analysis.
#[derive(Debug, Clone, PartialEq, Eq)]
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
    pub name: String,
    pub sha256: String,
}

/// Zero-based, half-open covered interval on the reference.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReferenceSegment {
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
}

/// Supported called-variant type.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum VariantKind {
    Snv,
    Ins,
    Del,
}

/// One reportable evidence-backed primary-sequence difference.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Variant {
    pub contig: String,
    pub position_1based: usize,
    pub reference: String,
    pub alternate: String,
    pub kind: VariantKind,
}

/// Cross-modality boundary for evidence-backed variants against one reference.
///
/// This type does not imply right/left alignment, nomenclature, VCF
/// normalization, or another representation policy.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CalledVariantSet {
    pub reference: ReferenceIdentity,
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

/// Immutable reference/configuration context reusable across independent Sanger traces.
///
/// This type owns no thread pool and performs no output publication. Callers may
/// share it across worker threads and choose their own outer parallelism.
#[derive(Debug, Clone)]
pub struct SangerAnalyzer {
    config: Config,
    reference: Reference,
    reference_identity: ReferenceIdentity,
}

impl SangerAnalyzer {
    /// Loads and validates one shared reference and scientific configuration.
    pub fn load(reference: &Path, config: &Path) -> Result<Self> {
        let context = sanger::load_analysis_context(reference, config)?;
        let reference_identity = reference_identity(&context.reference);
        Ok(Self {
            config: context.config,
            reference: context.reference,
            reference_identity,
        })
    }

    /// Analyzes one ABIF trace using the already loaded reference/configuration.
    pub fn analyze(&self, trace: &Path) -> Result<VariantAnalysisResult> {
        let trace = sanger::load_trace(trace)?;
        analyze_loaded(&trace, &self.reference, &self.config)
    }

    /// Returns the shared reference identity.
    #[must_use]
    pub const fn reference_identity(&self) -> &ReferenceIdentity {
        &self.reference_identity
    }

    /// Returns the shared validated configuration identity.
    #[must_use]
    pub fn configuration_sha256(&self) -> &str {
        &self.config.source_sha256
    }
}

/// Runs the current Sanger AB1 adapter through the canonical Variant Analysis
/// capability without CLI logging or JSON publication side effects.
pub fn analyze_sanger(
    trace: &Path,
    reference: &Path,
    config: &Path,
) -> Result<VariantAnalysisResult> {
    let prepared = sanger::prepare_analysis(trace, reference, config)?;
    let inputs = sanger::load_analysis(prepared)?;
    analyze_loaded(&inputs.trace, &inputs.reference, &inputs.config)
}

fn analyze_loaded(
    trace: &Chromatogram,
    reference: &Reference,
    config: &Config,
) -> Result<VariantAnalysisResult> {
    let identity = reference_identity(reference);
    let mut log = SilentStageLog;
    let mut stage = "read_processing";
    let completed = observation::build(trace, reference, config, &mut log, &mut stage)?;

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
        reference: identity,
        configuration_sha256: read.configuration_sha256,
        reference_segments,
        variants,
    })
}

fn reference_identity(reference: &Reference) -> ReferenceIdentity {
    ReferenceIdentity {
        name: reference.name.clone(),
        sha256: reference.sequence_sha256.clone(),
    }
}

fn project_variant(variant: &internal_variant::Variant) -> Variant {
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

struct SilentStageLog;

impl StageLog for SilentStageLog {
    fn info(&mut self, _module: &str, _line: u32, _message: fmt::Arguments<'_>) -> Result<()> {
        Ok(())
    }

    fn warn(&mut self, _module: &str, _line: u32, _message: fmt::Arguments<'_>) -> Result<()> {
        Ok(())
    }
}
