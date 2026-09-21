//! Compact serializable `dna.analysis/v7` contract.

use serde::Serialize;

use crate::model::alignment::Orientation;
use crate::model::reference::ReferenceTopology;
use crate::model::variant::{VariantCallRole, VariantKind};

/// Successful compact analysis document.
#[derive(Debug, Serialize)]
pub struct AnalysisResult {
    pub(crate) schema_version: &'static str,
    pub(crate) provenance: ProvenanceResult,
    pub(crate) read: ReadResult,
    pub(crate) signal_quality: DNAQualityResult,
    pub(crate) alignment: AlignmentResult,
    pub(crate) variants: Vec<VariantResult>,
    pub(crate) warnings: WarningSummaryResult,
}

/// Deterministic input identities retained for an analysis.
#[derive(Debug, Serialize)]
pub struct ProvenanceResult {
    pub(crate) input: InputResult,
    pub(crate) reference: ReferenceResult,
    pub(crate) configuration_sha256: String,
}

/// Input trace identity without an identifying filename or decoded bulk data.
#[derive(Debug, Serialize)]
pub struct InputResult {
    pub(crate) sha256: String,
}

/// Reference identity used by the selected alignment.
#[derive(Debug, Serialize)]
pub struct ReferenceResult {
    pub(crate) name: String,
    pub(crate) topology: ReferenceTopology,
    pub(crate) sha256: String,
}

/// Call count and retained interval without complete sequence strings.
#[derive(Debug, Serialize)]
pub struct ReadResult {
    pub(crate) call_count: usize,
    pub(crate) trim: IntervalResult,
}

/// A shared 0-based half-open result interval.
#[derive(Debug, Serialize)]
pub struct IntervalResult {
    pub(crate) start: usize,
    pub(crate) end: usize,
}

/// Compact observation-only integrity evidence shared by read result contracts.
#[derive(Debug, Clone, Serialize)]
pub struct TraceIntegrityResult {
    pub(crate) ploc_count: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) vendor_primary_count: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) vendor_quality_count: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) minimum_ploc_spacing: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) median_ploc_spacing: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) maximum_ploc_spacing: Option<usize>,
    pub(crate) clipped_channel_samples: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) maximum_to_median_event_signal_ratio: Option<f64>,
}

/// Shared merged observation-only signal-quality and trace-integrity evidence.
#[derive(Debug, Serialize)]
pub struct DNAQualityResult {
    pub(crate) integrity: TraceIntegrityResult,
    pub(crate) noisy_regions: Vec<NoisyRegionResult>,
}

/// A union of overlapping or adjacent candidate-noisy windows.
#[derive(Debug, Serialize)]
pub struct NoisyRegionResult {
    pub(crate) calls: IntervalResult,
    pub(crate) samples: IntervalResult,
    pub(crate) minimum_primary_snr: f64,
}

/// Concise summary of the selected alignment.
#[derive(Debug, Serialize)]
pub struct AlignmentResult {
    pub(crate) orientation: Orientation,
    pub(crate) callable_bases: usize,
    pub(crate) identity: f64,
    pub(crate) unresolved_bases: usize,
    pub(crate) gap_opens: usize,
    pub(crate) reference_segments: Vec<IntervalResult>,
    pub(crate) wraps_origin: bool,
}

/// Compact normalized variant with mapped trace calls.
#[derive(Debug, Serialize)]
pub struct VariantResult {
    pub(crate) position: usize,
    pub(crate) reference: String,
    pub(crate) alternate: String,
    pub(crate) kind: VariantKind,
    pub(crate) calls: Vec<VariantCallResult>,
}

/// Co-located reference-oriented A/C/G/T channel heights.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct PeakHeightsResult {
    #[serde(rename = "A")]
    pub(crate) a: i32,
    #[serde(rename = "C")]
    pub(crate) c: i32,
    #[serde(rename = "G")]
    pub(crate) g: i32,
    #[serde(rename = "T")]
    pub(crate) t: i32,
}

impl From<[i32; 4]> for PeakHeightsResult {
    fn from(value: [i32; 4]) -> Self {
        Self {
            a: value[0],
            c: value[1],
            g: value[2],
            t: value[3],
        }
    }
}

/// Reviewer-facing signal evidence for one variant-associated call.
#[derive(Debug, Serialize)]
pub struct VariantCallResult {
    pub(crate) role: VariantCallRole,
    pub(crate) base: char,
    pub(crate) peaks: PeakHeightsResult,
    /// Uncalibrated relative score exposed under the concise public name.
    pub(crate) quality: u8,
}

/// Public non-fatal analysis counts.
#[derive(Debug, Serialize)]
pub struct WarningSummaryResult {
    pub(crate) unresolved_primary_calls: usize,
    pub(crate) multi_channel_unresolved_calls: usize,
    pub(crate) ploc_vendor_length_mismatches: usize,
    pub(crate) clipped_channel_samples: usize,
    pub(crate) excluded_variant_candidates: usize,
}
