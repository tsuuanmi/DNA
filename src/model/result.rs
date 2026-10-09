//! Compact serializable `dna.analysis/v9` contract.

use serde::Serialize;

use dna_core::model::alignment::Orientation;
use dna_core::model::variant::{VariantCallRole, VariantKind};
use dna_kernel::model::reference::ReferenceTopology;
use dna_sanger::model::callability::PhaseState;

/// Successful compact analysis document.
#[derive(Debug, Serialize)]
pub(crate) struct AnalysisResult {
    pub(crate) schema_version: &'static str,
    pub(crate) provenance: ProvenanceResult,
    pub(crate) read: ReadResult,
    pub(crate) signal_quality: SignalQualityResult,
    pub(crate) alignment: AlignmentResult,
    pub(crate) variants: Vec<VariantResult>,
    pub(crate) warnings: WarningSummaryResult,
}

/// Deterministic input identities retained for an analysis.
#[derive(Debug, Serialize)]
pub(crate) struct ProvenanceResult {
    pub(crate) input: InputResult,
    pub(crate) reference: ReferenceResult,
    pub(crate) configuration_sha256: String,
    pub(crate) profile: ProfileResult,
    pub(crate) plugins: Vec<PluginResult>,
}

/// Identity of a plugin that took part in producing a result (ADR-0069).
#[derive(Debug, Serialize)]
pub(crate) struct PluginResult {
    pub(crate) id: &'static str,
    pub(crate) family: &'static str,
    pub(crate) version: u32,
}

/// Identity of the target profile a result was produced under.
#[derive(Debug, Serialize)]
pub(crate) struct ProfileResult {
    pub(crate) id: String,
    pub(crate) sha256: String,
}

/// Input trace identity without an identifying filename or decoded bulk data.
#[derive(Debug, Serialize)]
pub(crate) struct InputResult {
    pub(crate) sha256: String,
}

/// Reference identity used by the selected alignment.
#[derive(Debug, Serialize)]
pub(crate) struct ReferenceResult {
    pub(crate) name: String,
    pub(crate) topology: ReferenceTopology,
    pub(crate) sha256: String,
}

/// Call count, retained interval, and callability without complete sequence strings.
#[derive(Debug, Serialize)]
pub(crate) struct ReadResult {
    pub(crate) call_count: usize,
    pub(crate) trim: IntervalResult,
    pub(crate) callability: CallabilityResult,
}

/// Signal-derived callability view of one read (`dna.read_callability/v1`).
#[derive(Debug, Serialize)]
pub(crate) struct CallabilityResult {
    pub(crate) callable_span: IntervalResult,
    pub(crate) segments: Vec<CallabilitySegmentResult>,
    pub(crate) masked_calls: usize,
}

/// One maximal run of calls sharing a phase state.
#[derive(Debug, Serialize)]
pub(crate) struct CallabilitySegmentResult {
    pub(crate) calls: IntervalResult,
    pub(crate) state: PhaseState,
    pub(crate) after_repeat: bool,
    /// Reported shadow offsets of a dephased segment, ascending.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) shadow_offsets: Option<Vec<i8>>,
}

/// A shared 0-based half-open result interval.
#[derive(Debug, Serialize)]
pub(crate) struct IntervalResult {
    pub(crate) start: usize,
    pub(crate) end: usize,
}

/// Compact observation-only integrity evidence shared by read result contracts.
#[derive(Debug, Clone, Serialize)]
pub(crate) struct TraceIntegrityResult {
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
pub(crate) struct SignalQualityResult {
    pub(crate) integrity: TraceIntegrityResult,
    pub(crate) noisy_regions: Vec<NoisyRegionResult>,
}

/// A union of overlapping or adjacent candidate-noisy windows.
#[derive(Debug, Serialize)]
pub(crate) struct NoisyRegionResult {
    pub(crate) calls: IntervalResult,
    pub(crate) samples: IntervalResult,
    pub(crate) minimum_primary_snr: f64,
}

/// Concise summary of the selected alignment.
#[derive(Debug, Serialize)]
pub(crate) struct AlignmentResult {
    pub(crate) orientation: Orientation,
    pub(crate) callable_bases: usize,
    pub(crate) identity: f64,
    pub(crate) unresolved_bases: usize,
    pub(crate) masked_bases: usize,
    pub(crate) gap_opens: usize,
    pub(crate) reference_segments: Vec<IntervalResult>,
    pub(crate) callable_reference_segments: Vec<IntervalResult>,
    pub(crate) wraps_origin: bool,
}

/// Compact normalized variant with mapped trace calls.
#[derive(Debug, Serialize)]
pub(crate) struct VariantResult {
    pub(crate) position: usize,
    pub(crate) reference: String,
    pub(crate) alternate: String,
    pub(crate) kind: VariantKind,
    pub(crate) calls: Vec<VariantCallResult>,
}

/// Co-located reference-oriented A/C/G/T channel heights.
#[derive(Debug, Clone, Copy, Serialize)]
pub(crate) struct PeakHeightsResult {
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
pub(crate) struct VariantCallResult {
    pub(crate) role: VariantCallRole,
    pub(crate) base: char,
    pub(crate) peaks: PeakHeightsResult,
    /// Uncalibrated relative score exposed under the concise public name.
    pub(crate) quality: u8,
}

/// Public non-fatal analysis counts.
#[derive(Debug, Serialize)]
pub(crate) struct WarningSummaryResult {
    pub(crate) unresolved_primary_calls: usize,
    pub(crate) multi_channel_unresolved_calls: usize,
    pub(crate) ploc_vendor_length_mismatches: usize,
    pub(crate) clipped_channel_samples: usize,
    pub(crate) excluded_variant_candidates: usize,
}
