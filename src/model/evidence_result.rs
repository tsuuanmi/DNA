//! Serializable `dna.read_evidence/v1` contract: one read's evidence as the
//! modality hands it to the core (PROP-0011).

use serde::Serialize;

use crate::model::basecall_result::BasecallWarningSummaryResult;
use crate::model::result::{
    CallabilityResult, InputResult, IntervalResult, PeakHeightsResult, PlacementResult,
    PluginResult, SignalQualityResult,
};

/// One read's evidence document.
#[derive(Debug, Serialize)]
pub(crate) struct EvidenceResult {
    pub(crate) schema_version: &'static str,
    pub(crate) provenance: EvidenceProvenanceResult,
    pub(crate) read: EvidenceReadResult,
    /// Trace integrity and merged noisy regions, as in `dna.basecalls/v3`.
    pub(crate) signal_quality: SignalQualityResult,
    /// The modality's warning counts, as in `dna.basecalls/v3`.
    pub(crate) warnings: BasecallWarningSummaryResult,
    pub(crate) calls: Vec<CallEvidenceResult>,
}

/// Deterministic input identities of the evidence.
#[derive(Debug, Serialize)]
pub(crate) struct EvidenceProvenanceResult {
    pub(crate) input: InputResult,
    pub(crate) configuration_sha256: String,
    pub(crate) plugins: Vec<PluginResult>,
}

/// What the evidence asserts about the read as a whole.
#[derive(Debug, Serialize)]
pub(crate) struct EvidenceReadResult {
    pub(crate) call_count: usize,
    /// The relative-quality trim interval, as in `dna.basecalls/v3`.
    pub(crate) trim: IntervalResult,
    /// The calls the modality retains, 0-based half-open call indices.
    pub(crate) informative: IntervalResult,
    /// Whether the modality vouches for the calls up to the read's ends.
    pub(crate) ends: Ends,
    /// The modality's reason labels in the order the core reports them.
    pub(crate) reasons: Vec<&'static str>,
    /// The modality's callability, as in `dna.basecalls/v3`.
    pub(crate) callability: CallabilityResult,
    /// The stretches the modality decomposes into the read's own sequence at
    /// whole-call offsets (ADR-0082).
    pub(crate) mixtures: Vec<MixtureEvidenceResult>,
}

/// One phase mixture as the modality asserts it, over call indices.
#[derive(Debug, Serialize)]
pub(crate) struct MixtureEvidenceResult {
    /// The stretch, 0-based half-open call indices, from the run whose length
    /// varies.
    pub(crate) calls: IntervalResult,
    pub(crate) run_calls: usize,
    pub(crate) explained: bool,
    pub(crate) residual: f64,
    pub(crate) components: Vec<MixtureComponentResult>,
    pub(crate) unexplained: Vec<usize>,
}

/// One component: its offset from the dominant sequence and its fraction.
#[derive(Debug, Serialize)]
pub(crate) struct MixtureComponentResult {
    pub(crate) offset: i32,
    pub(crate) fraction: f64,
}

/// Whether the read's ends bound its evidence.
#[derive(Debug, Clone, Copy, Serialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum Ends {
    Unvouched,
    Vouched,
}

/// One call's evidence, in call order: the kernel contract the core reads,
/// then the Sanger signal the modality judged it on.
#[derive(Debug, Serialize)]
pub(crate) struct CallEvidenceResult {
    pub(crate) index: usize,
    pub(crate) base: char,
    /// Basecall-independent A/C/G/T weights, when the call has signal.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) profile: Option<[f64; 4]>,
    pub(crate) placement: PlacementResult,
    pub(crate) identity: VerdictResult,
    pub(crate) indel: VerdictResult,
    pub(crate) bound: BoundResult,
    pub(crate) post_repeat: PostRepeatResult,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) second: Option<SecondBaseResult>,
    pub(crate) second_unexplained: bool,
    pub(crate) reference_only: bool,
    pub(crate) signal: CallSignalResult,
}

/// The modality's verdict on one use of the call.
#[derive(Debug, Serialize)]
pub(crate) struct VerdictResult {
    pub(crate) edge: bool,
    pub(crate) withheld: Vec<&'static str>,
}

/// How the call shows where a run of another base ends.
#[derive(Debug, Clone, Copy, Serialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum BoundResult {
    InPhase,
    Anchored,
    Withheld,
}

/// The call's part in a run the read bounds before a post-repeat phase loss.
#[derive(Debug, Serialize)]
pub(crate) struct PostRepeatResult {
    pub(crate) call: bool,
    pub(crate) bound: bool,
}

/// The second base of the call's signal, on the read strand.
#[derive(Debug, Serialize)]
pub(crate) struct SecondBaseResult {
    pub(crate) base: char,
    pub(crate) share: f64,
    pub(crate) stands_out: bool,
}

/// The Sanger signal of the call: its trace position, its strongest peak,
/// its channel heights at the primary event, its relative quality, and the
/// cause of its second peak (SRS-CALL-021).
#[derive(Debug, Serialize)]
pub(crate) struct CallSignalResult {
    pub(crate) scan: usize,
    pub(crate) height: u32,
    /// A/C/G/T heights at the primary event on the read strand, when the call
    /// has one.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) peaks: Option<PeakHeightsResult>,
    pub(crate) relative_quality: u8,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) cause: Option<&'static str>,
}
