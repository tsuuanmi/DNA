//! Serializable `dna.read_evidence/v1` contract: one read's evidence as the
//! modality hands it to the core, or why the modality set the read aside
//! (PROP-0011).

use std::collections::BTreeMap;

use serde::Serialize;

use crate::model::basecall_result::BasecallWarningSummaryResult;
use crate::model::result::{
    CallabilityResult, IntervalResult, PeakHeightsResult, PlacementResult, PluginResult,
    SignalQualityResult, TraceIntegrityResult,
};
use crate::model::sample_result::ReadRejectionResult;

/// One read evidence document: a read the modality hands to the core, or one
/// it set aside.
#[derive(Debug, Serialize)]
#[serde(untagged)]
pub(crate) enum EvidenceDocument {
    Read(Box<EvidenceResult>),
    Rejected(Box<RejectedEvidenceResult>),
}

impl EvidenceDocument {
    pub(crate) const fn schema_version(&self) -> &'static str {
        match self {
            Self::Read(document) => document.schema_version,
            Self::Rejected(document) => document.schema_version,
        }
    }
}

/// The evidence of a read the modality set aside: its cause and, when it
/// reached callability, its integrity and callability, as the sample
/// documents report a rejected read.
#[derive(Debug, Serialize)]
pub(crate) struct RejectedEvidenceResult {
    pub(crate) schema_version: &'static str,
    pub(crate) provenance: EvidenceProvenanceResult,
    pub(crate) rejection: ReadRejectionResult,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) integrity: Option<TraceIntegrityResult>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) callability: Option<CallabilityResult>,
}

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
    pub(crate) input: EvidenceInputResult,
    /// The configuration file the evidence was written under.
    pub(crate) configuration_sha256: String,
    pub(crate) plugins: Vec<PluginResult>,
    /// Each plugin's configuration digest, by plugin id: what the evidence
    /// depends on of the configuration.
    pub(crate) plugin_configuration_sha256: BTreeMap<&'static str, String>,
}

/// The trace the evidence was written from: its file name, from which the
/// sample documents name the read, and its content identity.
#[derive(Debug, Serialize)]
pub(crate) struct EvidenceInputResult {
    pub(crate) name: String,
    pub(crate) sha256: String,
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
    /// The spans of calls whose differences the modality asserts are not to
    /// stand unchecked, with its reason (SRS-CALL-024).
    pub(crate) unreliable: Vec<UnreliableEvidenceResult>,
    /// The homopolymer tracts the modality decodes from the read's own
    /// signal, with their length marginals (PROP-0012).
    pub(crate) tracts: Vec<TractEvidenceResult>,
}

/// One tract decode as the modality asserts it, over call indices.
#[derive(Debug, Serialize)]
pub(crate) struct TractEvidenceResult {
    /// The call where the read leaves phase, after the tract.
    pub(crate) onset: usize,
    /// The calls the decode fits, from the tract's first call.
    pub(crate) calls: IntervalResult,
    /// Calls from the tract's first one that read the run or the mark.
    pub(crate) extent: usize,
    pub(crate) run: char,
    pub(crate) mark: Option<char>,
    /// Every structure fitted, single run first, then A, C, G, T.
    pub(crate) candidates: Vec<StructureFitResult>,
    pub(crate) residual: f64,
    pub(crate) sigma: f64,
    pub(crate) autocorrelation: f64,
    pub(crate) first: LengthMarginalResult,
    pub(crate) total: LengthMarginalResult,
    pub(crate) cleaned: Option<CleanedLengthsResult>,
}

/// One structure the decode fitted.
#[derive(Debug, Serialize)]
pub(crate) struct StructureFitResult {
    pub(crate) mark: Option<char>,
    pub(crate) residual: f64,
}

/// One length marginal: lengths ascending with their shares and standard
/// errors, the dominant and next lengths, and the margin between them.
#[derive(Debug, Serialize)]
pub(crate) struct LengthMarginalResult {
    pub(crate) lengths: Vec<usize>,
    pub(crate) shares: Vec<f64>,
    pub(crate) errors: Vec<f64>,
    pub(crate) dominant: Option<usize>,
    pub(crate) next: Option<usize>,
    pub(crate) z: Option<f64>,
    pub(crate) z_adjusted: Option<f64>,
    pub(crate) margin: &'static str,
}

/// The refit with the read's own signatures and leakage.
#[derive(Debug, Serialize)]
pub(crate) struct CleanedLengthsResult {
    pub(crate) leakage: Option<Vec<LeakageShareResult>>,
    pub(crate) residual: f64,
    pub(crate) first: LengthMarginalResult,
}

/// One leakage offset's share.
#[derive(Debug, Serialize)]
pub(crate) struct LeakageShareResult {
    pub(crate) offset: i32,
    pub(crate) share: f64,
}

/// One unreliable span as the modality asserts it, over call indices.
#[derive(Debug, Serialize)]
pub(crate) struct UnreliableEvidenceResult {
    /// The span, 0-based half-open call indices.
    pub(crate) calls: IntervalResult,
    pub(crate) reason: &'static str,
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

/// The Sanger signal of the call: its trace position, its selected channel
/// peaks, its channel heights at the primary event, its relative quality,
/// and the cause of its second peak (SRS-CALL-021).
#[derive(Debug, Serialize)]
pub(crate) struct CallSignalResult {
    pub(crate) scan: usize,
    /// The selected A/C/G/T peak heights in the call's window, on the read
    /// strand.
    pub(crate) channels: PeakHeightsResult,
    /// A/C/G/T heights at the primary event on the read strand, when the call
    /// has one.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) peaks: Option<PeakHeightsResult>,
    pub(crate) relative_quality: u8,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) cause: Option<&'static str>,
}
