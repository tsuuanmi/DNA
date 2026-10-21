//! Serializable `dna.review/v1` contract: one sample's review list.

use serde::Serialize;

use crate::model::consensus_result::StretchResult;

use crate::model::notation_result::EmpopResult;
use crate::model::result::{
    ColumnResult, IntervalResult, PluginResult, ProfileResult, ReferenceResult, ScopeResult,
};
use crate::model::sample_result::{ConventionResult, ReadRejectionResult};
use dna_align::model::alignment::Orientation;
use dna_sample::model::sample_evidence::LocusState;
use dna_sanger::model::callability::PhaseState;

/// One sample's review document.
#[derive(Debug, Serialize)]
pub(crate) struct ReviewResult {
    pub(crate) schema_version: &'static str,
    pub(crate) sample_id: String,
    pub(crate) status: ReviewStatus,
    /// Why the call failed; only for `call_failed`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) error: Option<String>,
    pub(crate) provenance: ReviewProvenanceResult,
    pub(crate) reads: Vec<ReviewReadResult>,
    /// The sample's calls; only for `called` under a profile with notation.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) calls: Option<Vec<ReviewCallResult>>,
    /// Merged callable ranges, 1-based inclusive; only for `called`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) callable: Option<Vec<IntervalResult>>,
    /// Reportable positions (1-based, increasing) the consensus decided as
    /// the reference by reference confirmations alone (ADR-0080); only for
    /// `called`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) reference_only: Option<Vec<usize>>,
    /// The stretches at the reportable regions' edges the call does not
    /// analyse, with their reason (SRS-CONS-012): no item and no zone holds
    /// them; only for `called`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) not_analysed: Option<Vec<StretchResult>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) empop: Option<EmpopResult>,
    /// The reportable positions as zones a reviewer opens whole, in
    /// reference order; only for `called`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) zones: Option<Vec<ZoneResult>>,
    pub(crate) items: Vec<ItemResult>,
    pub(crate) summary: SummaryResult,
}

/// How far the sample got.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum ReviewStatus {
    /// The analysis and the call succeeded.
    Called,
    /// Every read was rejected.
    Rejected,
    /// The analysis succeeded but the call failed.
    CallFailed,
}

/// Identities of the inputs, plugins, and labelling rules.
#[derive(Debug, Serialize)]
pub(crate) struct ReviewProvenanceResult {
    pub(crate) reference: ReferenceResult,
    pub(crate) configuration_sha256: String,
    pub(crate) profile: ProfileResult,
    pub(crate) plugins: Vec<PluginResult>,
    /// Version of the labelling rules, the `review` plugin version.
    pub(crate) rules: u32,
    /// The run signal the calls were judged against, when one was given.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) signal: Option<SignalProvenanceResult>,
}

/// The identity of a run signal a review read.
#[derive(Debug, Serialize)]
pub(crate) struct SignalProvenanceResult {
    pub(crate) run_id: String,
    pub(crate) sha256: String,
    pub(crate) selection: u32,
}

/// One read: admitted with its strand, or rejected with its cause.
#[derive(Debug, Serialize)]
pub(crate) struct ReviewReadResult {
    pub(crate) name: String,
    pub(crate) sha256: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) strand: Option<Orientation>,
    pub(crate) admitted: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) rejection: Option<ReadRejectionResult>,
    /// An admitted read's masked stretches over the reportable regions.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) masked: Option<Vec<MaskedSpanResult>>,
    /// An admitted read's unreliable stretches over the reportable regions,
    /// with the modality's reason (SRS-CALL-024).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) unreliable: Option<Vec<UnreliableSpanResult>>,
    /// An admitted read's phase mixtures, each component compared with the
    /// reference on its own.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) mixtures: Option<Vec<MixtureResult>>,
}

/// A stretch after a repeat run whose signal the modality asserts is the
/// read's own sequence superposed with copies offset by whole calls
/// (ADR-0082): its components, and the reference positions of the calls the
/// components do not explain.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub(crate) struct MixtureResult {
    /// The components explain the stretch at the read's own noise level;
    /// otherwise they are listed by offset and fraction alone.
    pub(crate) explained: bool,
    pub(crate) components: Vec<ComponentResult>,
    pub(crate) unexplained: Vec<usize>,
}

/// One sequence of a phase mixture, compared with the reference on its own:
/// the calls by which it lags the dominant sequence (`0` for the dominant
/// one), its share of the fitted signal, the reference interval it covers,
/// its bases there on the reference strand, and its calls in the profile
/// notation.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub(crate) struct ComponentResult {
    pub(crate) offset: i32,
    pub(crate) fraction: f64,
    #[serde(flatten, skip_serializing_if = "Option::is_none")]
    pub(crate) compared: Option<ComparedResult>,
}

/// A component of an explained mixture, compared with the reference.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub(crate) struct ComparedResult {
    pub(crate) first: usize,
    pub(crate) last: usize,
    pub(crate) sequence: String,
    pub(crate) calls: Vec<String>,
}

/// A stretch of reference positions an admitted read's calls cover while its
/// signal is masked, from the modality's callability segments: the read does
/// not show its bases there.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub(crate) struct MaskedSpanResult {
    pub(crate) first: usize,
    pub(crate) last: usize,
    pub(crate) state: PhaseState,
    pub(crate) after_repeat: bool,
}

/// A maximal stretch of reportable positions that a read's calls cover inside
/// one of its unreliable spans.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub(crate) struct UnreliableSpanResult {
    pub(crate) first: usize,
    pub(crate) last: usize,
    pub(crate) reason: &'static str,
}

/// A maximal stretch of reportable positions within one region of one zone
/// kind (and, for a repeat, one window), with the `check` items touching it.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub(crate) struct ZoneResult {
    pub(crate) kind: ZoneKind,
    pub(crate) first: usize,
    pub(crate) last: usize,
    /// The profile's name of the repeat window; only for `repeat`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) name: Option<String>,
    /// The `check` items touching the zone.
    pub(crate) checks: usize,
}

/// What the sample's reads show over a position, in precedence order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum ZoneKind {
    /// A repeat window of the profile's nomenclature.
    Repeat,
    /// No read gives a callable call there.
    NotCallable,
    /// A read reads past a repeat out of phase there.
    Dephased,
    /// A read's signal is masked there otherwise.
    Masked,
    /// Every read covering it is in phase or vouched for.
    Clean,
}

/// One call of the sample in the profile notation.
#[derive(Debug, Serialize)]
pub(crate) struct ReviewCallResult {
    pub(crate) call: String,
    pub(crate) position: usize,
    pub(crate) ordinal: usize,
    pub(crate) reference: Option<char>,
    pub(crate) alternate: Option<char>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) convention: Option<ConventionResult>,
}

/// One position or stretch a reviewer should look at.
#[derive(Debug, Serialize)]
pub(crate) struct ItemResult {
    pub(crate) first: ColumnResult,
    pub(crate) last: ColumnResult,
    pub(crate) kind: ItemKind,
    pub(crate) decision: Decision,
    /// The decided sequence over the item; `None` when undecided or not
    /// callable.
    pub(crate) sequence: Option<String>,
    pub(crate) calls: Vec<String>,
    pub(crate) labels: Vec<LabelResult>,
    pub(crate) evidence: Vec<EvidenceReadResult>,
}

/// Whether DNA calls a difference in the item.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum ItemKind {
    /// The item holds a call.
    Variant,
    /// DNA could not assure the item.
    Check,
    /// No call and no label that leaves it unassured, listed for a context
    /// label a reviewer may still want: an explained possible heteroplasmy.
    Note,
}

/// What DNA decided over the item.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum Decision {
    /// Decided, and every position is callable.
    Called,
    /// Some positions are callable.
    PartlyCallable,
    /// No position is callable.
    NotCallable,
}

/// One evidence label of an item.
#[derive(Debug, Clone, PartialEq, Serialize)]
#[serde(tag = "label", rename_all = "snake_case")]
pub(crate) enum LabelResult {
    StrandsDisagree {
        reads: Vec<SideResult>,
    },
    ReadDisagrees {
        read: String,
        strand: Orientation,
        sequence: String,
    },
    NotCalled,
    LongRun {
        runs: Vec<GrowthResult>,
    },
    EdgeDeletion {
        read: String,
    },
    RunDeletion {
        read: String,
    },
    PossibleHeteroplasmy {
        position: usize,
        base: char,
        /// A strand's qualifying reads all show a second peak whose shape
        /// explains it: the label is context, not a check.
        #[serde(skip_serializing_if = "std::ops::Not::not")]
        explained: bool,
        reads: Vec<ShareResult>,
    },
    /// A difference resting on one read whose call at `position` shows
    /// another base than the reference with the reference base as a second
    /// peak above its threshold (SRS-CALL-020): a check.
    ReferenceSecond {
        position: usize,
        read: String,
        strand: Orientation,
        share: f64,
        /// The threshold the share exceeded (SRS-CALL-020).
        threshold: f64,
        /// The key's calls in the run signal, the sample's own taken out.
        #[serde(skip_serializing_if = "Option::is_none")]
        calls: Option<u64>,
    },
    EstimatedLength {
        runs: Vec<RunResult>,
    },
    /// A run of an anchored repeat window whose reads' decodes do not
    /// confirm DNA's length (SRS-API-020): a check.
    LengthUnconfirmed {
        runs: Vec<LengthRunResult>,
    },
    /// The reads' decodes of an anchored window whose run lengths the
    /// profile does not record: context, never a check.
    LengthNotRecorded {
        runs: Vec<LengthRunResult>,
    },
    Undecided,
    OneRead {
        read: String,
    },
    OneStrand {
        strand: Orientation,
        reads: Vec<String>,
    },
    NoisySignal {
        read: String,
        reason: &'static str,
    },
    NotCovered,
    PhaseShift {
        read: String,
        from: usize,
    },
    DephasedEnd {
        runs: Vec<RunResult>,
    },
    ReadEdge {
        read: String,
    },
    PoorSignal {
        read: String,
    },
    BesideNotAnalysed,
}

/// A read on one side of a site's disagreement.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub(crate) struct SideResult {
    pub(crate) read: String,
    pub(crate) strand: Orientation,
    pub(crate) role: Role,
}

/// Which side of the consensus call a read takes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum Role {
    Supporting,
    Opposing,
}

/// A read showing the second base, with its peak share (SRS-CALL-015).
#[derive(Debug, Clone, PartialEq, Serialize)]
pub(crate) struct ShareResult {
    pub(crate) read: String,
    pub(crate) strand: Orientation,
    pub(crate) share: f64,
    /// The threshold the share exceeded (SRS-CALL-020).
    pub(crate) threshold: f64,
    /// The key's calls in the run signal, the sample's own taken out.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) calls: Option<u64>,
}

/// A decided run longer than its reference run by `run_growth` or more,
/// whose length rests on an anchored end or a single read.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub(crate) struct GrowthResult {
    pub(crate) base: char,
    pub(crate) length: usize,
    pub(crate) reference_length: usize,
    /// How the consensus knows the length (`dna.consensus/v1`); `None` when
    /// the site was decided position by position.
    pub(crate) length_evidence: Option<&'static str>,
}

/// A run whose length is an estimate, as `dna.consensus/v1` lists it.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub(crate) struct RunResult {
    pub(crate) base: char,
    pub(crate) length: usize,
    pub(crate) length_evidence: &'static str,
}

/// One run of an anchored repeat window: which run (1 before the anchor, 2
/// after it, 1 alone when the anchor is lost), its base, DNA's length, how
/// the reads' decodes agree, and each read's length.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub(crate) struct LengthRunResult {
    pub(crate) run: usize,
    pub(crate) base: char,
    pub(crate) called: Option<usize>,
    pub(crate) agreement: &'static str,
    /// The lab's convention reading of a near-equal pair, a suggestion only:
    /// never DNA's call.
    pub(crate) convention: Option<LengthConventionResult>,
    pub(crate) reads: Vec<LengthReadResult>,
}

/// A run length a reporting convention would read, with the convention.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub(crate) struct LengthConventionResult {
    pub(crate) length: usize,
    pub(crate) rule: &'static str,
}

/// One read's length of a run: decoded from its first-run marginal
/// (`decode`), derived from its total when the run is the next it reads
/// (`derived`), or read in phase, both ends shown (`in_phase`).
#[derive(Debug, Clone, PartialEq, Serialize)]
pub(crate) struct LengthReadResult {
    pub(crate) read: String,
    pub(crate) strand: Orientation,
    pub(crate) length: usize,
    pub(crate) share: f64,
    pub(crate) next: Option<usize>,
    pub(crate) margin: &'static str,
    pub(crate) source: &'static str,
}

/// What one admitted read shows over an item.
#[derive(Debug, Serialize)]
pub(crate) struct EvidenceReadResult {
    pub(crate) read: String,
    pub(crate) sha256: String,
    pub(crate) strand: Orientation,
    pub(crate) positions: Vec<EvidencePositionResult>,
}

/// One read's observation of one reference position.
#[derive(Debug, Serialize)]
pub(crate) struct EvidencePositionResult {
    pub(crate) position: usize,
    pub(crate) ordinal: usize,
    pub(crate) base: char,
    pub(crate) state: LocusState,
    pub(crate) call: Option<usize>,
    /// How far the core trusts the call: `full`, or `base` for a call the
    /// modality rescued that shows its base only; `null` for a deletion.
    pub(crate) scope: Option<ScopeResult>,
    pub(crate) scan: Option<usize>,
    pub(crate) quality: Option<u8>,
    pub(crate) second: Option<SecondResult>,
}

/// The call's second base on the reference strand and its peak share
/// (SRS-CALL-015).
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub(crate) struct SecondResult {
    pub(crate) base: char,
    pub(crate) share: f64,
    /// Why the call shows a second peak, from its signal: the modality's
    /// cause, refined by the run signal when one is given.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) cause: Option<&'static str>,
}

/// Counts over the items.
#[derive(Debug, Serialize)]
pub(crate) struct SummaryResult {
    pub(crate) items: usize,
    pub(crate) checks: usize,
    pub(crate) variants: usize,
    pub(crate) labels: LabelCounts,
    /// Reportable positions, and those a `check` item holds.
    pub(crate) positions: PositionCounts,
    /// Zones, and those a `check` item touches.
    pub(crate) zones: ZoneCounts,
    /// The call succeeded and no item is a `check`.
    pub(crate) clear: bool,
}

/// The profile's reportable positions, and those inside `check` items: what a
/// reviewer looks at, against what DNA decided alone.
#[derive(Debug, Default, PartialEq, Eq, Serialize)]
pub(crate) struct PositionCounts {
    pub(crate) reportable: usize,
    pub(crate) check: usize,
    /// Reportable positions decided by reference confirmations alone.
    pub(crate) reference_only: usize,
    /// Reportable positions the call does not analyse (SRS-CONS-012).
    pub(crate) not_analysed: usize,
}

/// The sample's zones, and those a `check` item touches: what a reviewer
/// opens.
#[derive(Debug, Default, PartialEq, Eq, Serialize)]
pub(crate) struct ZoneCounts {
    pub(crate) total: usize,
    pub(crate) check: usize,
}

/// Items carrying each label, in label priority order.
#[derive(Debug, Default, PartialEq, Eq, Serialize)]
pub(crate) struct LabelCounts {
    pub(crate) strands_disagree: usize,
    pub(crate) read_disagrees: usize,
    pub(crate) possible_heteroplasmy: usize,
    pub(crate) reference_second: usize,
    pub(crate) estimated_length: usize,
    pub(crate) length_unconfirmed: usize,
    pub(crate) undecided: usize,
    pub(crate) not_called: usize,
    pub(crate) long_run: usize,
    pub(crate) edge_deletion: usize,
    pub(crate) run_deletion: usize,
    pub(crate) beside_not_analysed: usize,
    pub(crate) one_read: usize,
    pub(crate) one_strand: usize,
    pub(crate) noisy_signal: usize,
    pub(crate) not_covered: usize,
    pub(crate) phase_shift: usize,
    pub(crate) dephased_end: usize,
    pub(crate) read_edge: usize,
    pub(crate) poor_signal: usize,
    pub(crate) length_not_recorded: usize,
}
