//! Serializable `dna.review/v1` contract: one sample's review list.

use serde::Serialize;

use crate::model::notation_result::EmpopResult;
use crate::model::result::{IntervalResult, PluginResult, ProfileResult, ReferenceResult};
use crate::model::sample_result::{ConventionResult, ReadRejectionResult};
use dna_align::model::alignment::Orientation;
use dna_sample::model::sample_evidence::LocusState;

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
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) empop: Option<EmpopResult>,
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

/// One reference column: a position and an insertion ordinal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize)]
pub(crate) struct ColumnResult {
    pub(crate) position: usize,
    pub(crate) ordinal: usize,
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
    PossibleHeteroplasmy {
        position: usize,
        base: char,
        reads: Vec<ShareResult>,
    },
    EstimatedLength {
        runs: Vec<RunResult>,
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
    ReadEdge {
        read: String,
    },
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

/// A read showing the second base, with its share.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub(crate) struct ShareResult {
    pub(crate) read: String,
    pub(crate) strand: Orientation,
    pub(crate) share: f64,
}

/// A run whose length is an estimate, as `dna.consensus/v1` lists it.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub(crate) struct RunResult {
    pub(crate) base: char,
    pub(crate) length: usize,
    pub(crate) length_evidence: &'static str,
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

/// How far the core trusts a call (`dna::sample_analysis::TrustScope`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum ScopeResult {
    /// The base, a run's length and ends, and the indels it shows.
    Full,
    /// The base only.
    Base,
}

/// The strongest profile channel other than the call.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub(crate) struct SecondResult {
    pub(crate) base: char,
    pub(crate) share: f64,
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
    /// The call succeeded and no item is a `check`.
    pub(crate) clear: bool,
}

/// The profile's reportable positions, and those inside `check` items: what a
/// reviewer looks at, against what DNA decided alone.
#[derive(Debug, Default, PartialEq, Eq, Serialize)]
pub(crate) struct PositionCounts {
    pub(crate) reportable: usize,
    pub(crate) check: usize,
}

/// Items carrying each label, in label priority order.
#[derive(Debug, Default, PartialEq, Eq, Serialize)]
pub(crate) struct LabelCounts {
    pub(crate) strands_disagree: usize,
    pub(crate) read_disagrees: usize,
    pub(crate) possible_heteroplasmy: usize,
    pub(crate) estimated_length: usize,
    pub(crate) undecided: usize,
    pub(crate) not_called: usize,
    pub(crate) one_read: usize,
    pub(crate) one_strand: usize,
    pub(crate) noisy_signal: usize,
    pub(crate) not_covered: usize,
    pub(crate) phase_shift: usize,
    pub(crate) read_edge: usize,
}
