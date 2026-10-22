//! Serializable `dna.review/v1` contract: one sample's review list.

use serde::Serialize;

use crate::model::consensus_result::StretchResult;

use crate::model::notation_result::EmpopResult;
use crate::model::result::{
    ColumnResult, IntervalResult, PhaseStateResult, PluginResult, ProfileResult, ReferenceResult,
    ScopeResult,
};
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
    /// The run context the consensus calibrated its length rule by, when
    /// one was given (PROP-0014).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) run_context: Option<super::run_context_result::ContextUseResult>,
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
    pub(crate) state: PhaseStateResult,
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

/// The rank of a context label: a label that leaves an item assured, and an
/// `explained` possible heteroplasmy.
pub(crate) const CONTEXT_RANK: u8 = 4;

/// The kinds of evidence label. [`LabelKind::spec`] is the one table of their
/// published names, ranks, and whether they make an item a check (PROP-0004);
/// the schema's summary lists them in [`LabelKind::ALL`] order with the same
/// ranks and checks (`x-rank`, `x-check`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum LabelKind {
    StrandsDisagree,
    ReadDisagrees,
    SameStrandDissent,
    PossibleHeteroplasmy,
    ReferenceSecond,
    EstimatedLength,
    LengthUnconfirmed,
    Undecided,
    NotCalled,
    LongRun,
    EdgeDeletion,
    RunDeletion,
    RunIndel,
    OneStrandInsertion,
    BesideNotAnalysed,
    NotAnalysed,
    OneRead,
    OneStrand,
    NoisySignal,
    NotCovered,
    PhaseShift,
    DephasedEnd,
    ReadEdge,
    PoorSignal,
    LengthNotRecorded,
    TractLength,
}

impl LabelKind {
    /// Every kind, in the order the summary counts them.
    pub(crate) const ALL: [Self; 26] = [
        Self::StrandsDisagree,
        Self::ReadDisagrees,
        Self::SameStrandDissent,
        Self::PossibleHeteroplasmy,
        Self::ReferenceSecond,
        Self::EstimatedLength,
        Self::LengthUnconfirmed,
        Self::Undecided,
        Self::NotCalled,
        Self::LongRun,
        Self::EdgeDeletion,
        Self::RunDeletion,
        Self::RunIndel,
        Self::OneStrandInsertion,
        Self::BesideNotAnalysed,
        Self::NotAnalysed,
        Self::OneRead,
        Self::OneStrand,
        Self::NoisySignal,
        Self::NotCovered,
        Self::PhaseShift,
        Self::DephasedEnd,
        Self::ReadEdge,
        Self::PoorSignal,
        Self::LengthNotRecorded,
        Self::TractLength,
    ];

    /// The kind's published name, rank (the lower, the earlier a reviewer
    /// looks), and whether it leaves its item unassured, so that the item is
    /// a check: reads that disagree, a possible heteroplasmy, an estimated
    /// length, an undecided position, or no read. The read labels (one read,
    /// one strand, noisy signal, read edge, phase shift) are context for the
    /// reviewer.
    const fn spec(self) -> (&'static str, u8, bool) {
        match self {
            Self::StrandsDisagree => ("strands_disagree", 0, true),
            Self::ReadDisagrees => ("read_disagrees", 0, true),
            Self::SameStrandDissent => ("same_strand_dissent", 0, true),
            Self::PossibleHeteroplasmy => ("possible_heteroplasmy", 1, true),
            Self::ReferenceSecond => ("reference_second", 1, true),
            Self::EstimatedLength => ("estimated_length", 2, true),
            Self::LengthUnconfirmed => ("length_unconfirmed", 2, true),
            Self::Undecided => ("undecided", 2, true),
            Self::NotCalled => ("not_called", 2, true),
            Self::LongRun => ("long_run", 2, true),
            Self::EdgeDeletion => ("edge_deletion", 2, true),
            Self::RunDeletion => ("run_deletion", 2, true),
            Self::RunIndel => ("run_indel", 2, true),
            Self::OneStrandInsertion => ("one_strand_insertion", 2, true),
            Self::BesideNotAnalysed => ("beside_not_analysed", 2, true),
            Self::NotAnalysed => ("not_analysed", 2, true),
            Self::OneRead => ("one_read", 3, false),
            Self::OneStrand => ("one_strand", 3, false),
            Self::NoisySignal => ("noisy_signal", CONTEXT_RANK, false),
            Self::NotCovered => ("not_covered", CONTEXT_RANK, true),
            Self::PhaseShift => ("phase_shift", CONTEXT_RANK, false),
            Self::DephasedEnd => ("dephased_end", CONTEXT_RANK, false),
            Self::ReadEdge => ("read_edge", CONTEXT_RANK, false),
            Self::PoorSignal => ("poor_signal", CONTEXT_RANK, false),
            Self::LengthNotRecorded => ("length_not_recorded", CONTEXT_RANK, false),
            Self::TractLength => ("tract_length", CONTEXT_RANK, false),
        }
    }

    /// The published name.
    pub(crate) const fn name(self) -> &'static str {
        self.spec().0
    }

    /// The rank.
    pub(crate) const fn rank(self) -> u8 {
        self.spec().1
    }

    /// Whether a label of this kind makes its item a check.
    pub(crate) const fn check(self) -> bool {
        self.spec().2
    }
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
    /// A site whose supporting and opposing reads all lie on one strand
    /// (rules version 31, ADR-0090): clean reads of one strand disagree, and
    /// their count alone cannot settle it; a check. Lists the supporting,
    /// then the opposing reads.
    SameStrandDissent {
        reads: Vec<SideResult>,
    },
    NotCalled,
    LongRun {
        runs: Vec<GrowthResult>,
    },
    /// A decided deletion near a reportable region's end whose supporting
    /// reads all lie on one strand (rules version 33, ADR-0090): a check.
    /// Lists every supporting read.
    EdgeDeletion {
        strand: Orientation,
        reads: Vec<String>,
    },
    /// A decided site at least `review.run_growth` bases shorter than its
    /// reference, over a reference run, whose supporting reads all lie on
    /// one strand (rules version 33, ADR-0090): a check. Lists every
    /// supporting read.
    RunDeletion {
        strand: Orientation,
        reads: Vec<String>,
    },
    /// A decided site outside every repeat window that changes a homopolymer
    /// run's length and whose supporting reads all lie on one strand (rules
    /// versions 28 and 33, ADR-0090): a check. Lists every supporting read.
    RunIndel {
        strand: Orientation,
        reads: Vec<String>,
    },
    /// A decided insertion outside every repeat window and length tract
    /// whose supporting reads all lie on one strand, while no read of the
    /// other strand shows it (rules version 29): a check. `other` is the
    /// other strand's best view of the site (`not_covered`, `read_end`,
    /// `unresolved`, `reads_through`), and `from_start` the fewest calls
    /// from a supporting read's first call to its first call over the site.
    OneStrandInsertion {
        strand: Orientation,
        reads: Vec<String>,
        other: &'static str,
        from_start: usize,
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
    /// A difference resting on one strand, at a supporting read whose call
    /// at `position` shows another base than the reference with the
    /// reference base as a second peak above its threshold (SRS-CALL-020;
    /// every read of the strand since rules version 33, ADR-0090): a check,
    /// one per position and read.
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
    /// Runs of a profile length tract whose dominant length the reads
    /// determine, or that DNA does not report: a length heteroplasmy noted,
    /// context and never a check (SRS-API-020).
    TractLength {
        tract: String,
        /// The window a database search disregards length variants in
        /// (EMPOP).
        ignored_in_search: [usize; 2],
        runs: Vec<TractRunResult>,
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
    /// A stretch at a reportable region's end that the call does not
    /// analyse (SRS-CONS-012), with its reason (rules version 36): a check,
    /// so that an allele the report omits there is never passed silently.
    NotAnalysed {
        reason: &'static str,
    },
}

impl LabelResult {
    /// The label's kind.
    pub(crate) const fn kind(&self) -> LabelKind {
        match self {
            Self::StrandsDisagree { .. } => LabelKind::StrandsDisagree,
            Self::ReadDisagrees { .. } => LabelKind::ReadDisagrees,
            Self::SameStrandDissent { .. } => LabelKind::SameStrandDissent,
            Self::PossibleHeteroplasmy { .. } => LabelKind::PossibleHeteroplasmy,
            Self::ReferenceSecond { .. } => LabelKind::ReferenceSecond,
            Self::EstimatedLength { .. } => LabelKind::EstimatedLength,
            Self::LengthUnconfirmed { .. } => LabelKind::LengthUnconfirmed,
            Self::Undecided => LabelKind::Undecided,
            Self::NotCalled => LabelKind::NotCalled,
            Self::LongRun { .. } => LabelKind::LongRun,
            Self::EdgeDeletion { .. } => LabelKind::EdgeDeletion,
            Self::RunDeletion { .. } => LabelKind::RunDeletion,
            Self::RunIndel { .. } => LabelKind::RunIndel,
            Self::OneStrandInsertion { .. } => LabelKind::OneStrandInsertion,
            Self::BesideNotAnalysed => LabelKind::BesideNotAnalysed,
            Self::NotAnalysed { .. } => LabelKind::NotAnalysed,
            Self::OneRead { .. } => LabelKind::OneRead,
            Self::OneStrand { .. } => LabelKind::OneStrand,
            Self::NoisySignal { .. } => LabelKind::NoisySignal,
            Self::NotCovered => LabelKind::NotCovered,
            Self::PhaseShift { .. } => LabelKind::PhaseShift,
            Self::DephasedEnd { .. } => LabelKind::DephasedEnd,
            Self::ReadEdge { .. } => LabelKind::ReadEdge,
            Self::PoorSignal { .. } => LabelKind::PoorSignal,
            Self::LengthNotRecorded { .. } => LabelKind::LengthNotRecorded,
            Self::TractLength { .. } => LabelKind::TractLength,
        }
    }

    /// The label's rank: its kind's, and [`CONTEXT_RANK`] for an `explained`
    /// possible heteroplasmy.
    pub(crate) const fn rank(&self) -> u8 {
        match self {
            Self::PossibleHeteroplasmy {
                explained: true, ..
            } => CONTEXT_RANK,
            _ => self.kind().rank(),
        }
    }

    /// Whether the label leaves its item unassured: its kind's, and never
    /// for an `explained` possible heteroplasmy, whose second peaks are
    /// context.
    pub(crate) const fn unassured(&self) -> bool {
        match self {
            Self::PossibleHeteroplasmy {
                explained: true, ..
            } => false,
            _ => self.kind().check(),
        }
    }
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
/// whose length rests on an anchored end or on one strand (ADR-0090).
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
/// after it, 1 alone when the anchor is lost), its base, DNA's length (none
/// when DNA does not report it), how the reads' decodes agree, and each
/// read's length.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub(crate) struct LengthRunResult {
    pub(crate) run: usize,
    pub(crate) base: char,
    pub(crate) called: Option<usize>,
    pub(crate) agreement: &'static str,
    /// The lab's convention reading of a near-equal pair, a suggestion only:
    /// never DNA's call.
    pub(crate) convention: Option<LengthConventionResult>,
    /// The reads' floor reading of a run DNA reports no length for: the
    /// least length any read gives as its top or next. A suggestion only.
    pub(crate) floor: Option<usize>,
    pub(crate) reads: Vec<LengthReadResult>,
}

/// One run of a length tract, noted: DNA's reported length, how the run
/// stands, the length a reviewer would record, and the reads' evidence.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub(crate) struct TractRunResult {
    pub(crate) run: usize,
    pub(crate) base: char,
    pub(crate) called: Option<usize>,
    /// `confirmed`, `dominant`, `convention`, or `unconfirmed`.
    pub(crate) status: &'static str,
    /// The length a reviewer would record: DNA's, the convention's, or the
    /// floor reading for an unconfirmed run; `None` where the profile
    /// records no length or no read reaches the run.
    pub(crate) reading: Option<usize>,
    /// How the reads' lengths agree, as `length_unconfirmed` gives it;
    /// `None` when no read reaches the run.
    pub(crate) agreement: Option<&'static str>,
    pub(crate) strands: StrandsResult,
    /// The reads that show the run in phase.
    pub(crate) in_phase: Vec<String>,
    pub(crate) floor: Option<usize>,
    pub(crate) convention: Option<LengthConventionResult>,
    pub(crate) reads: Vec<TractReadResult>,
}

/// Each strand's one clear length of a run, when it gives one.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize)]
pub(crate) struct StrandsResult {
    pub(crate) forward: Option<usize>,
    pub(crate) reverse: Option<usize>,
}

/// One read's length of a noted run with its decode's components.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub(crate) struct TractReadResult {
    pub(crate) read: String,
    pub(crate) strand: Orientation,
    pub(crate) length: usize,
    pub(crate) source: &'static str,
    /// The decode's first-run and total length distributions; `None` for
    /// a read in phase.
    pub(crate) first: Option<MarginalResult>,
    pub(crate) total: Option<MarginalResult>,
}

/// A decode's length distribution: signal shares, never heteroplasmy
/// levels, and the dominant's margin over the next in standard errors.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub(crate) struct MarginalResult {
    pub(crate) lengths: Vec<usize>,
    pub(crate) shares: Vec<f64>,
    pub(crate) z: Option<f64>,
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

/// Items carrying each label, by kind, serialized in [`LabelKind::ALL`]
/// order.
#[derive(Debug, Default, PartialEq, Eq)]
pub(crate) struct LabelCounts([usize; LabelKind::ALL.len()]);

impl LabelCounts {
    /// Counts one more item carrying a label of `kind`.
    pub(crate) const fn add(&mut self, kind: LabelKind) {
        self.0[kind as usize] += 1;
    }

    /// Items carrying a label of `kind`.
    #[cfg(test)]
    pub(crate) const fn get(&self, kind: LabelKind) -> usize {
        self.0[kind as usize]
    }
}

impl Serialize for LabelCounts {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        use serde::ser::SerializeMap;
        let mut map = serializer.serialize_map(Some(LabelKind::ALL.len()))?;
        for kind in LabelKind::ALL {
            map.serialize_entry(kind.name(), &self.0[kind as usize])?;
        }
        map.end()
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;
    use std::path::Path;

    use super::*;

    /// The label table agrees with the published schema: the summary lists
    /// every kind in `ALL` order with its rank and check, and the labels an
    /// item can carry are exactly those kinds.
    // SRS-API-014
    #[test]
    fn the_label_table_is_the_schemas() -> Result<(), Box<dyn std::error::Error>> {
        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("docs/reference/schemas/review-v1.schema.json");
        let schema: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(path)?)?;
        let summary = &schema["$defs"]["summary"]["properties"]["labels"];
        let required = summary["required"]
            .as_array()
            .ok_or("summary labels must list their names")?;
        assert_eq!(
            required
                .iter()
                .filter_map(serde_json::Value::as_str)
                .collect::<Vec<_>>(),
            LabelKind::ALL.map(LabelKind::name)
        );
        for (index, kind) in LabelKind::ALL.into_iter().enumerate() {
            assert_eq!(kind as usize, index);
            let entry = &summary["properties"][kind.name()];
            assert_eq!(entry["x-rank"], u64::from(kind.rank()), "{}", kind.name());
            assert_eq!(entry["x-check"], kind.check(), "{}", kind.name());
        }
        let labels = schema["$defs"]["label"]["oneOf"]
            .as_array()
            .ok_or("labels must be a oneOf")?
            .iter()
            .filter_map(|label| label["properties"]["label"]["const"].as_str())
            .collect::<BTreeSet<_>>();
        assert_eq!(
            labels,
            LabelKind::ALL.map(LabelKind::name).into_iter().collect()
        );
        Ok(())
    }
}
