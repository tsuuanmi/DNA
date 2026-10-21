//! An adjudicated sample consensus: sequence segments and the evidence behind
//! each decision (PROP-0003, ADR-0073).

use std::collections::{BTreeMap, BTreeSet};

use dna_kernel::read_evidence::EvidenceReason;

/// How the rule decided a site or stretch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiteState {
    /// Every clean observation agrees, or a strict majority decided.
    Called,
    /// Clean observations tie: a one-against-one tie takes the reference
    /// sequence; any other tie leaves none.
    Contested,
    /// No read observes the site cleanly.
    Unresolved,
}

impl SiteState {
    /// Stable published label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::Called => "called",
            Self::Contested => "contested",
            Self::Unresolved => "unresolved",
        }
    }
}

/// How a decided run length is known (ADR-0074).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RunEvidence {
    /// A read shows both ends of the run with clean calls.
    InPhase,
    /// A read shows both ends, one of them only through a masked call that
    /// still anchors the alignment (`Placement::Anchoring`) and reads
    /// the next base. The length is at most the dominant one of a mixture.
    AnchoredEnd,
    /// No read shows both ends, no reads locate both, and a read that shows
    /// one end loses phase at the other: the length is the longest run read
    /// in phase before losing phase. Within a length mixture this lies near
    /// the dominant length; it is an estimate.
    PhaseLoss,
    /// No read shows both ends, but a read enters the run in phase from each
    /// end, every call of the run unmasked, and both lose phase at the other
    /// end, on a call that does not read the run's base, after the same
    /// number of calls: each read's count stands for the end the other read
    /// cannot show.
    ConfirmedEnd,
    /// No read shows both ends, but reads locate each end on the reference
    /// and read every base between them cleanly: the length is the number of
    /// reference positions between those ends. It assumes no length change
    /// and is not observed.
    ReferenceFrame,
    /// No read shows both ends, but reads of both strands decode the run's
    /// length from their own signal, each beyond the strands' spread and
    /// its strand's replicate noise (SRS-CONS-013): the length
    /// replaces an estimate, never a length a read shows.
    TwoStrand,
}

impl RunEvidence {
    /// Stable published label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::InPhase => "in_phase",
            Self::AnchoredEnd => "anchored_end",
            Self::PhaseLoss => "phase_loss",
            Self::ConfirmedEnd => "confirmed_end",
            Self::ReferenceFrame => "reference_frame",
            Self::TwoStrand => "two_strand",
        }
    }

    /// The length is an estimate, not observed: `phase_loss` or
    /// `reference_frame`. It keeps its label even where it equals the
    /// reference.
    #[must_use]
    pub const fn is_estimate(self) -> bool {
        matches!(self, Self::PhaseLoss | Self::ReferenceFrame)
    }

    /// A read shows both of the run's ends, so its bounds place the run:
    /// `in_phase` or `anchored_end`.
    #[must_use]
    pub const fn shown_by_read(self) -> bool {
        matches!(self, Self::InPhase | Self::AnchoredEnd)
    }

    /// The two-strand rule may decide the length instead (SRS-CONS-013): an
    /// `anchored_end` length, at most the dominant one of a mixture, or an
    /// estimate.
    #[must_use]
    pub const fn yields_to_two_strand(self) -> bool {
        matches!(self, Self::AnchoredEnd) || self.is_estimate()
    }
}

/// One run of a stretch decided run by run: its base, its length, and how
/// the length is known.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ConsensusRun {
    /// The run's base.
    pub base: char,
    /// The run's length.
    pub length: usize,
    /// How the length is known.
    pub evidence: RunEvidence,
}

/// One decided reference interval: a single position, or a stretch of whole
/// reference runs around reads' differences, insertions at its edges
/// included.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConsensusSite {
    /// First reference position (0-based, inclusive).
    pub first_0based: usize,
    /// Last reference position (0-based, inclusive).
    pub last_0based: usize,
    /// Reference sequence over the interval.
    pub reference: String,
    /// Decided sequence, or `None` when undecided.
    pub call: Option<String>,
    /// How the sequence was decided.
    pub state: SiteState,
    /// Reads whose clean sequence matches the call. Every read list names a
    /// read by its identity, its input's content SHA-256, never by its file
    /// name, which two inputs can share.
    pub supporting: Vec<String>,
    /// Reads whose clean sequence differs from the call.
    pub opposing: Vec<String>,
    /// Reads that cover the interval without a clean sequence.
    pub uninformative: Vec<String>,
    /// The call's runs when the stretch was decided run by run; empty for a
    /// single position, an undecided stretch, or one decided position by
    /// position.
    pub runs: Vec<ConsensusRun>,
}

/// A run of covered reference positions and its consensus sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConsensusSegment {
    /// First reference position (0-based, inclusive).
    pub start_0based: usize,
    /// Last reference position (0-based, inclusive).
    pub end_0based: usize,
    /// Consensus bases: `N` where undecided, insertions spliced in,
    /// deletions left out.
    pub sequence: String,
}

/// Counts of decisions over the consensus segments.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ConsensusSummary {
    /// Reference positions decided without contest.
    pub called_positions: usize,
    /// Positions and stretches decided by a tie.
    pub contested_sites: usize,
    /// Reference positions written as `N`.
    pub unresolved_positions: usize,
    /// Runs whose length is the reference frame's, not observed.
    pub reference_frame_runs: usize,
    /// Runs whose length is an in-phase estimate where reads lose phase.
    pub phase_loss_runs: usize,
}

/// A sample consensus over the admitted reads.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Consensus {
    /// Consensus segments in reference order.
    pub segments: Vec<ConsensusSegment>,
    /// Sites whose call differs from the reference, that were not called
    /// outright, or that had disagreeing clean observations.
    pub sites: Vec<ConsensusSite>,
    /// Decision counts.
    pub summary: ConsensusSummary,
    /// Reference positions (0-based) that only reference confirmations
    /// observe cleanly (ADR-0080): decided as the reference, they bound no
    /// neighbour's run when the decisions are called.
    pub reference_only: BTreeSet<usize>,
    /// Reference positions (0-based) that reads observe only inside their
    /// unreliable spans (SRS-CONS-012): some read's call there lies in such
    /// a span, and no read shows the position cleanly outside its own. Each
    /// with the reason of the first such span.
    pub unreliable_only: BTreeMap<usize, EvidenceReason>,
}
