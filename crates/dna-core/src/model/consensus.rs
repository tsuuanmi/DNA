//! An adjudicated sample consensus: sequence segments and the evidence behind
//! each decision (PROP-0003, ADR-0073).

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
    /// Reads whose clean sequence matches the call.
    pub supporting: Vec<String>,
    /// Reads whose clean sequence differs from the call.
    pub opposing: Vec<String>,
    /// Reads that cover the interval without a clean sequence.
    pub uninformative: Vec<String>,
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
}
