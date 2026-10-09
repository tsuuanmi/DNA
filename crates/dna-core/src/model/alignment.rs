//! Pairwise alignment records with explicit strand and coordinates.

use serde::Serialize;

/// Query orientation relative to the supplied reference strand.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Orientation {
    /// Native retained query matches the reference strand.
    Forward,
    /// Reverse-complemented retained query matches the reference strand.
    Reverse,
}

impl Orientation {
    /// Projects one trace-strand canonical base onto the reference strand.
    pub(crate) const fn reference_base(self, base: char) -> char {
        match self {
            Self::Forward => base,
            Self::Reverse => match base {
                'A' => 'T',
                'C' => 'G',
                'G' => 'C',
                'T' => 'A',
                other => other,
            },
        }
    }
}

/// One half-open segment on the original reference.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReferenceSegment {
    /// First covered reference position (0-based, inclusive).
    pub start_0based: usize,
    /// Position after the last covered reference position (0-based, exclusive).
    pub end_0based_exclusive: usize,
}

/// Alignment quality metrics.
#[derive(Debug, Clone, Serialize)]
pub struct AlignmentMetrics {
    pub(crate) exact_matches: usize,
    pub(crate) mismatches: usize,
    /// Number of gap openings in the selected traceback.
    pub gap_opens: usize,
    /// Columns pairing a canonical read base with a canonical reference base.
    pub callable_columns: usize,
    /// Fraction of callable columns whose bases agree.
    pub callable_identity: f64,
    /// Unresolved query bases, excluding masked calls.
    pub unresolved_query_bases: usize,
    /// Query bases from masked calls: dephased calls align with their call
    /// and profile, every other masked call as unresolved.
    pub masked_query_bases: usize,
}

/// One column of the selected alignment.
#[derive(Debug, Clone)]
pub(crate) struct AlignmentColumn {
    pub(crate) query_base: char,
    pub(crate) reference_base: char,
    pub(crate) original_call_index_0based: Option<usize>,
    pub(crate) reference_index_0based: Option<usize>,
}

/// Selected alignment and both orientation summaries.
#[derive(Debug, Clone)]
pub struct Alignment {
    /// Read orientation selected against the reference.
    pub orientation: Orientation,
    pub(crate) score: i64,
    /// Reference segments covered by the alignment, in alignment order.
    pub reference_segments: Vec<ReferenceSegment>,
    /// Reference segments observed by unmasked calls, and by deletions between
    /// them, in alignment order; a subset of `reference_segments`.
    pub callable_segments: Vec<ReferenceSegment>,
    /// Whether the alignment crosses the origin of a circular reference.
    pub wraps_origin: bool,
    /// Summary metrics of the selected alignment.
    pub metrics: AlignmentMetrics,
    pub(crate) columns: Vec<AlignmentColumn>,
}
