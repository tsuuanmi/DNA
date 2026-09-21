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

    /// Projects A/C/G/T channel heights from trace strand to reference strand.
    pub(crate) const fn reference_peak_heights(self, peaks: [i32; 4]) -> [i32; 4] {
        match self {
            Self::Forward => peaks,
            Self::Reverse => [peaks[3], peaks[2], peaks[1], peaks[0]],
        }
    }

    /// Projects A/C/G/T signal metrics from trace strand to reference strand.
    pub(crate) const fn reference_signal_values(self, values: [f64; 4]) -> [f64; 4] {
        match self {
            Self::Forward => values,
            Self::Reverse => [values[3], values[2], values[1], values[0]],
        }
    }
}

/// One half-open segment on the original reference.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReferenceSegment {
    pub(crate) start_0based: usize,
    pub(crate) end_0based_exclusive: usize,
}

/// Alignment quality metrics.
#[derive(Debug, Clone, Serialize)]
pub struct AlignmentMetrics {
    pub(crate) exact_matches: usize,
    pub(crate) mismatches: usize,
    pub(crate) gap_opens: usize,
    pub(crate) callable_columns: usize,
    pub(crate) callable_identity: f64,
    pub(crate) unresolved_query_bases: usize,
}

/// One column of the selected alignment.
#[derive(Debug, Clone)]
pub struct AlignmentColumn {
    pub(crate) query_base: char,
    pub(crate) reference_base: char,
    pub(crate) original_call_index_0based: Option<usize>,
    pub(crate) reference_index_0based: Option<usize>,
}

/// Selected alignment and both orientation summaries.
#[derive(Debug, Clone)]
pub struct Alignment {
    pub(crate) orientation: Orientation,
    pub(crate) score: i64,
    pub(crate) reference_segments: Vec<ReferenceSegment>,
    pub(crate) wraps_origin: bool,
    pub(crate) metrics: AlignmentMetrics,
    pub(crate) columns: Vec<AlignmentColumn>,
}
