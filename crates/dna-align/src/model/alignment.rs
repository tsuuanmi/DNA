//! Pairwise alignment records with explicit strand and coordinates.

use std::collections::BTreeSet;

use serde::{Deserialize, Serialize};

use dna_kernel::error::{AlignmentError, Result};
use dna_kernel::model::nucleotide::is_canonical;

/// Query orientation relative to the supplied reference strand.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Orientation {
    /// Native retained query matches the reference strand.
    Forward,
    /// Reverse-complemented retained query matches the reference strand.
    Reverse,
}

impl Orientation {
    /// Projects one trace-strand base onto the reference strand: a reverse
    /// read's canonical base is complemented; any other symbol is kept.
    #[must_use]
    pub const fn reference_base(self, base: char) -> char {
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
    /// Query bases from masked calls: a call whose mask anchors the alignment
    /// aligns with its call and profile, every other masked call as
    /// unresolved.
    pub masked_query_bases: usize,
}

impl AlignmentMetrics {
    /// Callable columns whose bases agree.
    #[must_use]
    pub const fn exact_matches(&self) -> usize {
        self.exact_matches
    }

    /// Callable columns whose bases differ.
    #[must_use]
    pub const fn mismatches(&self) -> usize {
        self.mismatches
    }
}

/// One column of an alignment: a read base, or `-` where the read has none,
/// against a reference base, or `-` where the reference has none.
#[derive(Debug, Clone)]
pub struct AlignmentColumn {
    pub(crate) query_base: char,
    pub(crate) reference_base: char,
    pub(crate) original_call_index_0based: Option<usize>,
    pub(crate) reference_index_0based: Option<usize>,
}

impl AlignmentColumn {
    /// A column pairing `query_base`, the read's call at
    /// `original_call_index_0based` in the read's own order, with
    /// `reference_base` at `reference_index_0based`; a gap has no index on
    /// its side.
    #[must_use]
    pub const fn new(
        query_base: char,
        reference_base: char,
        original_call_index_0based: Option<usize>,
        reference_index_0based: Option<usize>,
    ) -> Self {
        Self {
            query_base,
            reference_base,
            original_call_index_0based,
            reference_index_0based,
        }
    }

    /// The read's base on the reference strand, or `-` in a deletion.
    #[must_use]
    pub const fn query_base(&self) -> char {
        self.query_base
    }

    /// The reference base, or `-` in an insertion.
    #[must_use]
    pub const fn reference_base(&self) -> char {
        self.reference_base
    }

    /// The read's call, in the read's own order, or `None` in a deletion.
    #[must_use]
    pub const fn original_call_index_0based(&self) -> Option<usize> {
        self.original_call_index_0based
    }

    /// The reference position (0-based), or `None` in an insertion.
    #[must_use]
    pub const fn reference_index_0based(&self) -> Option<usize> {
        self.reference_index_0based
    }
}

/// Selected alignment and both orientation summaries.
#[derive(Debug, Clone)]
///
/// Only [`align_best`](crate::align_best), [`Alignment::placed`], and the
/// placed-read contract's reader ([`placed::read`](crate::placed::read), of a
/// document [`align_best`](crate::align_best)'s caller wrote) build one,
/// and only its own methods narrow its callable segments, so its invariants
/// hold wherever it is read.
pub struct Alignment {
    pub(crate) orientation: Orientation,
    pub(crate) score: i64,
    pub(crate) reference_segments: Vec<ReferenceSegment>,
    pub(crate) callable_segments: Vec<ReferenceSegment>,
    pub(crate) wraps_origin: bool,
    pub(crate) metrics: AlignmentMetrics,
    pub(crate) columns: Vec<AlignmentColumn>,
}

impl Alignment {
    /// An alignment its source places rather than one [`align_best`]
    /// searches: a forward read whose columns are given, on a reference it
    /// does not wrap. Its score is zero, and its metrics follow from the
    /// columns: every non-canonical read base is unresolved and none is
    /// masked.
    ///
    /// [`align_best`]: crate::align_best
    ///
    /// # Errors
    ///
    /// Returns `AlignmentError::Inconsistent` unless the read's calls and the
    /// reference positions both increase along the columns, every reference
    /// position lies in `reference_segments`, every callable segment lies in
    /// one of them, and the segments are ordered and disjoint.
    pub fn placed(
        reference_segments: Vec<ReferenceSegment>,
        callable_segments: Vec<ReferenceSegment>,
        columns: Vec<AlignmentColumn>,
    ) -> Result<Self> {
        let ordered = |segments: &[ReferenceSegment]| {
            segments
                .iter()
                .all(|segment| segment.start_0based < segment.end_0based_exclusive)
                && segments
                    .windows(2)
                    .all(|pair| pair[0].end_0based_exclusive <= pair[1].start_0based)
        };
        if !ordered(&reference_segments) || !ordered(&callable_segments) {
            return Err(AlignmentError::Inconsistent(
                "placed alignment segments are not ordered and disjoint",
            )
            .into());
        }
        let within = |start: usize, end: usize| {
            reference_segments
                .iter()
                .any(|segment| segment.start_0based <= start && end <= segment.end_0based_exclusive)
        };
        if !callable_segments
            .iter()
            .all(|segment| within(segment.start_0based, segment.end_0based_exclusive))
        {
            return Err(AlignmentError::Inconsistent(
                "placed alignment has a callable segment outside its reference segments",
            )
            .into());
        }
        if !increasing(
            columns
                .iter()
                .filter_map(|column| column.original_call_index_0based),
        ) || !increasing(
            columns
                .iter()
                .filter_map(|column| column.reference_index_0based),
        ) {
            return Err(AlignmentError::Inconsistent(
                "placed alignment calls or reference positions do not increase",
            )
            .into());
        }
        if !columns
            .iter()
            .filter_map(|column| column.reference_index_0based)
            .all(|position| within(position, position + 1))
        {
            return Err(AlignmentError::Inconsistent(
                "placed alignment reads a position outside its reference segments",
            )
            .into());
        }
        let metrics = placed_metrics(&columns);
        Ok(Self {
            orientation: Orientation::Forward,
            score: 0,
            reference_segments,
            callable_segments,
            wraps_origin: false,
            metrics,
            columns,
        })
    }

    /// Read orientation selected against the reference.
    #[must_use]
    pub const fn orientation(&self) -> Orientation {
        self.orientation
    }

    /// Reference segments covered by the alignment, in alignment order.
    #[must_use]
    pub fn reference_segments(&self) -> &[ReferenceSegment] {
        &self.reference_segments
    }

    /// Reference segments observed by unmasked calls, and by deletions between
    /// them, in alignment order; a subset of the reference segments.
    #[must_use]
    pub fn callable_segments(&self) -> &[ReferenceSegment] {
        &self.callable_segments
    }

    /// Whether the alignment crosses the origin of a circular reference.
    #[must_use]
    pub const fn wraps_origin(&self) -> bool {
        self.wraps_origin
    }

    /// Summary metrics of the alignment.
    #[must_use]
    pub const fn metrics(&self) -> &AlignmentMetrics {
        &self.metrics
    }

    /// The alignment's columns, in reference order.
    #[must_use]
    pub fn columns(&self) -> &[AlignmentColumn] {
        &self.columns
    }

    /// The same alignment with the reference positions `excluded` (0-based)
    /// left out of its callable segments, each segment split around them.
    #[must_use]
    pub fn without_callable(mut self, excluded: &BTreeSet<usize>) -> Self {
        self.callable_segments = self
            .callable_segments
            .iter()
            .flat_map(|segment| {
                let mut pieces = Vec::new();
                let mut start = segment.start_0based;
                for position in excluded
                    .range(segment.start_0based..segment.end_0based_exclusive)
                    .copied()
                {
                    if position > start {
                        pieces.push(ReferenceSegment {
                            start_0based: start,
                            end_0based_exclusive: position,
                        });
                    }
                    start = position + 1;
                }
                if segment.end_0based_exclusive > start {
                    pieces.push(ReferenceSegment {
                        start_0based: start,
                        end_0based_exclusive: segment.end_0based_exclusive,
                    });
                }
                pieces
            })
            .collect();
        self
    }

    /// The same alignment with no callable segment: the placement itself is
    /// in doubt, so nothing it reads is callable.
    #[must_use]
    pub fn uncallable(mut self) -> Self {
        self.callable_segments.clear();
        self
    }

    /// The alignment's profile score; zero for a placed alignment.
    #[must_use]
    pub const fn score(&self) -> i64 {
        self.score
    }
}

#[cfg(any(test, feature = "test-support"))]
impl Alignment {
    /// An alignment of any shape, unchecked, so tests can reach the guards of
    /// the stages that read alignments. Its score is zero and its metrics
    /// follow from the columns as a placed alignment's do.
    #[must_use]
    pub fn fixture(
        orientation: Orientation,
        reference_segments: Vec<ReferenceSegment>,
        callable_segments: Vec<ReferenceSegment>,
        columns: Vec<AlignmentColumn>,
    ) -> Self {
        Self {
            orientation,
            score: 0,
            reference_segments,
            callable_segments,
            wraps_origin: false,
            metrics: placed_metrics(&columns),
            columns,
        }
    }
}

/// Whether every index exceeds the one before it.
fn increasing(indices: impl Iterator<Item = usize>) -> bool {
    let mut previous = None;
    for index in indices {
        if previous.is_some_and(|before| before >= index) {
            return false;
        }
        previous = Some(index);
    }
    true
}

/// The metrics of a placed alignment's columns.
fn placed_metrics(columns: &[AlignmentColumn]) -> AlignmentMetrics {
    let mut metrics = AlignmentMetrics {
        exact_matches: 0,
        mismatches: 0,
        gap_opens: 0,
        callable_columns: 0,
        callable_identity: 0.0,
        unresolved_query_bases: 0,
        masked_query_bases: 0,
    };
    let mut previous = None;
    for column in columns {
        let gap = match (column.query_base, column.reference_base) {
            ('-', _) => Some('D'),
            (_, '-') => Some('I'),
            _ => None,
        };
        if gap.is_some() && gap != previous {
            metrics.gap_opens += 1;
        }
        previous = gap;
        if column.query_base != '-' && !is_canonical(column.query_base) {
            metrics.unresolved_query_bases += 1;
        }
        if is_canonical(column.query_base) && is_canonical(column.reference_base) {
            metrics.callable_columns += 1;
            if column.query_base == column.reference_base {
                metrics.exact_matches += 1;
            } else {
                metrics.mismatches += 1;
            }
        }
    }
    if metrics.callable_columns > 0 {
        #[allow(clippy::cast_precision_loss)]
        let identity = metrics.exact_matches as f64 / metrics.callable_columns as f64;
        metrics.callable_identity = identity;
    }
    metrics
}

#[cfg(test)]
mod tests {
    use super::*;

    fn segment(start_0based: usize, end_0based_exclusive: usize) -> ReferenceSegment {
        ReferenceSegment {
            start_0based,
            end_0based_exclusive,
        }
    }

    fn placed() -> Result<Alignment> {
        let columns = (10..20)
            .enumerate()
            .map(|(call, position)| AlignmentColumn::new('A', 'A', Some(call), Some(position)))
            .collect();
        Alignment::placed(
            vec![segment(10, 20)],
            vec![segment(10, 15), segment(16, 20)],
            columns,
        )
    }

    #[test]
    fn cutting_positions_splits_the_callable_segments_only() -> Result<()> {
        let cut = placed()?.without_callable(&BTreeSet::from([10, 12, 15, 19, 30]));
        assert_eq!(
            cut.callable_segments(),
            [segment(11, 12), segment(13, 15), segment(16, 19)]
        );
        assert_eq!(cut.reference_segments(), [segment(10, 20)]);
        assert!(placed()?.uncallable().callable_segments().is_empty());
        Ok(())
    }

    #[test]
    fn a_placed_alignment_rejects_a_callable_segment_outside_its_reach() {
        let columns = vec![AlignmentColumn::new('A', 'A', Some(0), Some(10))];
        assert!(Alignment::placed(vec![segment(10, 11)], vec![segment(10, 12)], columns).is_err());
    }
}
