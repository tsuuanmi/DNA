//! Run-structure re-expression of repeat stretches (ADR-0072).
//!
//! Masked calls keep anchoring the optimal alignment, but their mixed signal
//! can make it split one run-length change into a substitution next to an
//! insertion, or merge two of them into one gap. A clean read is merged the
//! same way whenever a gap opening costs more than a mismatch and a lost
//! match: two run-length changes become one gap beside a substitution. Where
//! the read's calls spell the reference's runs in the same order, every run
//! only longer or every run only shorter, and a masked call or a substitution
//! lies among them, the selected placement is re-expressed as one length edit
//! at the 3' end of each changed run, so the read shows every reference base
//! it carries. The placement, its score, and every column outside such a
//! stretch are kept.

use crate::traceback::{self, RawAlignment, RawColumn};
use dna_kernel::model::reference::Reference;

/// Columns between two anchors: `left` and `right` are the anchor columns and
/// `first..=last` the working-reference positions strictly between them.
#[derive(Debug, Clone, Copy)]
struct Stretch {
    left: usize,
    right: usize,
    first: usize,
    last: usize,
}

/// Re-expresses every qualifying stretch of `alignment` in place; `masked`
/// tells whether the call at an oriented query index is masked. A stretch
/// spans at most `max_stretch` reference bases, which bounds the work per
/// disturbed block and keeps a stretch local to one repeat.
pub(super) fn reexpress(
    alignment: &mut RawAlignment,
    reference: &Reference,
    masked: &dyn Fn(usize) -> bool,
    max_stretch: usize,
) {
    let stretches = stretches(&alignment.columns, reference, max_stretch);
    let mut columns = Vec::with_capacity(alignment.columns.len());
    let mut cursor = 0;
    let mut changed = false;
    for stretch in stretches {
        columns.extend_from_slice(&alignment.columns[cursor..=stretch.left]);
        let inside = &alignment.columns[stretch.left + 1..stretch.right];
        if let Some(rebuilt) = rebuild(inside, stretch, reference, masked) {
            columns.extend(rebuilt);
            changed = true;
        } else {
            columns.extend_from_slice(inside);
        }
        cursor = stretch.right;
    }
    if changed {
        columns.extend_from_slice(&alignment.columns[cursor..]);
        alignment.columns = columns;
        alignment.metrics = traceback::metrics(&alignment.columns);
    }
}

/// The reference base at a working-reference position.
fn base(reference: &Reference, index: usize) -> u8 {
    let bytes = reference.sequence.as_bytes();
    bytes[index % bytes.len()]
}

/// Whether two working-reference positions lie in the same copy of a circular
/// reference, so a stretch never crosses the origin seam.
fn same_copy(reference: &Reference, left: usize, right: usize) -> bool {
    left / reference.len() == right / reference.len()
}

fn run_start(reference: &Reference, mut index: usize) -> usize {
    while index > 0
        && same_copy(reference, index - 1, index)
        && base(reference, index - 1) == base(reference, index)
    {
        index -= 1;
    }
    index
}

fn run_end(reference: &Reference, mut index: usize) -> usize {
    while same_copy(reference, index, index + 1)
        && base(reference, index + 1) == base(reference, index)
    {
        index += 1;
    }
    index
}

const fn exact(column: &RawColumn) -> bool {
    column.query_base == column.reference_base && matches!(column.query_base, 'A' | 'C' | 'G' | 'T')
}

/// Disturbed blocks grown to whole reference runs and bounded by anchors:
/// exact matches whose base the read's adjacent call inside does not repeat.
fn stretches(columns: &[RawColumn], reference: &Reference, max_stretch: usize) -> Vec<Stretch> {
    let Some(origin) = columns.iter().find_map(|column| column.reference_index) else {
        return Vec::new();
    };
    let mut at = Vec::new();
    for (index, column) in columns.iter().enumerate() {
        if let Some(position) = column.reference_index {
            at.resize(position - origin + 1, None);
            at[position - origin] = Some(index);
        }
    }
    let column_of = |position: usize| {
        position
            .checked_sub(origin)
            .and_then(|offset| at.get(offset).copied().flatten())
    };
    let mut found: Vec<Stretch> = Vec::new();
    let mut index = 0;
    while index < columns.len() {
        if exact(&columns[index]) {
            index += 1;
            continue;
        }
        let start = index;
        while index < columns.len() && !exact(&columns[index]) {
            index += 1;
        }
        // An insertion also touches the reference positions around it.
        let mut footprint = Vec::new();
        for offset in start..index {
            if let Some(position) = columns[offset].reference_index {
                footprint.push(position);
            } else {
                let mut previous = columns[..offset].iter().rev();
                footprint.extend(previous.find_map(|column| column.reference_index));
                let mut following = columns[offset + 1..].iter();
                footprint.extend(following.find_map(|column| column.reference_index));
            }
        }
        let (Some(&low), Some(&high)) = (footprint.iter().min(), footprint.iter().max()) else {
            continue;
        };
        let Some(stretch) = anchored(columns, reference, &column_of, low, high, max_stretch) else {
            continue;
        };
        // Stretches whose anchors are each other's first and last columns
        // share no column between their anchors: each stays its own.
        match found.last_mut() {
            Some(previous) if stretch.left + 1 < previous.right => {
                previous.right = previous.right.max(stretch.right);
                previous.last = previous.last.max(stretch.last);
            }
            _ => found.push(stretch),
        }
        index = index.max(stretch.right);
    }
    found
}

/// Grows `low..=high` run by run until anchors bound it on both sides.
fn anchored(
    columns: &[RawColumn],
    reference: &Reference,
    column_of: &dyn Fn(usize) -> Option<usize>,
    low: usize,
    high: usize,
    max_stretch: usize,
) -> Option<Stretch> {
    let (mut first, mut last) = (run_start(reference, low), run_end(reference, high));
    while last - first < max_stretch && same_copy(reference, first, last) {
        let left = column_of(first.checked_sub(1)?)?;
        let right = column_of(last + 1)?;
        let mut inside = columns[left + 1..right]
            .iter()
            .filter(|column| column.query_index.is_some())
            .map(|column| column.query_base);
        let (head, tail) = (inside.next(), inside.next_back());
        let tail = tail.or(head);
        let left_holds =
            exact(&columns[left]) && head.is_some_and(|call| call != columns[left].query_base);
        let right_holds =
            exact(&columns[right]) && tail.is_some_and(|call| call != columns[right].query_base);
        if left_holds && right_holds {
            return Some(Stretch {
                left,
                right,
                first,
                last,
            });
        }
        if !left_holds {
            first = run_start(reference, first.checked_sub(1)?);
        }
        if !right_holds {
            last = run_end(reference, last + 1);
        }
    }
    None
}

/// Runs of equal symbols as `(symbol, length)`.
fn runs(symbols: impl Iterator<Item = char>) -> Vec<(char, usize)> {
    let mut runs: Vec<(char, usize)> = Vec::new();
    for symbol in symbols {
        match runs.last_mut() {
            Some((previous, length)) if *previous == symbol => *length += 1,
            _ => runs.push((symbol, 1)),
        }
    }
    runs
}

/// A stretch's run-length changes, read run against reference run. Where one
/// run is shorter and its neighbour longer, the boundary between them moves:
/// the shorter run's positions beside it are substitutions by the longer
/// run's base (`16182C 16183C`). What remains of each run's change is one
/// length edit at its 3' end.
struct Edits {
    /// Per boundary between runs `i` and `i + 1`: how far read run `i`
    /// reaches into reference run `i + 1` (positive), or read run `i + 1`
    /// back into reference run `i` (negative).
    shifts: Vec<isize>,
    /// Per run: the inserted (positive) or deleted (negative) bases left
    /// after the boundary moves.
    residual: Vec<isize>,
}

impl Edits {
    /// Moves boundaries from 5' to 3', each as far as both of its runs'
    /// changes allow. `None` when the remaining changes still both lengthen
    /// and shorten runs: substitutions elsewhere explain such a stretch.
    fn of(read: &[usize], expected: &[usize]) -> Option<Self> {
        let mut residual = read
            .iter()
            .zip(expected)
            .map(|(&observed, &wanted)| {
                Some(isize::try_from(observed).ok()? - isize::try_from(wanted).ok()?)
            })
            .collect::<Option<Vec<_>>>()?;
        let mut shifts = vec![0; residual.len().saturating_sub(1)];
        for boundary in 0..shifts.len() {
            let (left, right) = (residual[boundary], residual[boundary + 1]);
            if left.signum() * right.signum() < 0 {
                let moved = left.abs().min(right.abs());
                shifts[boundary] = left.signum() * moved;
                residual[boundary] -= left.signum() * moved;
                residual[boundary + 1] -= right.signum() * moved;
            }
        }
        let longer = residual.iter().any(|&change| change > 0);
        let shorter = residual.iter().any(|&change| change < 0);
        (!(longer && shorter)).then_some(Self { shifts, residual })
    }

    /// Reference positions at the 5' end of `run` the previous read run reads.
    fn head(&self, run: usize) -> usize {
        run.checked_sub(1)
            .map_or(0, |boundary| self.shifts[boundary].max(0).unsigned_abs())
    }

    /// Reference positions at the 3' end of `run` the next read run reads.
    fn tail(&self, run: usize) -> usize {
        self.shifts
            .get(run)
            .map_or(0, |&shift| shift.min(0).unsigned_abs())
    }

    /// Reference positions of the next run that read run `run` reads.
    fn reach(&self, run: usize) -> usize {
        self.shifts
            .get(run)
            .map_or(0, |&shift| shift.max(0).unsigned_abs())
    }
}

/// The stretch as its moved run boundaries and one length edit at the 3' end
/// of each changed run, or `None` when it does not qualify or already reads
/// that way.
fn rebuild(
    inside: &[RawColumn],
    stretch: Stretch,
    reference: &Reference,
    masked: &dyn Fn(usize) -> bool,
) -> Option<Vec<RawColumn>> {
    let calls = inside
        .iter()
        .filter_map(|column| column.query_index.map(|index| (index, column.query_base)))
        .collect::<Vec<_>>();
    // In a stretch whose read spells the reference's runs, a substitution is
    // a run-length change the affine optimum merged into a neighbouring gap.
    let substituted = inside.iter().any(|column| {
        column.query_index.is_some()
            && column.reference_index.is_some()
            && column.query_base != column.reference_base
    });
    if !calls
        .iter()
        .all(|(_, call)| matches!(call, 'A' | 'C' | 'G' | 'T'))
        || !(substituted || calls.iter().any(|(index, _)| masked(*index)))
    {
        return None;
    }
    let read = runs(calls.iter().map(|(_, call)| *call));
    let expected =
        runs((stretch.first..=stretch.last).map(|index| char::from(base(reference, index))));
    if read.len() != expected.len()
        || read
            .iter()
            .zip(&expected)
            .any(|(observed, wanted)| observed.0 != wanted.0)
    {
        return None;
    }
    let lengths = |runs: &[(char, usize)]| runs.iter().map(|run| run.1).collect::<Vec<_>>();
    let edits = Edits::of(&lengths(&read), &lengths(&expected))?;
    let symbols = read.iter().map(|run| run.0).collect::<Vec<_>>();
    let mut rebuilt = Vec::with_capacity(inside.len());
    let mut call = calls.iter().map(|(index, _)| *index);
    let mut position = stretch.first;
    let mut push = |query: Option<char>, wanted: Option<char>| {
        let query_index = query.and_then(|_| call.next());
        let reference_index = wanted.map(|_| {
            position += 1;
            position - 1
        });
        rebuilt.push(RawColumn {
            query_base: query.unwrap_or('-'),
            reference_base: wanted.unwrap_or('-'),
            query_index,
            reference_index,
        });
    };
    for (run, &symbol) in symbols.iter().enumerate() {
        let wanted = expected[run].1;
        let (head, tail) = (edits.head(run), edits.tail(run));
        let residual = edits.residual[run];
        if head > 0 {
            // The previous read run reaches into this one, then takes its
            // own insertion at its 3' end.
            let previous = symbols[run - 1];
            (0..head).for_each(|_| push(Some(previous), Some(symbol)));
            (0..edits.residual[run - 1].max(0)).for_each(|_| push(Some(previous), None));
        }
        let core = wanted - head - tail;
        let deleted = residual.min(0).unsigned_abs();
        (0..core - deleted).for_each(|_| push(Some(symbol), Some(symbol)));
        (0..deleted).for_each(|_| push(None, Some(symbol)));
        if tail > 0 {
            let next = symbols[run + 1];
            (0..tail).for_each(|_| push(Some(next), Some(symbol)));
        } else if edits.reach(run) == 0 {
            (0..residual.max(0)).for_each(|_| push(Some(symbol), None));
        }
    }
    let same = rebuilt.len() == inside.len()
        && rebuilt.iter().zip(inside).all(|(new, old)| {
            new.query_index == old.query_index && new.reference_index == old.reference_index
        });
    (!same).then_some(rebuilt)
}

// SRS-ALN-016
#[cfg(test)]
mod tests {
    use dna_kernel::model::reference::ReferenceTopology;

    use super::*;

    fn reference(sequence: &str) -> Reference {
        Reference {
            name: "ref".into(),
            sequence: sequence.into(),
            topology: ReferenceTopology::Linear,
            sequence_sha256: String::new(),
        }
    }

    /// An alignment from equal-length gapped reference and query rows.
    fn alignment(reference_row: &str, query_row: &str) -> RawAlignment {
        let (mut reference_index, mut query_index) = (0, 0);
        let columns = reference_row
            .chars()
            .zip(query_row.chars())
            .map(|(reference_base, query_base)| {
                let column = RawColumn {
                    query_base,
                    reference_base,
                    query_index: (query_base != '-').then_some(query_index),
                    reference_index: (reference_base != '-').then_some(reference_index),
                };
                query_index += usize::from(query_base != '-');
                reference_index += usize::from(reference_base != '-');
                column
            })
            .collect::<Vec<_>>();
        RawAlignment {
            score: 0,
            start_reference: 0,
            end_reference: reference_index,
            metrics: traceback::metrics(&columns),
            columns,
        }
    }

    fn rows(alignment: &RawAlignment) -> (String, String) {
        alignment
            .columns
            .iter()
            .map(|column| (column.reference_base, column.query_base))
            .unzip()
    }

    fn reexpressed(reference_row: &str, query_row: &str, masked: &[usize]) -> RawAlignment {
        let mut raw = alignment(reference_row, query_row);
        let sequence = reference_row.replace('-', "");
        reexpress(
            &mut raw,
            &reference(&sequence),
            &|query| masked.contains(&query),
            64,
        );
        raw
    }

    // The reference rows read `GA CCCCCCC T CCCCC GTA`: two C runs around
    // an anchor T.

    #[test]
    fn rewrites_a_split_run_change_as_one_insertion_per_run() {
        // C8 T C6: the 8th C sits on T and the masked T is inserted.
        let raw = reexpressed("GACCCCCCCT-CCCCC-GTA", "GACCCCCCCCTCCCCCCGTA", &[10]);
        assert_eq!(
            rows(&raw),
            (
                "GACCCCCCC-TCCCCC-GTA".to_owned(),
                "GACCCCCCCCTCCCCCCGTA".to_owned()
            )
        );
        assert_eq!(raw.metrics.mismatches, 0);
        assert_eq!(raw.metrics.gap_opens, 2);
        let queries = raw
            .columns
            .iter()
            .filter_map(|column| column.query_index)
            .collect::<Vec<_>>();
        assert_eq!(queries, (0..20).collect::<Vec<_>>());
    }

    #[test]
    fn rewrites_shorter_runs_as_deletions_at_their_three_prime_ends() {
        // C6 T C4: the masked T sits on the 7th C and T310 is deleted.
        let raw = reexpressed("GACCCCCCCTCCCCCGTA", "GACCCCCCT-CCCC-GTA", &[8]);
        assert_eq!(
            rows(&raw),
            (
                "GACCCCCCCTCCCCCGTA".to_owned(),
                "GACCCCCC-TCCCC-GTA".to_owned()
            )
        );
    }

    #[test]
    fn rewrites_an_unmasked_merged_gap_beside_a_substitution() {
        // A clean C9 T C6: one gap opening saved for a substitution, the 8th
        // C on T, the 9th on C, and `TCC` inserted. The read shows the T, so the T
        // stays matched.
        let raw = reexpressed("GACCCCCCCTC---CCCCGTA", "GACCCCCCCCCTCCCCCCGTA", &[]);
        assert_eq!(
            rows(&raw),
            (
                "GACCCCCCC--TCCCCC-GTA".to_owned(),
                "GACCCCCCCCCTCCCCCCGTA".to_owned()
            )
        );
        assert_eq!(raw.metrics.mismatches, 0);
        assert_eq!(raw.metrics.gap_opens, 2);
    }

    #[test]
    fn keeps_an_unmasked_stretch_without_a_substitution() {
        // One merged two-base gap without a substitution already reads as
        // length changes only; a masked call alone re-expresses it.
        let raw = reexpressed("GACC--TGA", "GACCCTTGA", &[]);
        assert_eq!(rows(&raw), ("GACC--TGA".to_owned(), "GACCCTTGA".to_owned()));
        let masked = reexpressed("GACC--TGA", "GACCCTTGA", &[4]);
        assert_eq!(
            rows(&masked),
            ("GACC-T-GA".to_owned(), "GACCCTTGA".to_owned())
        );
    }

    #[test]
    fn moves_a_run_boundary_by_substitution_then_edits_lengths() {
        for masked in [&[4][..], &[]] {
            // A2 C8 against A4 C5: the C run starts two bases early, then
            // gains one base at its 3' end.
            let shifted = reexpressed("GAAAACC-CCCTG", "GAACCCCCCCCTG", masked);
            assert_eq!(
                rows(&shifted),
                ("GAAAACCCCC-TG".to_owned(), "GAACCCCCCCCTG".to_owned())
            );
        }
        // A read whose last in-phase call reads the C run one base early
        // (A2 C10 T C6 against A3 C7 T C5): the boundary moves by one
        // substitution and each C run gains its bases at its 3' end, so the
        // T stays matched.
        let raw = reexpressed(
            "GCACCAAACCCCCCCTC---CCCCGCT",
            "GCACCAACCCCCCCCCCTCCCCCCGCT",
            &[],
        );
        assert_eq!(
            rows(&raw),
            (
                "GCACCAAACCCCCCC--TCCCCC-GCT".to_owned(),
                "GCACCAACCCCCCCCCCTCCCCCCGCT".to_owned()
            )
        );
        assert_eq!(raw.metrics.mismatches, 1);
    }

    #[test]
    fn keeps_a_changed_run_structure_or_unpaired_mixed_changes() {
        for masked in [&[12][..], &[]] {
            // A substitution inside the second run gives the read another run.
            let substituted = reexpressed("GACCCCCCCTCC-CCCGTA", "GACCCCCCCTCCGCCCGTA", masked);
            assert_eq!(substituted.metrics.gap_opens, 1);
            assert_eq!(rows(&substituted).0, "GACCCCCCCTCC-CCCGTA");
        }
        for masked in [&[3][..], &[]] {
            // A shorter and T longer with C unchanged between them: no
            // boundary pairs the two changes, so the substitutions stay.
            let shifted = reexpressed("GAAACCCTTTG", "GAACCCTTTTG", masked);
            assert_eq!(
                rows(&shifted),
                ("GAAACCCTTTG".to_owned(), "GAACCCTTTTG".to_owned())
            );
        }
    }

    #[test]
    fn needs_an_anchor_on_both_sides() {
        // The read starts inside the first run, so nothing bounds it.
        let raw = reexpressed("CCCCCT-CCCCC-GTA", "CCCCCCTCCCCCCGTA", &[6]);
        assert_eq!(rows(&raw).0, "CCCCCT-CCCCC-GTA");
    }
}
