//! Run-structure re-expression of masked repeat stretches (ADR-0072).
//!
//! Masked calls keep anchoring the optimal alignment, but their mixed signal
//! can make it split one run-length change into a substitution next to an
//! insertion, or merge two of them into one gap. Where the read's calls spell
//! the reference's runs in the same order, every run only longer or every run
//! only shorter, and a masked call lies among them, the selected placement is
//! re-expressed as one length edit at the 3' end of each changed run. The
//! placement, its score, and every column outside such a stretch are kept.

use crate::alignment::traceback::{self, RawAlignment, RawColumn};
use dna_kernel::model::reference::Reference;

/// Longest reference stretch, in bases, a re-expression may span; it bounds
/// the work per disturbed block and keeps a stretch local to one repeat.
const MAX_STRETCH: usize = 64;

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
/// tells whether the call at an oriented query index is masked.
pub(super) fn reexpress(
    alignment: &mut RawAlignment,
    reference: &Reference,
    masked: &dyn Fn(usize) -> bool,
) {
    let stretches = stretches(&alignment.columns, reference);
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
fn stretches(columns: &[RawColumn], reference: &Reference) -> Vec<Stretch> {
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
        let Some(stretch) = anchored(columns, reference, &column_of, low, high) else {
            continue;
        };
        match found.last_mut() {
            Some(previous) if stretch.left < previous.right => {
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
) -> Option<Stretch> {
    let (mut first, mut last) = (run_start(reference, low), run_end(reference, high));
    while last - first < MAX_STRETCH && same_copy(reference, first, last) {
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

/// The stretch as one length edit at the 3' end of each changed run, or
/// `None` when it does not qualify or already reads that way.
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
    if !calls
        .iter()
        .all(|(_, call)| matches!(call, 'A' | 'C' | 'G' | 'T'))
        || !calls.iter().any(|(index, _)| masked(*index))
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
    let longer = read
        .iter()
        .zip(&expected)
        .any(|(observed, wanted)| observed.1 > wanted.1);
    let shorter = read
        .iter()
        .zip(&expected)
        .any(|(observed, wanted)| observed.1 < wanted.1);
    if longer == shorter {
        return None;
    }
    let mut rebuilt = Vec::with_capacity(inside.len());
    let mut call = calls.iter();
    let mut position = stretch.first;
    for (&(symbol, observed), &(_, wanted)) in read.iter().zip(&expected) {
        for paired in 0..observed.max(wanted) {
            let query_index = (paired < observed)
                .then(|| call.next().map(|(index, _)| *index))
                .flatten();
            let reference_index = (paired < wanted).then_some(position);
            if reference_index.is_some() {
                position += 1;
            }
            rebuilt.push(RawColumn {
                query_base: if query_index.is_some() { symbol } else { '-' },
                reference_base: if reference_index.is_some() {
                    symbol
                } else {
                    '-'
                },
                query_index,
                reference_index,
            });
        }
    }
    let same = rebuilt.len() == inside.len()
        && rebuilt.iter().zip(inside).all(|(new, old)| {
            new.query_index == old.query_index && new.reference_index == old.reference_index
        });
    (!same).then_some(rebuilt)
}

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
        reexpress(&mut raw, &reference(&sequence), &|query| {
            masked.contains(&query)
        });
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
    fn keeps_stretches_without_a_masked_call() {
        let raw = reexpressed("GACCCCCCCT-CCCCC-GTA", "GACCCCCCCCTCCCCCCGTA", &[]);
        assert_eq!(
            rows(&raw),
            (
                "GACCCCCCCT-CCCCC-GTA".to_owned(),
                "GACCCCCCCCTCCCCCCGTA".to_owned()
            )
        );
    }

    #[test]
    fn keeps_a_changed_run_structure_or_mixed_length_changes() {
        // A substitution inside the second run gives the read another run.
        let substituted = reexpressed("GACCCCCCCTCC-CCCGTA", "GACCCCCCCTCCGCCCGTA", &[12]);
        assert_eq!(substituted.metrics.gap_opens, 1);
        assert_eq!(rows(&substituted).0, "GACCCCCCCTCC-CCCGTA");
        // One run shorter and the next longer: substitutions explain it too.
        let shifted = reexpressed("GAAAACC-CCCTG", "GAACCCCCCCCTG", &[4]);
        assert_eq!(rows(&shifted).0, "GAAAACC-CCCTG");
    }

    #[test]
    fn needs_an_anchor_on_both_sides() {
        // The read starts inside the first run, so nothing bounds it.
        let raw = reexpressed("CCCCCT-CCCCC-GTA", "CCCCCCTCCCCCCGTA", &[6]);
        assert_eq!(rows(&raw).0, "CCCCCT-CCCCC-GTA");
    }
}
