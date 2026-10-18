//! Bounded semi-global Gotoh dynamic programming.

use std::ops::Range;

use crate::canonical;
use crate::scoring::{
    NEGATIVE_INFINITY, State, add, scaled, substitution_index, substitution_scores,
};
use crate::traceback::{RawAlignment, TracebackInput, decode};
use crate::{AlignmentConfig, MAX_ALIGNMENT_CELLS};
use dna_kernel::error::{AlignmentError, Result};
use dna_kernel::read_evidence::EvidenceProfile;

/// Returns up to two distinct equally scoring placements.
pub(crate) fn align(
    query: &str,
    profiles: &[Option<EvidenceProfile>],
    reference: &str,
    config: &AlignmentConfig,
    modulo_length: Option<usize>,
) -> Result<Vec<RawAlignment>> {
    align_window(
        query,
        profiles,
        reference,
        0..reference.len(),
        config,
        modulo_length,
    )
}

/// Returns up to two distinct equally scoring placements inside one
/// reference window, in whole-reference coordinates.
///
/// Placements start and end inside the window. Canonicalization and the
/// circular span bound see whole-reference coordinates, so a window gives the
/// same placements as the whole reference whenever it contains every
/// optimal one.
pub(crate) fn align_window(
    query: &str,
    profiles: &[Option<EvidenceProfile>],
    reference: &str,
    window: Range<usize>,
    config: &AlignmentConfig,
    modulo_length: Option<usize>,
) -> Result<Vec<RawAlignment>> {
    let reference = reference
        .get(window.clone())
        .ok_or(AlignmentError::Inconsistent(
            "alignment window is outside the reference",
        ))?;
    if query.is_empty() || reference.is_empty() {
        return Err(AlignmentError::EmptyInput.into());
    }
    let query_bytes = query.as_bytes();
    if profiles.len() != query_bytes.len() {
        return Err(AlignmentError::ProfileLengthMismatch {
            bases: query_bytes.len(),
            profiles: profiles.len(),
        }
        .into());
    }
    let reference_bytes = reference.as_bytes();
    let rows = query_bytes
        .len()
        .checked_add(1)
        .ok_or(AlignmentError::Overflow("query length overflow"))?;
    let width = reference_bytes
        .len()
        .checked_add(1)
        .ok_or(AlignmentError::Overflow("reference length overflow"))?;
    let cells = rows
        .checked_mul(width)
        .ok_or(AlignmentError::Overflow("alignment cell count overflow"))?;
    if cells > MAX_ALIGNMENT_CELLS {
        return Err(AlignmentError::TooManyCells {
            cells,
            maximum: MAX_ALIGNMENT_CELLS,
        }
        .into());
    }
    let mut trace = vec![0_u8; cells];
    let last_row = fill(profiles, reference_bytes, config, width, &mut trace);

    let mut endpoints = Vec::with_capacity(width);
    for (column, cell) in last_row.iter().enumerate() {
        let (score, state) = [
            (cell.matched, State::Match),
            (cell.deletion, State::Deletion),
            (cell.insertion, State::Insertion),
        ]
        .into_iter()
        .max_by_key(|(score, state)| (*score, state_priority(*state)))
        .ok_or(AlignmentError::Inconsistent(
            "alignment endpoint state is missing",
        ))?;
        endpoints.push((score, column, state));
    }
    endpoints
        .sort_unstable_by(|left, right| right.0.cmp(&left.0).then_with(|| left.1.cmp(&right.1)));

    let mut placements = Vec::new();
    let mut bounded_best_score = None;
    for (score, column, state) in endpoints {
        if bounded_best_score.is_some_and(|best| score < best) {
            break;
        }
        let mut raw = decode(TracebackInput {
            query: query_bytes,
            reference: reference_bytes,
            trace: &trace,
            row_width: width,
            endpoint: column,
            state,
            score,
        })?;
        raw.shift_reference(window.start);
        canonical::right_align(&mut raw, profiles, config, modulo_length)?;
        if let Some(length) = modulo_length
            && raw.end_reference - raw.start_reference > length
        {
            continue;
        }
        bounded_best_score.get_or_insert(score);
        let duplicate = placements
            .iter()
            .any(|existing: &RawAlignment| existing.same_placement(&raw, modulo_length));
        if !duplicate {
            placements.push(raw);
            if placements.len() == 2 {
                return Ok(placements);
            }
        }
    }
    if placements.is_empty() {
        return Err(AlignmentError::NoTraceback.into());
    }
    Ok(placements)
}

/// Scores of the three Gotoh states at one cell.
#[derive(Debug, Clone, Copy)]
struct Cell {
    matched: i64,
    insertion: i64,
    deletion: i64,
}

/// Fills the traceback matrix row by row and returns the last row's scores.
///
/// Each cell's packed traceback byte is written once: bits 0-1 hold the
/// diagonal predecessor state, bit 2 records insertion extension, and bit 3
/// records deletion extension. Equal scores keep the documented state order
/// Match > Deletion > Insertion and prefer extension over opening.
fn fill(
    profiles: &[Option<EvidenceProfile>],
    reference: &[u8],
    config: &AlignmentConfig,
    width: usize,
    trace: &mut [u8],
) -> Vec<Cell> {
    let reference_slots = reference
        .iter()
        .map(|base| substitution_index(*base))
        .collect::<Vec<_>>();
    let gap_extension = scaled(config.gap_extension_score);
    let open_and_extend = scaled(config.gap_open_score) + gap_extension;
    let unreachable = Cell {
        matched: NEGATIVE_INFINITY,
        insertion: NEGATIVE_INFINITY,
        deletion: NEGATIVE_INFINITY,
    };
    let mut previous = vec![
        Cell {
            matched: 0,
            ..unreachable
        };
        width
    ];
    let mut current = vec![unreachable; width];

    for (profile, trace_row) in profiles.iter().zip(trace.chunks_exact_mut(width).skip(1)) {
        let substitution = substitution_scores(*profile, config);
        let (first_trace, trace_cells) = trace_row.split_at_mut(1);
        let (insertion, extended) = gap(
            previous[0].matched,
            previous[0].insertion,
            open_and_extend,
            gap_extension,
        );
        current[0] = Cell {
            insertion,
            ..unreachable
        };
        first_trace[0] = if extended { 0b100 } else { 0 };

        let mut left = current[0];
        for (((diagonal_and_above, cell), slot), packed) in previous
            .windows(2)
            .zip(current.iter_mut().skip(1))
            .zip(&reference_slots)
            .zip(trace_cells.iter_mut())
        {
            let diagonal = diagonal_and_above[0];
            let above = diagonal_and_above[1];
            let (best_diagonal, predecessor) = if diagonal.matched >= diagonal.deletion
                && diagonal.matched >= diagonal.insertion
            {
                (diagonal.matched, State::Match)
            } else if diagonal.deletion >= diagonal.insertion {
                (diagonal.deletion, State::Deletion)
            } else {
                (diagonal.insertion, State::Insertion)
            };
            let matched = add(best_diagonal, substitution[*slot]);
            let (insertion, insertion_extended) = gap(
                above.matched,
                above.insertion,
                open_and_extend,
                gap_extension,
            );
            let (deletion, deletion_extended) =
                gap(left.matched, left.deletion, open_and_extend, gap_extension);
            *cell = Cell {
                matched,
                insertion,
                deletion,
            };
            left = *cell;
            *packed = predecessor as u8
                | if insertion_extended { 0b100 } else { 0 }
                | if deletion_extended { 0b1000 } else { 0 };
        }
        std::mem::swap(&mut previous, &mut current);
    }
    previous
}

/// Scores entering a gap state by extension or by opening from Match,
/// preferring extension on equal scores.
fn gap(matched: i64, gapped: i64, open_and_extend: i64, gap_extension: i64) -> (i64, bool) {
    let open = add(matched, open_and_extend);
    let extend = add(gapped, gap_extension);
    if extend >= open {
        (extend, true)
    } else {
        (open, false)
    }
}

const fn state_priority(state: State) -> u8 {
    match state {
        State::Match => 2,
        State::Deletion => 1,
        State::Insertion => 0,
    }
}

#[cfg(test)]
mod tests {
    use crate::scoring::SCORE_SCALE;

    use super::*;

    type TestResult<T = ()> = std::result::Result<T, Box<dyn std::error::Error>>;

    fn config() -> AlignmentConfig {
        AlignmentConfig {
            match_score: 3,
            mismatch_score: -5,
            ambiguous_score: 0,
            gap_open_score: -10,
            gap_extension_score: -4,
            minimum_callable_bases: 1,
            minimum_identity: 0.8,
        }
    }

    fn profiles(query: &str) -> Vec<Option<EvidenceProfile>> {
        query
            .bytes()
            .map(|base| {
                let weights = match base {
                    b'A' => [1.0, 0.0, 0.0, 0.0],
                    b'C' => [0.0, 1.0, 0.0, 0.0],
                    b'G' => [0.0, 0.0, 1.0, 0.0],
                    b'T' => [0.0, 0.0, 0.0, 1.0],
                    _ => return None,
                };
                Some(EvidenceProfile { weights })
            })
            .collect()
    }

    #[test]
    fn permits_free_reference_flanks() -> Result<()> {
        let query = "ACGT";
        let alignments = align(query, &profiles(query), "TTACGTGG", &config(), None)?;
        assert_eq!(alignments[0].score, 12 * SCORE_SCALE);
        assert_eq!(alignments[0].start_reference, 2);
        assert_eq!(alignments[0].end_reference, 6);
        Ok(())
    }

    #[test]
    fn accepts_one_full_circular_reference_span() -> Result<()> {
        let reference = "ACGTCAGTACGATCGTACCTGAGTACGA";
        let query = format!("{}{}", &reference[18..], &reference[..18]);
        let working_reference = format!("{reference}{reference}");
        let alignments = align(
            &query,
            &profiles(&query),
            &working_reference,
            &config(),
            Some(reference.len()),
        )?;
        assert_eq!(alignments[0].score, 84 * SCORE_SCALE);
        assert_eq!(
            alignments[0].end_reference - alignments[0].start_reference,
            reference.len()
        );
        Ok(())
    }

    #[test]
    fn scores_one_base_gap_as_open_plus_extension() -> Result<()> {
        let query = "ACGTT";
        let alignments = align(query, &profiles(query), "ACGT", &config(), None)?;
        assert_eq!(alignments[0].metrics.gap_opens, 1);
        assert_eq!(alignments[0].score, -2 * SCORE_SCALE);
        Ok(())
    }

    #[test]
    fn profile_can_score_unresolved_query_character() -> Result<()> {
        let profile = EvidenceProfile {
            weights: [1.0, 0.0, 0.0, 0.0],
        };
        let alignments = align("N", &[Some(profile)], "A", &config(), None)?;
        assert_eq!(alignments[0].score, 3 * SCORE_SCALE);
        assert_eq!(alignments[0].metrics.unresolved_query_bases, 1);
        Ok(())
    }

    #[test]
    fn rejects_query_profile_length_mismatch() {
        assert!(align("AC", &[None], "AC", &config(), None).is_err());
    }

    #[test]
    fn canonicalizes_homopolymer_deletion_to_rightmost_reference_base() -> TestResult {
        let query = "GCCAAAGTT";
        let alignments = align(query, &profiles(query), "GCCAAAAGTT", &config(), None)?;
        assert_eq!(alignments.len(), 1);
        let deletion = alignments[0]
            .columns
            .iter()
            .find(|column| column.query_base == '-')
            .ok_or("expected canonical deletion")?;
        assert_eq!(deletion.reference_index, Some(6));
        assert_eq!(deletion.reference_base, 'A');
        Ok(())
    }

    #[test]
    fn canonicalizes_homopolymer_insertion_to_rightmost_boundary() -> TestResult {
        let query = "CAAAAAG";
        let alignments = align(query, &profiles(query), "CAAAAG", &config(), None)?;
        assert_eq!(alignments.len(), 1);
        let insertion_index = alignments[0]
            .columns
            .iter()
            .position(|column| column.reference_base == '-')
            .ok_or("expected canonical insertion")?;
        assert_eq!(
            alignments[0].columns[insertion_index - 1].reference_index,
            Some(4)
        );
        assert_eq!(
            alignments[0].columns[insertion_index + 1].reference_index,
            Some(5)
        );
        Ok(())
    }

    #[test]
    fn preserves_distinct_repeat_placement_ambiguity_without_an_indel() -> Result<()> {
        let query = "AAA";
        let alignments = align(query, &profiles(query), "AAAAA", &config(), None)?;
        assert_eq!(alignments.len(), 2);
        assert_ne!(alignments[0].start_reference, alignments[1].start_reference);
        Ok(())
    }

    #[test]
    fn preserves_scores_beyond_i32_range() -> Result<()> {
        let mut scoring = config();
        scoring.match_score = i32::MAX;
        let query = "AA";
        let alignments = align(query, &profiles(query), "AA", &scoring, None)?;
        assert_eq!(alignments[0].score, 2 * i64::from(i32::MAX) * SCORE_SCALE);
        Ok(())
    }
}
