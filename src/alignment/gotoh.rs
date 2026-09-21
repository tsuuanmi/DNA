//! Bounded semi-global Gotoh dynamic programming.

use crate::alignment::canonical;
use crate::alignment::scoring::{NEGATIVE_INFINITY, State, add, scaled, substitution};
use crate::alignment::traceback::{RawAlignment, TracebackInput, decode};
use crate::config::{AlignmentConfig, MAX_ALIGNMENT_CELLS};
use crate::error::{Error, Result};
use crate::model::locus_evidence::EvidenceProfile;

/// Returns up to two distinct equally scoring placements.
pub(crate) fn align(
    query: &str,
    profiles: &[Option<EvidenceProfile>],
    reference: &str,
    config: &AlignmentConfig,
    modulo_length: Option<usize>,
) -> Result<Vec<RawAlignment>> {
    if query.is_empty() || reference.is_empty() {
        return Err(Error::Alignment(
            "query and reference must both be non-empty".into(),
        ));
    }
    let query_bytes = query.as_bytes();
    if profiles.len() != query_bytes.len() {
        return Err(Error::Alignment(format!(
            "query/profile length mismatch: {} query bases, {} profiles",
            query_bytes.len(),
            profiles.len()
        )));
    }
    let reference_bytes = reference.as_bytes();
    let rows = query_bytes
        .len()
        .checked_add(1)
        .ok_or_else(|| Error::Alignment("query length overflow".into()))?;
    let width = reference_bytes
        .len()
        .checked_add(1)
        .ok_or_else(|| Error::Alignment("reference length overflow".into()))?;
    let cells = rows
        .checked_mul(width)
        .ok_or_else(|| Error::Alignment("alignment cell count overflow".into()))?;
    if cells > MAX_ALIGNMENT_CELLS {
        return Err(Error::Alignment(format!(
            "alignment requires {cells} cells; cap is {MAX_ALIGNMENT_CELLS}"
        )));
    }
    let mut trace = vec![0_u8; cells];
    let mut previous_match = vec![0_i64; width];
    let mut previous_insertion = vec![NEGATIVE_INFINITY; width];
    let mut previous_deletion = vec![NEGATIVE_INFINITY; width];
    let mut current_match = vec![NEGATIVE_INFINITY; width];
    let mut current_insertion = vec![NEGATIVE_INFINITY; width];
    let mut current_deletion = vec![NEGATIVE_INFINITY; width];
    let gap_extension = scaled(config.gap_extension_score);
    let open_and_extend = scaled(config.gap_open_score) + gap_extension;

    for row in 1..rows {
        current_match.fill(NEGATIVE_INFINITY);
        current_insertion.fill(NEGATIVE_INFINITY);
        current_deletion.fill(NEGATIVE_INFINITY);
        let open = add(previous_match[0], open_and_extend);
        let extend = add(previous_insertion[0], gap_extension);
        if extend >= open {
            current_insertion[0] = extend;
            trace[row * width] |= 0b100;
        } else {
            current_insertion[0] = open;
        }

        for column in 1..width {
            let diagonal = [
                (previous_match[column - 1], State::Match),
                (previous_deletion[column - 1], State::Deletion),
                (previous_insertion[column - 1], State::Insertion),
            ];
            let (best_diagonal, predecessor) = diagonal
                .into_iter()
                .max_by_key(|(score, state)| (*score, state_priority(*state)))
                .ok_or_else(|| Error::Alignment("missing diagonal state".into()))?;
            current_match[column] = add(
                best_diagonal,
                substitution(profiles[row - 1], reference_bytes[column - 1], config),
            );
            trace[row * width + column] |= predecessor as u8;

            let open = add(previous_match[column], open_and_extend);
            let extend = add(previous_insertion[column], gap_extension);
            if extend >= open {
                current_insertion[column] = extend;
                trace[row * width + column] |= 0b100;
            } else {
                current_insertion[column] = open;
            }

            let open = add(current_match[column - 1], open_and_extend);
            let extend = add(current_deletion[column - 1], gap_extension);
            if extend >= open {
                current_deletion[column] = extend;
                trace[row * width + column] |= 0b1000;
            } else {
                current_deletion[column] = open;
            }
        }
        std::mem::swap(&mut previous_match, &mut current_match);
        std::mem::swap(&mut previous_insertion, &mut current_insertion);
        std::mem::swap(&mut previous_deletion, &mut current_deletion);
    }

    let mut endpoints = Vec::with_capacity(width);
    for column in 0..width {
        let (score, state) = [
            (previous_match[column], State::Match),
            (previous_deletion[column], State::Deletion),
            (previous_insertion[column], State::Insertion),
        ]
        .into_iter()
        .max_by_key(|(score, state)| (*score, state_priority(*state)))
        .ok_or_else(|| Error::Alignment("alignment endpoint state is missing".into()))?;
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
        canonical::right_align(&mut raw, profiles, config, modulo_length)?;
        if let Some(length) = modulo_length
            && raw.end_reference - raw.start_reference > length
        {
            continue;
        }
        bounded_best_score.get_or_insert(score);
        let key_start = modulo_length
            .map(|length| raw.start_reference % length)
            .unwrap_or(raw.start_reference);
        let duplicate = placements.iter().any(|existing: &RawAlignment| {
            let existing_start = modulo_length
                .map(|length| existing.start_reference % length)
                .unwrap_or(existing.start_reference);
            existing_start == key_start
                && existing.columns.len() == raw.columns.len()
                && existing
                    .columns
                    .iter()
                    .zip(&raw.columns)
                    .all(|(left, right)| {
                        left.query_base == right.query_base
                            && left.reference_base == right.reference_base
                    })
        });
        if !duplicate {
            placements.push(raw);
            if placements.len() == 2 {
                return Ok(placements);
            }
        }
    }
    if placements.is_empty() {
        return Err(Error::Alignment(
            "no valid bounded alignment traceback was found".into(),
        ));
    }
    Ok(placements)
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
    use crate::alignment::scoring::SCORE_SCALE;

    use super::*;

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
    fn canonicalizes_homopolymer_deletion_to_rightmost_reference_base() -> Result<()> {
        let query = "GCCAAAGTT";
        let alignments = align(query, &profiles(query), "GCCAAAAGTT", &config(), None)?;
        assert_eq!(alignments.len(), 1);
        let deletion = alignments[0]
            .columns
            .iter()
            .find(|column| column.query_base == '-')
            .ok_or_else(|| Error::Alignment("expected canonical deletion".into()))?;
        assert_eq!(deletion.reference_index, Some(6));
        assert_eq!(deletion.reference_base, 'A');
        Ok(())
    }

    #[test]
    fn canonicalizes_homopolymer_insertion_to_rightmost_boundary() -> Result<()> {
        let query = "CAAAAAG";
        let alignments = align(query, &profiles(query), "CAAAAG", &config(), None)?;
        assert_eq!(alignments.len(), 1);
        let insertion_index = alignments[0]
            .columns
            .iter()
            .position(|column| column.reference_base == '-')
            .ok_or_else(|| Error::Alignment("expected canonical insertion".into()))?;
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
