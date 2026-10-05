//! Provably exact alignment acceleration.

use memchr::memmem::Finder;

use crate::alignment::gotoh;
use crate::alignment::scoring::{scaled, substitution_scores};
use crate::alignment::traceback::{RawAlignment, RawColumn, metrics};
use crate::config::{AlignmentConfig, MAX_ALIGNMENT_CELLS};
use crate::error::Result;
use crate::model::locus_evidence::EvidenceProfile;

const CANONICAL_BASES: &[u8; 4] = b"ACGT";
const PROBE_COUNT: usize = 4;
const PROBE_MAX_SEED_LENGTH: usize = 16;
const MAX_SEED_HITS: usize = 128;
const MAX_CERTIFIED_SEEDS: usize = 32;

#[derive(Debug)]
pub(crate) enum UpperBoundPlacement {
    Unproven,
    Unattained,
    Unique(RawAlignment),
    Ambiguous,
}

#[derive(Debug)]
struct ProfileBound {
    optimal_sequence: Vec<u8>,
    upper_bound: i64,
    edit_loss_floor: i64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Window {
    start: usize,
    end: usize,
}

/// Classifies whether one orientation can attain its theoretical profile-score upper bound.
///
/// The proof is deliberately conservative. Every retained locus must have one strictly
/// best canonical reference base, that substitution must beat extending a query gap, and
/// the ordinary Gotoh problem must remain within the configured matrix cap. Under those
/// conditions, the concatenated per-locus best bases define the only possible gapless
/// sequence that can attain the upper bound. Exact occurrences of that sequence therefore
/// identify every maximum-score placement for this orientation.
pub(crate) fn classify(
    query: &str,
    profiles: &[Option<EvidenceProfile>],
    reference: &str,
    config: &AlignmentConfig,
    modulo_length: Option<usize>,
) -> UpperBoundPlacement {
    if !problem_is_provable(query, profiles, reference, config, modulo_length) {
        return UpperBoundPlacement::Unproven;
    }
    let Some(bound) = profile_bound(profiles, config) else {
        return UpperBoundPlacement::Unproven;
    };
    let Some(limit) =
        occurrence_start_limit(reference, bound.optimal_sequence.len(), modulo_length)
    else {
        return UpperBoundPlacement::Unproven;
    };

    let haystack_end = limit + bound.optimal_sequence.len() - 1;
    let haystack = &reference.as_bytes()[..haystack_end];
    let finder = Finder::new(&bound.optimal_sequence);
    let mut search_start = 0;
    let mut found = None;
    while search_start < limit {
        let Some(relative) = finder.find(&haystack[search_start..]) else {
            break;
        };
        let start = search_start + relative;
        if start >= limit {
            break;
        }
        if found.is_some() {
            return UpperBoundPlacement::Ambiguous;
        }
        found = Some(start);
        search_start = start + 1;
    }

    let Some(start_reference) = found else {
        return UpperBoundPlacement::Unattained;
    };
    UpperBoundPlacement::Unique(gapless_alignment(
        query,
        &bound.optimal_sequence,
        start_reference,
        bound.upper_bound,
    ))
}

/// Attempts an exact seeded-pruning alignment for one orientation.
///
/// A cheap seed-derived local alignment supplies only a valid lower-bound score.
/// That score is converted to a maximum edit budget from the profile upper bound.
/// The q-gram lemma then proves that every alignment able to tie or beat the lower
/// bound contains at least one exact seed, allowing Gotoh to run only on the union
/// of certified windows. Returning `None` means the proof or cost bound was not
/// strong enough and the caller must use full-reference Gotoh.
pub(crate) fn align_pruned(
    query: &str,
    profiles: &[Option<EvidenceProfile>],
    reference: &str,
    config: &AlignmentConfig,
    modulo_length: Option<usize>,
) -> Result<Option<Vec<RawAlignment>>> {
    if !problem_is_provable(query, profiles, reference, config, modulo_length) {
        return Ok(None);
    }
    let Some(bound) = profile_bound(profiles, config) else {
        return Ok(None);
    };
    let Some(lower_bound) =
        seed_lower_bound(query, profiles, reference, config, modulo_length, &bound)?
    else {
        return Ok(None);
    };
    let result = align_at_or_above_with_bound(
        query,
        profiles,
        reference,
        config,
        modulo_length,
        lower_bound,
        &bound,
    )?;
    Ok(result.filter(|placements| !placements.is_empty()))
}

/// Proves and evaluates every placement in one orientation that can reach `threshold`.
///
/// `Some(empty)` is a proof that this orientation cannot reach the threshold.
/// `None` means the proof is too weak or expensive, so the caller must use Gotoh.
pub(crate) fn align_at_or_above(
    query: &str,
    profiles: &[Option<EvidenceProfile>],
    reference: &str,
    config: &AlignmentConfig,
    modulo_length: Option<usize>,
    threshold: i64,
) -> Result<Option<Vec<RawAlignment>>> {
    if !problem_is_provable(query, profiles, reference, config, modulo_length) {
        return Ok(None);
    }
    let Some(bound) = profile_bound(profiles, config) else {
        return Ok(None);
    };
    align_at_or_above_with_bound(
        query,
        profiles,
        reference,
        config,
        modulo_length,
        threshold,
        &bound,
    )
}

fn align_at_or_above_with_bound(
    query: &str,
    profiles: &[Option<EvidenceProfile>],
    reference: &str,
    config: &AlignmentConfig,
    modulo_length: Option<usize>,
    threshold: i64,
    bound: &ProfileBound,
) -> Result<Option<Vec<RawAlignment>>> {
    if threshold > bound.upper_bound {
        return Ok(Some(Vec::new()));
    }
    let Some(deficit) = bound.upper_bound.checked_sub(threshold) else {
        return Ok(None);
    };
    let edit_budget_i64 = deficit / bound.edit_loss_floor;
    let Ok(edit_budget) = usize::try_from(edit_budget_i64) else {
        return Ok(None);
    };
    if edit_budget >= query.len() {
        return Ok(None);
    }
    let seed_count = edit_budget + 1;
    if seed_count > MAX_CERTIFIED_SEEDS {
        return Ok(None);
    }

    let Some(windows) = certified_windows(
        &bound.optimal_sequence,
        reference,
        modulo_length,
        edit_budget,
        seed_count,
    ) else {
        return Ok(None);
    };
    if windows.is_empty() {
        return Ok(Some(Vec::new()));
    }
    if !windows_are_economical(
        &windows,
        canonical_reference_length(reference, modulo_length),
    ) {
        return Ok(None);
    }

    let placements = align_windows(query, profiles, reference, config, &windows)?;
    Ok(Some(best_placements_at_or_above(
        placements,
        threshold,
        modulo_length,
    )))
}

fn profile_bound(
    profiles: &[Option<EvidenceProfile>],
    config: &AlignmentConfig,
) -> Option<ProfileBound> {
    let gap_extension = scaled(config.gap_extension_score);
    if gap_extension >= 0 {
        return None;
    }

    let mut optimal_sequence = Vec::with_capacity(profiles.len());
    let mut upper_bound = 0_i64;
    let mut edit_loss_floor = -gap_extension;

    for profile in profiles {
        let scores = substitution_scores(*profile, config);
        let (best_index, best_score) = unique_canonical_maximum(scores)?;
        if best_score <= gap_extension {
            return None;
        }
        let second_best = scores
            .iter()
            .enumerate()
            .filter(|(index, _)| *index != best_index)
            .map(|(_, score)| *score)
            .max()?;
        let substitution_loss = best_score.checked_sub(second_best)?;
        let insertion_loss = best_score.checked_sub(gap_extension)?;
        if substitution_loss <= 0 || insertion_loss <= 0 {
            return None;
        }

        upper_bound = upper_bound.checked_add(best_score)?;
        edit_loss_floor = edit_loss_floor.min(substitution_loss).min(insertion_loss);
        optimal_sequence.push(CANONICAL_BASES[best_index]);
    }

    (edit_loss_floor > 0).then_some(ProfileBound {
        optimal_sequence,
        upper_bound,
        edit_loss_floor,
    })
}

fn seed_lower_bound(
    query: &str,
    profiles: &[Option<EvidenceProfile>],
    reference: &str,
    config: &AlignmentConfig,
    modulo_length: Option<usize>,
    bound: &ProfileBound,
) -> Result<Option<i64>> {
    if query.len() < PROBE_COUNT {
        return Ok(None);
    }
    let margin = 32.min((query.len() / 4).max(4));
    let probe_ranges = partition_ranges(query.len(), PROBE_COUNT);
    let mut windows = Vec::new();
    let mut hit_count = 0;

    for (block_start, block_end) in probe_ranges {
        let block_length = block_end - block_start;
        let seed_length = block_length.min(PROBE_MAX_SEED_LENGTH);
        let seed_start = block_start + (block_length - seed_length) / 2;
        let seed_end = seed_start + seed_length;
        let seed = &bound.optimal_sequence[seed_start..seed_end];

        for hit in seed_hits(reference, seed, modulo_length, false) {
            hit_count += 1;
            if hit_count > MAX_SEED_HITS {
                return Ok(None);
            }
            if let Some(window) = window_from_seed(
                hit,
                seed_start,
                seed_end,
                query.len(),
                margin,
                reference,
                modulo_length,
            ) {
                windows.push(window);
            }
        }
    }

    let windows = merge_windows(windows);
    if windows.is_empty() {
        return Ok(None);
    }

    let mut best = None;
    for window in windows {
        let slice = &reference[window.start..window.end];
        let placements = gotoh::align(query, profiles, slice, config, None)?;
        if let Some(first) = placements.first() {
            best = Some(best.map_or(first.score, |score: i64| score.max(first.score)));
        }
    }
    Ok(best)
}

fn certified_windows(
    optimal_sequence: &[u8],
    reference: &str,
    modulo_length: Option<usize>,
    edit_budget: usize,
    seed_count: usize,
) -> Option<Vec<Window>> {
    let mut windows = Vec::new();
    let mut hit_count = 0;
    for (seed_start, seed_end) in partition_ranges(optimal_sequence.len(), seed_count) {
        let seed = &optimal_sequence[seed_start..seed_end];
        let hits = seed_hits(reference, seed, modulo_length, true);
        for hit in hits {
            hit_count += 1;
            if hit_count > MAX_SEED_HITS {
                return None;
            }
            let window = window_from_seed(
                hit,
                seed_start,
                seed_end,
                optimal_sequence.len(),
                edit_budget,
                reference,
                modulo_length,
            )?;
            windows.push(window);
        }
    }
    Some(merge_windows(windows))
}

fn seed_hits(
    reference: &str,
    seed: &[u8],
    modulo_length: Option<usize>,
    include_seam: bool,
) -> Vec<usize> {
    if seed.is_empty() {
        return Vec::new();
    }
    let canonical_length = canonical_reference_length(reference, modulo_length);
    let haystack_end = if include_seam {
        modulo_length
            .and_then(|length| length.checked_add(seed.len().saturating_sub(1)))
            .map_or(reference.len(), |end| end.min(reference.len()))
    } else {
        canonical_length
    };
    let haystack = &reference.as_bytes()[..haystack_end];
    let finder = Finder::new(seed);
    let mut hits = Vec::new();
    let mut search_start = 0;

    while search_start <= haystack.len().saturating_sub(seed.len()) {
        let Some(relative) = finder.find(&haystack[search_start..]) else {
            break;
        };
        let hit = search_start + relative;
        if modulo_length.is_some_and(|length| hit >= length) {
            break;
        }
        hits.push(hit);
        search_start = hit + 1;
    }
    hits
}

fn window_from_seed(
    hit: usize,
    seed_start: usize,
    seed_end: usize,
    query_length: usize,
    drift: usize,
    reference: &str,
    modulo_length: Option<usize>,
) -> Option<Window> {
    let canonical_length = canonical_reference_length(reference, modulo_length);
    let prefix = seed_start.checked_add(drift)?;
    let suffix = query_length.checked_sub(seed_end)?.checked_add(drift)?;
    let start = hit.saturating_sub(prefix);
    let end = hit
        .checked_add(seed_end - seed_start)?
        .checked_add(suffix)?
        .min(canonical_length);

    if let Some(length) = modulo_length {
        let crosses_left = hit < prefix;
        let seed_crosses_seam = hit.checked_add(seed_end - seed_start)? > length;
        let crosses_right = hit
            .checked_add(seed_end - seed_start)?
            .checked_add(suffix)?
            > length;
        if crosses_left || seed_crosses_seam || crosses_right {
            return None;
        }
    }

    (start < end).then_some(Window { start, end })
}

fn align_windows(
    query: &str,
    profiles: &[Option<EvidenceProfile>],
    reference: &str,
    config: &AlignmentConfig,
    windows: &[Window],
) -> Result<Vec<RawAlignment>> {
    let mut placements = Vec::new();
    for window in windows {
        let slice = &reference[window.start..window.end];
        let local = gotoh::align(query, profiles, slice, config, None)?;
        placements.extend(local.into_iter().map(|mut alignment| {
            alignment.start_reference += window.start;
            alignment.end_reference += window.start;
            for column in &mut alignment.columns {
                if let Some(index) = &mut column.reference_index {
                    *index += window.start;
                }
            }
            alignment
        }));
    }
    Ok(placements)
}

fn best_placements_at_or_above(
    mut placements: Vec<RawAlignment>,
    threshold: i64,
    modulo_length: Option<usize>,
) -> Vec<RawAlignment> {
    let Some(best_score) = placements
        .iter()
        .filter(|alignment| alignment.score >= threshold)
        .map(|alignment| alignment.score)
        .max()
    else {
        return Vec::new();
    };
    placements.retain(|alignment| alignment.score == best_score);
    placements.sort_unstable_by_key(|alignment| {
        (
            modulo_length
                .map(|length| alignment.start_reference % length)
                .unwrap_or(alignment.start_reference),
            alignment.end_reference,
        )
    });

    let mut unique = Vec::new();
    for alignment in placements {
        if unique
            .iter()
            .any(|existing: &RawAlignment| same_placement(existing, &alignment, modulo_length))
        {
            continue;
        }
        unique.push(alignment);
        if unique.len() == 2 {
            break;
        }
    }
    unique
}

fn same_placement(left: &RawAlignment, right: &RawAlignment, modulo_length: Option<usize>) -> bool {
    let left_start = modulo_length
        .map(|length| left.start_reference % length)
        .unwrap_or(left.start_reference);
    let right_start = modulo_length
        .map(|length| right.start_reference % length)
        .unwrap_or(right.start_reference);
    left_start == right_start
        && left.columns.len() == right.columns.len()
        && left
            .columns
            .iter()
            .zip(&right.columns)
            .all(|(left, right)| {
                left.query_base == right.query_base && left.reference_base == right.reference_base
            })
}

fn merge_windows(mut windows: Vec<Window>) -> Vec<Window> {
    windows.sort_unstable_by_key(|window| (window.start, window.end));
    let mut merged: Vec<Window> = Vec::new();
    for window in windows {
        if let Some(last) = merged.last_mut()
            && window.start <= last.end
        {
            last.end = last.end.max(window.end);
            continue;
        }
        merged.push(window);
    }
    merged
}

fn windows_are_economical(windows: &[Window], canonical_length: usize) -> bool {
    if canonical_length < 1_024 {
        return true;
    }
    let total = windows
        .iter()
        .map(|window| window.end - window.start)
        .sum::<usize>();
    total.saturating_mul(4) <= canonical_length.saturating_mul(3)
}

fn partition_ranges(length: usize, count: usize) -> Vec<(usize, usize)> {
    (0..count)
        .map(|index| (index * length / count, (index + 1) * length / count))
        .collect()
}

fn canonical_reference_length(reference: &str, modulo_length: Option<usize>) -> usize {
    modulo_length.unwrap_or(reference.len())
}

fn gapless_alignment(
    query: &str,
    optimal_sequence: &[u8],
    start_reference: usize,
    score: i64,
) -> RawAlignment {
    let columns = query
        .bytes()
        .zip(optimal_sequence.iter().copied())
        .enumerate()
        .map(|(index, (query_base, reference_base))| RawColumn {
            query_base: char::from(query_base),
            reference_base: char::from(reference_base),
            query_index: Some(index),
            reference_index: Some(start_reference + index),
        })
        .collect::<Vec<_>>();
    RawAlignment {
        score,
        start_reference,
        end_reference: start_reference + query.len(),
        metrics: metrics(&columns),
        columns,
    }
}

fn problem_is_provable(
    query: &str,
    profiles: &[Option<EvidenceProfile>],
    reference: &str,
    config: &AlignmentConfig,
    modulo_length: Option<usize>,
) -> bool {
    if query.is_empty()
        || reference.is_empty()
        || profiles.len() != query.len()
        || config.gap_open_score >= 0
        || config.gap_extension_score >= 0
    {
        return false;
    }
    if let Some(length) = modulo_length
        && (length == 0 || query.len() > length)
    {
        return false;
    }

    let Some(rows) = query.len().checked_add(1) else {
        return false;
    };
    let Some(width) = reference.len().checked_add(1) else {
        return false;
    };
    rows.checked_mul(width)
        .is_some_and(|cells| cells <= MAX_ALIGNMENT_CELLS)
}

fn unique_canonical_maximum(scores: [i64; 5]) -> Option<(usize, i64)> {
    let best_score = *scores.iter().max()?;
    let mut best = scores
        .iter()
        .enumerate()
        .filter(|(_, score)| **score == best_score);
    let (index, score) = best.next()?;
    if best.next().is_some() || index >= CANONICAL_BASES.len() {
        return None;
    }
    Some((index, *score))
}

fn occurrence_start_limit(
    reference: &str,
    query_length: usize,
    modulo_length: Option<usize>,
) -> Option<usize> {
    match modulo_length {
        Some(length) => {
            let required = length.checked_add(query_length)?.checked_sub(1)?;
            (reference.len() >= required).then_some(length)
        }
        None => reference
            .len()
            .checked_sub(query_length)
            .and_then(|last_start| last_start.checked_add(1)),
    }
}

#[cfg(test)]
mod tests {
    use crate::alignment::scoring::SCORE_SCALE;
    use crate::model::locus_evidence::EvidenceProfile;

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

    fn profiles(sequence: &str) -> Vec<Option<EvidenceProfile>> {
        sequence
            .bytes()
            .map(|base| {
                let weights = match base {
                    b'A' => [1.0, 0.0, 0.0, 0.0],
                    b'C' => [0.0, 1.0, 0.0, 0.0],
                    b'G' => [0.0, 0.0, 1.0, 0.0],
                    b'T' => [0.0, 0.0, 0.0, 1.0],
                    _ => [0.0; 4],
                };
                Some(EvidenceProfile { weights })
            })
            .collect()
    }

    fn unique(result: UpperBoundPlacement) -> RawAlignment {
        match result {
            UpperBoundPlacement::Unique(alignment) => alignment,
            other => panic!("expected unique proven placement, got {other:?}"),
        }
    }

    #[test]
    fn proves_unique_gapless_upper_bound_placement() {
        let query = "ACGT";
        let alignment = unique(classify(
            query,
            &profiles(query),
            "TTACGTGG",
            &config(),
            None,
        ));

        assert_eq!(alignment.score, 12 * SCORE_SCALE);
        assert_eq!(alignment.start_reference, 2);
        assert_eq!(alignment.end_reference, 6);
        assert_eq!(alignment.metrics.exact_matches, 4);
        assert_eq!(alignment.metrics.gap_opens, 0);
    }

    #[test]
    fn proof_follows_profile_optimum_not_primary_sequence() {
        let alignment = unique(classify(
            "TCGT",
            &profiles("ACGT"),
            "TTACGTGG",
            &config(),
            None,
        ));

        assert_eq!(alignment.score, 12 * SCORE_SCALE);
        assert_eq!(alignment.start_reference, 2);
        assert_eq!(alignment.metrics.exact_matches, 3);
        assert_eq!(alignment.metrics.mismatches, 1);
    }

    #[test]
    fn reports_multiple_upper_bound_placements_as_ambiguous() {
        assert!(matches!(
            classify("AAA", &profiles("AAA"), "AAAAA", &config(), None),
            UpperBoundPlacement::Ambiguous
        ));
    }

    #[test]
    fn distinguishes_unattained_bound_from_unprovable_bound() {
        assert!(matches!(
            classify("ACGT", &profiles("ACGT"), "TTTT", &config(), None),
            UpperBoundPlacement::Unattained
        ));

        let tied = vec![
            Some(EvidenceProfile {
                weights: [0.5, 0.5, 0.0, 0.0],
            }),
            profiles("C")[0],
        ];
        assert!(matches!(
            classify("AC", &tied, "AC", &config(), None),
            UpperBoundPlacement::Unproven
        ));

        let mut weak = config();
        weak.ambiguous_score = -5;
        weak.gap_extension_score = -1;
        let weak_profiles = vec![Some(EvidenceProfile {
            weights: [0.4, 0.3, 0.2, 0.1],
        })];
        assert!(matches!(
            classify("A", &weak_profiles, "A", &weak, None),
            UpperBoundPlacement::Unproven
        ));
    }

    #[test]
    fn over_limit_problem_falls_back_to_authoritative_gotoh() {
        let query = "A".repeat(2_000);
        let reference = "A".repeat(50_000);
        assert!(matches!(
            classify(&query, &profiles(&query), &reference, &config(), None),
            UpperBoundPlacement::Unproven
        ));
    }

    type ColumnSignature = (char, char, Option<usize>, Option<usize>);
    type AlignmentSignature = (i64, usize, usize, Vec<ColumnSignature>);

    fn alignment_signature(alignment: &RawAlignment) -> AlignmentSignature {
        (
            alignment.score,
            alignment.start_reference,
            alignment.end_reference,
            alignment
                .columns
                .iter()
                .map(|column| {
                    (
                        column.query_base,
                        column.reference_base,
                        column.query_index,
                        column.reference_index,
                    )
                })
                .collect(),
        )
    }

    fn assert_same_alignments(actual: &[RawAlignment], expected: &[RawAlignment]) {
        let actual = actual.iter().map(alignment_signature).collect::<Vec<_>>();
        let expected = expected.iter().map(alignment_signature).collect::<Vec<_>>();
        assert_eq!(actual, expected);
    }

    fn require_pruned(
        value: Option<Vec<RawAlignment>>,
        scenario: &str,
    ) -> crate::error::Result<Vec<RawAlignment>> {
        value.ok_or_else(|| {
            crate::error::Error::Alignment(format!("expected exact seeded pruning for {scenario}"))
        })
    }

    #[test]
    fn seeded_pruning_matches_gotoh_for_one_snv() -> crate::error::Result<()> {
        let reference_query = "ACGTCAGTACGATCGTACCTGAGTACGA";
        let mut query = reference_query.to_owned();
        query.replace_range(10..11, "T");
        let reference = format!("TTTT{reference_query}CCCC");
        let query_profiles = profiles(&query);

        let expected =
            crate::alignment::gotoh::align(&query, &query_profiles, &reference, &config(), None)?;
        let actual = require_pruned(
            align_pruned(&query, &query_profiles, &reference, &config(), None)?,
            "one-SNV alignment",
        )?;

        assert_same_alignments(&actual, &expected);
        Ok(())
    }

    #[test]
    fn seeded_pruning_matches_gotoh_for_one_insertion() -> crate::error::Result<()> {
        let reference_query = "ACGTCAGTACGATCGTACCTGAGTACGA";
        let query = format!("{}T{}", &reference_query[..12], &reference_query[12..]);
        let reference = format!("TTTT{reference_query}CCCC");
        let query_profiles = profiles(&query);

        let expected =
            crate::alignment::gotoh::align(&query, &query_profiles, &reference, &config(), None)?;
        let actual = require_pruned(
            align_pruned(&query, &query_profiles, &reference, &config(), None)?,
            "one-insertion alignment",
        )?;

        assert_same_alignments(&actual, &expected);
        Ok(())
    }

    #[test]
    fn seeded_pruning_preserves_rightmost_homopolymer_deletion() -> crate::error::Result<()> {
        let reference_query = "GCCAAAAGTTACGTCAGTACGATCGTAC";
        let query = reference_query.replacen("AAAA", "AAA", 1);
        let reference = format!("TTTT{reference_query}CCCC");
        let query_profiles = profiles(&query);

        let expected =
            crate::alignment::gotoh::align(&query, &query_profiles, &reference, &config(), None)?;
        let actual = require_pruned(
            align_pruned(&query, &query_profiles, &reference, &config(), None)?,
            "homopolymer deletion alignment",
        )?;

        assert_same_alignments(&actual, &expected);
        let deleted = actual[0]
            .columns
            .iter()
            .find(|column| column.query_base == '-')
            .and_then(|column| column.reference_index);
        assert_eq!(deleted, Some(10));
        Ok(())
    }

    #[test]
    fn seeded_pruning_preserves_distinct_equal_placements() -> crate::error::Result<()> {
        let motif = "ACGTCAGTACGATCGTACCTGAGTACGA";
        let mut query = motif.to_owned();
        query.replace_range(10..11, "T");
        let reference = format!("GG{motif}TT{motif}CC");
        let query_profiles = profiles(&query);

        let expected =
            crate::alignment::gotoh::align(&query, &query_profiles, &reference, &config(), None)?;
        let actual = require_pruned(
            align_pruned(&query, &query_profiles, &reference, &config(), None)?,
            "repeated SNV alignment",
        )?;

        assert_eq!(expected.len(), 2);
        assert_same_alignments(&actual, &expected);
        Ok(())
    }

    #[test]
    fn seeded_pruning_supports_circular_non_origin_windows() -> crate::error::Result<()> {
        let reference = "ACGTCAGTACGATCGTACCTGAGTACGATTTTGGGGCCCCAAAATTTT";
        let reference_query = &reference[8..36];
        let mut query = reference_query.to_owned();
        query.replace_range(10..11, "A");
        let working_reference = format!("{reference}{reference}");
        let query_profiles = profiles(&query);

        let expected = crate::alignment::gotoh::align(
            &query,
            &query_profiles,
            &working_reference,
            &config(),
            Some(reference.len()),
        )?;
        let actual = require_pruned(
            align_pruned(
                &query,
                &query_profiles,
                &working_reference,
                &config(),
                Some(reference.len()),
            )?,
            "non-origin circular alignment",
        )?;

        assert_same_alignments(&actual, &expected);
        Ok(())
    }

    fn reverse_complement(sequence: &str) -> String {
        sequence
            .bytes()
            .rev()
            .map(|base| match base {
                b'A' => 'T',
                b'C' => 'G',
                b'G' => 'C',
                b'T' => 'A',
                _ => 'N',
            })
            .collect()
    }

    #[test]
    fn threshold_proof_can_exclude_opposite_orientation() -> crate::error::Result<()> {
        let reference_query = "ACGTCAGTACGATCGTACCTGAGTACGA";
        let mut query = reference_query.to_owned();
        query.replace_range(10..11, "T");
        let reference = format!("TTTT{reference_query}CCCC");
        let query_profiles = profiles(&query);
        let forward = require_pruned(
            align_pruned(&query, &query_profiles, &reference, &config(), None)?,
            "forward threshold source",
        )?;
        let threshold = forward[0].score;

        let reverse_query = reverse_complement(&query);
        let reverse_profiles = profiles(&reverse_query);
        let reverse = align_at_or_above(
            &reverse_query,
            &reverse_profiles,
            &reference,
            &config(),
            None,
            threshold,
        )?;

        assert!(reverse.is_some_and(|placements| placements.is_empty()));
        Ok(())
    }

    #[test]
    fn seeded_pruning_falls_back_for_circular_origin_crossing() -> crate::error::Result<()> {
        let reference = "ACGTCAGTACGATCGTACCTGAGTACGA";
        let mut query = format!("{}{}", &reference[18..], &reference[..18]);
        query.replace_range(5..6, "A");
        let working_reference = format!("{reference}{reference}");
        let query_profiles = profiles(&query);

        assert!(
            align_pruned(
                &query,
                &query_profiles,
                &working_reference,
                &config(),
                Some(reference.len()),
            )?
            .is_none()
        );
        Ok(())
    }

    #[test]
    fn proves_unique_circular_origin_crossing_placement() {
        let reference = "ACGT";
        let working_reference = format!("{reference}{reference}");
        let alignment = unique(classify(
            "GTAC",
            &profiles("GTAC"),
            &working_reference,
            &config(),
            Some(reference.len()),
        ));

        assert_eq!(alignment.start_reference, 2);
        assert_eq!(alignment.end_reference, 6);
        assert_eq!(alignment.metrics.exact_matches, 4);
    }
}
