//! Score-bounded q-gram pruning of Gotoh to certified reference windows.

use memchr::memmem::Finder;

use crate::alignment::AlignmentConfig;
use crate::alignment::gotoh;
use crate::alignment::traceback::RawAlignment;
use dna_kernel::error::Result;
use dna_kernel::read_evidence::EvidenceProfile;

use super::{ProfileBound, canonical_reference_length, problem_is_provable, profile_bound};

const PROBE_COUNT: usize = 4;
const PROBE_MAX_SEED_LENGTH: usize = 16;
const MAX_SEED_HITS: usize = 128;
const MAX_CERTIFIED_SEEDS: usize = 32;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Window {
    start: usize,
    end: usize,
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
            modulo_length.map_or(alignment.start_reference, |length| {
                alignment.start_reference % length
            }),
            alignment.end_reference,
        )
    });

    let mut unique = Vec::new();
    for alignment in placements {
        if unique
            .iter()
            .any(|existing: &RawAlignment| existing.same_placement(&alignment, modulo_length))
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
