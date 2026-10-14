//! Score-bounded q-gram pruning of Gotoh to certified reference windows.

use std::cmp::Reverse;

use memchr::memmem::Finder;

use crate::alignment::AlignmentConfig;
use crate::alignment::gotoh;
use crate::alignment::traceback::RawAlignment;
use dna_kernel::error::Result;
use dna_kernel::model::nucleotide::Nucleotide;
use dna_kernel::read_evidence::EvidenceProfile;

use super::{
    ProfileBound, WILDCARD, canonical_reference_length, problem_is_provable, profile_bound,
};

const PROBE_COUNT: usize = 4;
const PROBE_MAX_SEED_LENGTH: usize = 16;
const MAX_SEED_HITS: usize = 512;
const MAX_CERTIFIED_SEEDS: usize = 128;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Window {
    start: usize,
    end: usize,
}

/// Outcome of seeded pruning for one orientation.
#[derive(Debug)]
pub(crate) enum Pruned {
    /// The orientation's exact best placements.
    Proven(Vec<RawAlignment>),
    /// No proof; `lower_bound` is the score of a valid alignment when the
    /// seed probes found one.
    Unproven { lower_bound: Option<i64> },
}

/// Attempts an exact seeded-pruning alignment for one orientation.
///
/// A cheap seed-derived local alignment supplies only a valid lower-bound score.
/// That score is converted to a maximum edit budget from the profile upper bound.
/// The q-gram lemma then proves that every alignment able to tie or beat the lower
/// bound contains at least one exact seed, allowing Gotoh to run only on the union
/// of certified windows. `Unproven` means the proof or cost bound was not
/// strong enough and the caller must use full-reference Gotoh.
pub(crate) fn align_pruned(
    query: &str,
    profiles: &[Option<EvidenceProfile>],
    reference: &str,
    config: &AlignmentConfig,
    modulo_length: Option<usize>,
) -> Result<Pruned> {
    let unproven = Pruned::Unproven { lower_bound: None };
    if !problem_is_provable(query, profiles, reference, config, modulo_length) {
        return Ok(unproven);
    }
    let Some(bound) = profile_bound(profiles, config) else {
        return Ok(unproven);
    };
    let Some(lower_bound) =
        seed_lower_bound(query, profiles, reference, config, modulo_length, &bound)?
    else {
        return Ok(unproven);
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
    Ok(match result {
        Some(placements) if !placements.is_empty() => Pruned::Proven(placements),
        _ => Pruned::Unproven {
            lower_bound: Some(lower_bound),
        },
    })
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
    let Ok(edit_budget) = usize::try_from(deficit / bound.edit_loss_floor) else {
        return Ok(None);
    };
    let Ok(drift) = usize::try_from(deficit / bound.deletion_loss_floor) else {
        return Ok(None);
    };
    // Each non-canonical reference symbol can break one seed at no edit cost.
    let Some(seed_count) = edit_budget
        .checked_add(1)
        .and_then(|count| count.checked_add(ambiguous_reference_bases(reference, modulo_length)))
    else {
        return Ok(None);
    };
    if seed_count > query.len() {
        return Ok(None);
    }
    if seed_count > MAX_CERTIFIED_SEEDS {
        return Ok(None);
    }

    let Some(windows) = certified_windows(
        &bound.optimal_sequence,
        reference,
        modulo_length,
        drift,
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

    let placements = align_windows(query, profiles, reference, config, modulo_length, &windows)?;
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
        let Some((run_start, run_end)) =
            longest_seed(&bound.optimal_sequence, block_start, block_end)
        else {
            continue;
        };
        let run_length = run_end - run_start;
        let seed_length = run_length.min(PROBE_MAX_SEED_LENGTH);
        let seed_start = run_start + (run_length - seed_length) / 2;
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
        let placements = gotoh::align_window(
            query,
            profiles,
            reference,
            window.start..window.end,
            config,
            modulo_length,
        )?;
        if let Some(first) = placements.first() {
            best = Some(best.map_or(first.score, |score: i64| score.max(first.score)));
        }
    }
    Ok(best)
}

/// Certified windows for `seed_count` disjoint wildcard-free seeds.
///
/// An alignment within the score threshold leaves at least one seed exact, so
/// that seed occurs at its aligned reference position; `drift` bounds the
/// deleted reference bases on either side.
fn certified_windows(
    optimal_sequence: &[u8],
    reference: &str,
    modulo_length: Option<usize>,
    drift: usize,
    seed_count: usize,
) -> Option<Vec<Window>> {
    let mut windows = Vec::new();
    let mut hit_count = 0;
    for (seed_start, seed_end) in disjoint_seeds(optimal_sequence, seed_count)? {
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
                drift,
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

/// The reference window holding every alignment whose query locus
/// `seed_start` lies on reference position `hit`, given at most `drift`
/// deleted reference bases.
///
/// A circular reference is doubled, so a window that would start before the
/// origin is taken one circle later instead; either way it starts inside the
/// first circle or covers the origin from the left, as the whole-reference
/// search does.
fn window_from_seed(
    hit: usize,
    seed_start: usize,
    seed_end: usize,
    query_length: usize,
    drift: usize,
    reference: &str,
    modulo_length: Option<usize>,
) -> Option<Window> {
    let prefix = seed_start.checked_add(drift)?;
    let suffix = query_length.checked_sub(seed_end)?.checked_add(drift)?;
    let seed_length = seed_end.checked_sub(seed_start)?;
    let window = match modulo_length {
        None => Window {
            start: hit.saturating_sub(prefix),
            end: hit
                .checked_add(seed_length)?
                .checked_add(suffix)?
                .min(reference.len()),
        },
        Some(length) => {
            let hit = if hit < prefix {
                hit.checked_add(length)?
            } else {
                hit
            };
            let end = hit.checked_add(seed_length)?.checked_add(suffix)?;
            if end > reference.len() {
                return None;
            }
            Window {
                start: hit - prefix,
                end,
            }
        }
    };
    (window.start < window.end).then_some(window)
}

fn align_windows(
    query: &str,
    profiles: &[Option<EvidenceProfile>],
    reference: &str,
    config: &AlignmentConfig,
    modulo_length: Option<usize>,
    windows: &[Window],
) -> Result<Vec<RawAlignment>> {
    let mut placements = Vec::new();
    for window in windows {
        let local = gotoh::align_window(
            query,
            profiles,
            reference,
            window.start..window.end,
            config,
            modulo_length,
        )?;
        placements.extend(local.into_iter().map(|mut alignment| {
            // A placement starting in the second circle repeats one starting
            // in the first, which is how the whole-reference search reports it.
            if let Some(length) = modulo_length
                && alignment.start_reference >= length
            {
                alignment.unshift_reference(length);
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

/// Non-canonical symbols one alignment can cover: an alignment consumes at
/// most one circle of a circular reference.
fn ambiguous_reference_bases(reference: &str, modulo_length: Option<usize>) -> usize {
    reference.as_bytes()[..canonical_reference_length(reference, modulo_length)]
        .iter()
        .filter(|base| Nucleotide::from_char(char::from(**base)).is_none())
        .count()
}

/// The longest wildcard-free run of `optimal_sequence[start..end]`, the
/// leftmost on ties.
fn longest_seed(optimal_sequence: &[u8], start: usize, end: usize) -> Option<(usize, usize)> {
    let mut best: Option<(usize, usize)> = None;
    let mut run_start = start;
    for index in start..=end {
        if index == end || optimal_sequence.get(index) == Some(&WILDCARD) {
            if index > run_start
                && best
                    .is_none_or(|(best_start, best_end)| index - run_start > best_end - best_start)
            {
                best = Some((run_start, index));
            }
            run_start = index + 1;
        }
    }
    best
}

/// `count` disjoint wildcard-free seeds whose shortest is as long as possible,
/// in query order; `None` when the query has too few seed loci.
fn disjoint_seeds(optimal_sequence: &[u8], count: usize) -> Option<Vec<(usize, usize)>> {
    let runs = optimal_sequence
        .split(|base| *base == WILDCARD)
        .scan(0, |start, run| {
            let range = (*start, *start + run.len());
            *start += run.len() + 1;
            Some(range)
        })
        .filter(|(start, end)| end > start)
        .collect::<Vec<_>>();
    let capacity = |length: usize| {
        runs.iter()
            .map(|(start, end)| (end - start) / length)
            .sum::<usize>()
    };
    if count == 0 || capacity(1) < count {
        return None;
    }
    // The largest seed length that still yields `count` seeds.
    let (mut shortest, mut longest) = (1, runs.iter().map(|(start, end)| end - start).max()?);
    while shortest < longest {
        let middle = shortest + (longest - shortest).div_ceil(2);
        if capacity(middle) >= count {
            shortest = middle;
        } else {
            longest = middle - 1;
        }
    }
    let mut seeds = runs
        .iter()
        .flat_map(|(start, end)| {
            partition_ranges(end - start, (end - start) / shortest)
                .into_iter()
                .map(move |(from, to)| (start + from, start + to))
        })
        .collect::<Vec<_>>();
    seeds.sort_by_key(|(start, end)| (Reverse(end - start), *start));
    seeds.truncate(count);
    seeds.sort_unstable();
    Some(seeds)
}

fn partition_ranges(length: usize, count: usize) -> Vec<(usize, usize)> {
    (0..count)
        .map(|index| (index * length / count, (index + 1) * length / count))
        .collect()
}
