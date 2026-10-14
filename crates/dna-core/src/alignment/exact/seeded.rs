//! Score-bounded q-gram pruning of Gotoh to certified reference windows.

use std::cmp::Reverse;

use memchr::memmem::Finder;

use crate::alignment::gotoh;
use crate::alignment::traceback::RawAlignment;
use dna_kernel::error::{AlignmentError, Result};

use super::{Problem, WILDCARD, WorkingReference};

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

impl Pruned {
    /// The proven placements, if any.
    pub(crate) fn into_proven(self) -> Option<Vec<RawAlignment>> {
        match self {
            Self::Proven(placements) => Some(placements),
            Self::Unproven { .. } => None,
        }
    }
}

/// Attempts an exact seeded-pruning alignment for one orientation.
///
/// A cheap seed-derived local alignment supplies only a valid lower-bound score.
/// That score is converted to a maximum edit budget from the profile upper bound.
/// The q-gram lemma then proves that every alignment able to tie or beat the lower
/// bound contains at least one exact seed, allowing Gotoh to run only on the union
/// of certified windows. `Unproven` means the proof or cost bound was not
/// strong enough and the caller must use full-reference Gotoh.
pub(crate) fn align_pruned(problem: &Problem<'_>) -> Result<Pruned> {
    let Some(lower_bound) = seed_lower_bound(problem)? else {
        return Ok(Pruned::Unproven { lower_bound: None });
    };
    proven_at_lower_bound(align_at_or_above(problem, lower_bound)?, lower_bound)
}

/// Reads the certified search at the seed lower bound, which a valid
/// alignment attains: a complete proof that finds no placement there
/// contradicts that alignment, so it is an internal invariant failure.
pub(super) fn proven_at_lower_bound(
    certified: Option<Vec<RawAlignment>>,
    lower_bound: i64,
) -> Result<Pruned> {
    match certified {
        Some(placements) if placements.is_empty() => Err(AlignmentError::Inconsistent(
            "seeded proof found no placement at its own lower bound",
        )
        .into()),
        Some(placements) => Ok(Pruned::Proven(placements)),
        None => Ok(Pruned::Unproven {
            lower_bound: Some(lower_bound),
        }),
    }
}

/// Proves and evaluates every placement in one orientation that can reach `threshold`.
///
/// `Some(empty)` is a proof that this orientation cannot reach the threshold.
/// `None` means the proof is too weak or expensive, so the caller must use Gotoh.
pub(crate) fn align_at_or_above(
    problem: &Problem<'_>,
    threshold: i64,
) -> Result<Option<Vec<RawAlignment>>> {
    let bound = &problem.bound;
    let reference = problem.reference;
    if threshold > bound.upper_bound {
        return Ok(Some(Vec::new()));
    }
    let Some(deficit) = bound.upper_bound.checked_sub(threshold) else {
        return Ok(None);
    };
    let Some(edit_budget) = deficit
        .checked_div(bound.edit_loss_floor)
        .and_then(|budget| usize::try_from(budget).ok())
    else {
        return Ok(None);
    };
    let Some(drift) = deficit
        .checked_div(bound.deletion_loss_floor)
        .and_then(|drift| usize::try_from(drift).ok())
    else {
        return Ok(None);
    };
    // Each non-canonical reference symbol can break one seed at no edit cost.
    let Some(seed_count) = edit_budget
        .checked_add(1)
        .and_then(|count| count.checked_add(reference.ambiguous_bases))
    else {
        return Ok(None);
    };
    if seed_count > MAX_CERTIFIED_SEEDS {
        return Ok(None);
    }

    let Some(windows) = certified_windows(&bound.optimal_sequence, reference, drift, seed_count)
    else {
        return Ok(None);
    };
    if windows.is_empty() {
        return Ok(Some(Vec::new()));
    }
    if !windows_are_economical(&windows, reference.canonical_length()) {
        return Ok(None);
    }

    let placements = align_windows(problem, &windows)?;
    Ok(Some(best_placements_at_or_above(
        placements,
        threshold,
        reference.modulo_length,
    )))
}

/// The best score Gotoh finds inside windows around probe-seed hits: a valid
/// alignment's score, so a lower bound on the orientation's optimum.
fn seed_lower_bound(problem: &Problem<'_>) -> Result<Option<i64>> {
    let query_length = problem.query.len();
    if query_length < PROBE_COUNT {
        return Ok(None);
    }
    let optimal_sequence = &problem.bound.optimal_sequence;
    let reference = problem.reference;
    let margin = 32.min((query_length / 4).max(4));
    let mut windows = Vec::new();
    let mut hit_count = 0;

    for (block_start, block_end) in partition_ranges(query_length, PROBE_COUNT) {
        let Some((run_start, run_end)) = longest_seed(optimal_sequence, block_start, block_end)
        else {
            continue;
        };
        let run_length = run_end - run_start;
        let seed_length = run_length.min(PROBE_MAX_SEED_LENGTH);
        let seed_start = run_start + (run_length - seed_length) / 2;
        let seed_end = seed_start + seed_length;
        let seed = &optimal_sequence[seed_start..seed_end];

        for hit in seed_hits(reference, seed, false) {
            hit_count += 1;
            if hit_count > MAX_SEED_HITS {
                return Ok(None);
            }
            if let Some(window) =
                window_from_seed(hit, seed_start, seed_end, query_length, margin, reference)
            {
                windows.push(window);
            }
        }
    }

    let windows = merge_windows(windows);
    // An uneconomical probe would cost about as much as the full-reference
    // Gotoh the caller runs anyway.
    if windows.is_empty() || !windows_are_economical(&windows, reference.canonical_length()) {
        return Ok(None);
    }

    let mut best = None;
    for window in windows {
        let placements = gotoh::align_window(
            problem.query,
            problem.profiles,
            reference.sequence,
            window.start..window.end,
            problem.config,
            reference.modulo_length,
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
    reference: &WorkingReference<'_>,
    drift: usize,
    seed_count: usize,
) -> Option<Vec<Window>> {
    let mut windows = Vec::new();
    let mut hit_count = 0;
    for (seed_start, seed_end) in disjoint_seeds(optimal_sequence, seed_count)? {
        let seed = &optimal_sequence[seed_start..seed_end];
        for hit in seed_hits(reference, seed, true) {
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
            )?;
            windows.push(window);
        }
    }
    Some(merge_windows(windows))
}

fn seed_hits(reference: &WorkingReference<'_>, seed: &[u8], include_seam: bool) -> Vec<usize> {
    if seed.is_empty() {
        return Vec::new();
    }
    let sequence = reference.sequence.as_bytes();
    let haystack_end = match reference.modulo_length {
        Some(length) if include_seam => length
            .checked_add(seed.len() - 1)
            .map_or(sequence.len(), |end| end.min(sequence.len())),
        _ => reference.canonical_length(),
    };
    let haystack = &sequence[..haystack_end];
    let finder = Finder::new(seed);
    let mut hits = Vec::new();
    let mut search_start = 0;

    while search_start <= haystack.len().saturating_sub(seed.len()) {
        let Some(relative) = finder.find(&haystack[search_start..]) else {
            break;
        };
        let hit = search_start + relative;
        if reference.modulo_length.is_some_and(|length| hit >= length) {
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
/// search does. A window that would start before the doubled reference itself
/// starts at its beginning: no alignment starts earlier, and alignments that
/// place the seed on the first copy of `hit` end inside the window too, so the
/// window stays complete.
fn window_from_seed(
    hit: usize,
    seed_start: usize,
    seed_end: usize,
    query_length: usize,
    drift: usize,
    reference: &WorkingReference<'_>,
) -> Option<Window> {
    let prefix = seed_start.checked_add(drift)?;
    let suffix = query_length.checked_sub(seed_end)?.checked_add(drift)?;
    let seed_length = seed_end.checked_sub(seed_start)?;
    let sequence_length = reference.sequence.len();
    let hit = match reference.modulo_length {
        Some(length) if hit < prefix => hit.checked_add(length)?,
        _ => hit,
    };
    let end = hit.checked_add(seed_length)?.checked_add(suffix)?;
    let end = match reference.modulo_length {
        None => end.min(sequence_length),
        Some(_) if end > sequence_length => return None,
        Some(_) => end,
    };
    let window = Window {
        start: hit.saturating_sub(prefix),
        end,
    };
    (window.start < window.end).then_some(window)
}

fn align_windows(problem: &Problem<'_>, windows: &[Window]) -> Result<Vec<RawAlignment>> {
    let reference = problem.reference;
    let mut placements = Vec::new();
    for window in windows {
        let local = gotoh::align_window(
            problem.query,
            problem.profiles,
            reference.sequence,
            window.start..window.end,
            problem.config,
            reference.modulo_length,
        )?;
        placements.extend(local.into_iter().map(|mut alignment| {
            // A placement starting in the second circle repeats one starting
            // in the first, which is how the whole-reference search reports it.
            if let Some(length) = reference.modulo_length
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
