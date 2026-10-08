//! Direct placement of reads that attain the profile-score upper bound.

use memchr::memmem::Finder;

use crate::alignment::traceback::{RawAlignment, RawColumn, metrics};
use crate::config::AlignmentConfig;
use crate::model::locus_evidence::EvidenceProfile;

use super::{problem_is_provable, profile_bound};

/// Outcome of the upper-bound proof for one orientation.
#[derive(Debug)]
pub(crate) enum UpperBoundPlacement {
    /// The proof preconditions do not hold; use Gotoh.
    Unproven,
    /// No placement attains the upper bound.
    Unattained,
    /// Exactly one placement attains the upper bound.
    Unique(RawAlignment),
    /// More than one placement attains the upper bound.
    Ambiguous,
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
