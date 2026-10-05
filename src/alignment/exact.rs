//! Provably exact upper-bound alignment fast path.

use memchr::memmem::Finder;

use crate::alignment::scoring::{scaled, substitution_scores};
use crate::alignment::traceback::{RawAlignment, RawColumn, metrics};
use crate::config::{AlignmentConfig, MAX_ALIGNMENT_CELLS};
use crate::model::locus_evidence::EvidenceProfile;

const CANONICAL_BASES: &[u8; 4] = b"ACGT";

#[derive(Debug)]
pub(crate) enum UpperBoundPlacement {
    Unproven,
    Unattained,
    Unique(RawAlignment),
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

    let gap_extension = scaled(config.gap_extension_score);
    let mut optimal_sequence = Vec::with_capacity(query.len());
    let mut upper_bound = 0_i64;

    for profile in profiles {
        let scores = substitution_scores(*profile, config);
        let Some((best_index, best_score)) = unique_canonical_maximum(scores) else {
            return UpperBoundPlacement::Unproven;
        };
        if best_score <= gap_extension {
            return UpperBoundPlacement::Unproven;
        }
        let Some(score) = upper_bound.checked_add(best_score) else {
            return UpperBoundPlacement::Unproven;
        };
        upper_bound = score;
        optimal_sequence.push(CANONICAL_BASES[best_index]);
    }

    let Some(limit) = occurrence_start_limit(reference, optimal_sequence.len(), modulo_length)
    else {
        return UpperBoundPlacement::Unproven;
    };

    let haystack_end = limit + optimal_sequence.len() - 1;
    let haystack = &reference.as_bytes()[..haystack_end];
    let finder = Finder::new(&optimal_sequence);
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
    let columns = query
        .bytes()
        .zip(optimal_sequence)
        .enumerate()
        .map(|(index, (query_base, reference_base))| RawColumn {
            query_base: char::from(query_base),
            reference_base: char::from(reference_base),
            query_index: Some(index),
            reference_index: Some(start_reference + index),
        })
        .collect::<Vec<_>>();

    UpperBoundPlacement::Unique(RawAlignment {
        score: upper_bound,
        start_reference,
        end_reference: start_reference + query.len(),
        metrics: metrics(&columns),
        columns,
    })
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
