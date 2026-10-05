//! Score-bounded seeded alignment proof.

use std::collections::BTreeSet;

use memchr::memmem::Finder;

use crate::alignment::scoring::{scaled, substitution_index, substitution_scores};
use crate::alignment::traceback::{RawAlignment, RawColumn, metrics};
use crate::config::{AlignmentConfig, MAX_ALIGNMENT_CELLS};
use crate::model::locus_evidence::EvidenceProfile;

const CANONICAL_BASES: &[u8; 4] = b"ACGT";
const MAX_SUBSTITUTION_EDITS: usize = 2;
const SEED_COUNT: usize = MAX_SUBSTITUTION_EDITS + 1;
const MAX_CANDIDATE_STARTS: usize = 1_024;

pub(crate) struct SeededSearch {
    upper_bound: i64,
    minimum_substitution_deficit: i64,
    minimum_gap_deficit: i64,
    best_score: Option<i64>,
    placements: Vec<RawAlignment>,
}

impl SeededSearch {
    pub(crate) const fn best_score(&self) -> Option<i64> {
        self.best_score
    }

    pub(crate) fn complete_at(&self, threshold: i64) -> bool {
        if threshold > self.upper_bound {
            return true;
        }
        let Some(deficit) = self.upper_bound.checked_sub(threshold) else {
            return false;
        };
        let Some(substitution_limit) = self
            .minimum_substitution_deficit
            .checked_mul(SEED_COUNT as i64)
        else {
            return false;
        };
        deficit < self.minimum_gap_deficit && deficit < substitution_limit
    }

    pub(crate) fn into_placements(self) -> Vec<RawAlignment> {
        self.placements
    }
}

pub(crate) fn search(
    query: &str,
    profiles: &[Option<EvidenceProfile>],
    reference: &str,
    config: &AlignmentConfig,
    modulo_length: Option<usize>,
) -> Option<SeededSearch> {
    if !problem_is_searchable(query, profiles, reference, config, modulo_length) {
        return None;
    }

    let gap_extension = scaled(config.gap_extension_score);
    let open_and_extend = scaled(config.gap_open_score).checked_add(gap_extension)?;
    let deletion_deficit = open_and_extend.checked_neg()?;
    let mut minimum_gap_deficit = deletion_deficit;
    let mut minimum_substitution_deficit = i64::MAX;
    let mut upper_bound = 0_i64;
    let mut optimal_sequence = Vec::with_capacity(query.len());
    let mut row_scores = Vec::with_capacity(query.len());

    for profile in profiles {
        let scores = substitution_scores(*profile, config);
        let (best_index, best_score, second_best) = unique_canonical_maximum(scores)?;
        if best_score <= gap_extension {
            return None;
        }

        upper_bound = upper_bound.checked_add(best_score)?;
        minimum_substitution_deficit =
            minimum_substitution_deficit.min(best_score.checked_sub(second_best)?);
        minimum_gap_deficit =
            minimum_gap_deficit.min(best_score.checked_sub(open_and_extend)?);
        optimal_sequence.push(CANONICAL_BASES[best_index]);
        row_scores.push(scores);
    }

    if minimum_substitution_deficit <= 0 || minimum_gap_deficit <= 0 {
        return None;
    }

    let starts = candidate_starts(
        &optimal_sequence,
        reference,
        modulo_length,
        MAX_CANDIDATE_STARTS,
    )?;
    let mut best_score = None;
    let mut best_starts = Vec::new();

    for start in starts {
        let score = score_gapless(start, &row_scores, reference)?;
        match best_score {
            None => {
                best_score = Some(score);
                best_starts.push(start);
            }
            Some(best) if score > best => {
                best_score = Some(score);
                best_starts.clear();
                best_starts.push(start);
            }
            Some(best) if score == best => best_starts.push(start),
            Some(_) => {}
        }
    }

    let placements = match best_score {
        Some(score) => best_starts
            .into_iter()
            .map(|start| gapless_alignment(query, reference, start, score))
            .collect::<Option<Vec<_>>>()?,
        None => Vec::new(),
    };

    Some(SeededSearch {
        upper_bound,
        minimum_substitution_deficit,
        minimum_gap_deficit,
        best_score,
        placements,
    })
}

fn problem_is_searchable(
    query: &str,
    profiles: &[Option<EvidenceProfile>],
    reference: &str,
    config: &AlignmentConfig,
    modulo_length: Option<usize>,
) -> bool {
    if query.len() < SEED_COUNT
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

fn unique_canonical_maximum(scores: [i64; 5]) -> Option<(usize, i64, i64)> {
    let (best_index, best_score) = scores
        .iter()
        .copied()
        .enumerate()
        .max_by_key(|(_, score)| *score)?;
    if best_index >= CANONICAL_BASES.len() {
        return None;
    }

    let mut second_best = i64::MIN;
    let mut best_count = 0;
    for (index, score) in scores.into_iter().enumerate() {
        if score == best_score {
            best_count += 1;
        }
        if index != best_index {
            second_best = second_best.max(score);
        }
    }
    (best_count == 1).then_some((best_index, best_score, second_best))
}

fn candidate_starts(
    optimal_sequence: &[u8],
    reference: &str,
    modulo_length: Option<usize>,
    maximum_candidates: usize,
) -> Option<BTreeSet<usize>> {
    let start_limit = valid_start_limit(reference, optimal_sequence.len(), modulo_length)?;
    let haystack_end = match modulo_length {
        Some(_) => start_limit
            .checked_add(optimal_sequence.len())?
            .checked_sub(1)?,
        None => reference.len(),
    };
    let haystack = reference.as_bytes().get(..haystack_end)?;
    let mut starts = BTreeSet::new();

    for (seed_start, seed_end) in seed_ranges(optimal_sequence.len())? {
        let seed = optimal_sequence.get(seed_start..seed_end)?;
        let finder = Finder::new(seed);
        let mut search_start = 0;
        while search_start < haystack.len() {
            let Some(relative) = finder.find(&haystack[search_start..]) else {
                break;
            };
            let hit = search_start.checked_add(relative)?;
            if hit >= seed_start {
                let start = hit - seed_start;
                if start < start_limit
                    && start.checked_add(optimal_sequence.len())? <= reference.len()
                {
                    starts.insert(start);
                    if starts.len() > maximum_candidates {
                        return None;
                    }
                }
            }
            search_start = hit.checked_add(1)?;
        }
    }

    Some(starts)
}

fn seed_ranges(length: usize) -> Option<[(usize, usize); SEED_COUNT]> {
    if length < SEED_COUNT {
        return None;
    }
    let base = length / SEED_COUNT;
    let remainder = length % SEED_COUNT;
    let mut ranges = [(0, 0); SEED_COUNT];
    let mut start = 0;
    for (index, range) in ranges.iter_mut().enumerate() {
        let size = base + usize::from(index < remainder);
        let end = start.checked_add(size)?;
        *range = (start, end);
        start = end;
    }
    (start == length).then_some(ranges)
}

fn valid_start_limit(
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

fn score_gapless(start: usize, row_scores: &[[i64; 5]], reference: &str) -> Option<i64> {
    let bases = reference
        .as_bytes()
        .get(start..start.checked_add(row_scores.len())?)?;
    row_scores
        .iter()
        .zip(bases)
        .try_fold(0_i64, |score, (scores, reference_base)| {
            score.checked_add(scores[substitution_index(*reference_base)])
        })
}

fn gapless_alignment(
    query: &str,
    reference: &str,
    start_reference: usize,
    score: i64,
) -> Option<RawAlignment> {
    let reference_bases = reference
        .as_bytes()
        .get(start_reference..start_reference.checked_add(query.len())?)?;
    let columns = query
        .bytes()
        .zip(reference_bases.iter().copied())
        .enumerate()
        .map(|(index, (query_base, reference_base))| RawColumn {
            query_base: char::from(query_base),
            reference_base: char::from(reference_base),
            query_index: Some(index),
            reference_index: Some(start_reference + index),
        })
        .collect::<Vec<_>>();

    Some(RawAlignment {
        score,
        start_reference,
        end_reference: start_reference + query.len(),
        metrics: metrics(&columns),
        columns,
    })
}

#[cfg(test)]
#[derive(Debug)]
pub(crate) struct SeededProof {
    pub(crate) score: i64,
    pub(crate) placements: Vec<RawAlignment>,
}

#[cfg(test)]
pub(crate) fn classify(
    query: &str,
    profiles: &[Option<EvidenceProfile>],
    reference: &str,
    config: &AlignmentConfig,
    modulo_length: Option<usize>,
) -> Option<SeededProof> {
    let search = search(query, profiles, reference, config, modulo_length)?;
    let score = search.best_score()?;
    if !search.complete_at(score) {
        return None;
    }
    Some(SeededProof {
        score,
        placements: search.into_placements(),
    })
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

    #[test]
    fn proves_unique_one_substitution_placement() {
        let Some(proof) = classify("ACGT", &profiles("ACGT"), "TTACATGG", &config(), None) else {
            panic!("one-substitution placement should be provable");
        };

        assert_eq!(proof.score, 4 * SCORE_SCALE);
        assert_eq!(proof.placements.len(), 1);
        assert_eq!(proof.placements[0].start_reference, 2);
        assert_eq!(proof.placements[0].end_reference, 6);
        assert_eq!(proof.placements[0].metrics.exact_matches, 3);
        assert_eq!(proof.placements[0].metrics.mismatches, 1);
        assert_eq!(proof.placements[0].metrics.gap_opens, 0);
    }

    #[test]
    fn proves_all_equally_best_repeated_substitution_placements() {
        let Some(proof) = classify("ACGT", &profiles("ACGT"), "ACATGGACAT", &config(), None) else {
            panic!("repeated best placements should be completely enumerated");
        };

        assert_eq!(proof.score, 4 * SCORE_SCALE);
        assert_eq!(
            proof
                .placements
                .iter()
                .map(|alignment| alignment.start_reference)
                .collect::<Vec<_>>(),
            vec![0, 6]
        );
    }

    #[test]
    fn refuses_candidate_when_a_gap_can_rival_its_score() {
        let mut weak_gap = config();
        weak_gap.gap_open_score = -1;
        weak_gap.gap_extension_score = -1;

        assert!(classify("ACGT", &profiles("ACGT"), "TTACATGG", &weak_gap, None,).is_none());
    }

    #[test]
    fn refuses_candidate_outside_two_substitution_seed_bound() {
        assert!(classify("ACGTAC", &profiles("ACGTAC"), "TTATGTTGG", &config(), None,).is_none());
    }

    #[test]
    fn proves_one_substitution_across_circular_origin() {
        let reference = "ACGT";
        let working_reference = format!("{reference}{reference}");
        let Some(proof) = classify(
            "GTTC",
            &profiles("GTTC"),
            &working_reference,
            &config(),
            Some(reference.len()),
        ) else {
            panic!("origin-crossing substitution should be provable");
        };

        assert_eq!(proof.placements.len(), 1);
        assert_eq!(proof.placements[0].start_reference, 2);
        assert_eq!(proof.placements[0].end_reference, 6);
        assert_eq!(proof.placements[0].metrics.mismatches, 1);
    }
}
