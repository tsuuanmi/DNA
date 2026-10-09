//! Provably exact alignment acceleration.
//!
//! Two proof tiers share one per-locus profile bound: `upper_bound` places reads
//! that attain the theoretical profile-score maximum directly, and `seeded`
//! prunes Gotoh to q-gram-certified reference windows for low-edit reads. Every
//! entry point returns "unproven" rather than guessing, so callers fall back to
//! full-reference Gotoh.

mod seeded;
mod upper_bound;

#[cfg(test)]
mod tests;

use crate::alignment::scoring::{scaled, substitution_scores};
use crate::alignment::{AlignmentConfig, MAX_ALIGNMENT_CELLS};
use crate::model::nucleotide::Nucleotide;
use crate::read_evidence::EvidenceProfile;

pub(crate) use seeded::{align_at_or_above, align_pruned};
pub(crate) use upper_bound::{UpperBoundPlacement, classify};

#[derive(Debug)]
struct ProfileBound {
    optimal_sequence: Vec<u8>,
    upper_bound: i64,
    edit_loss_floor: i64,
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
        optimal_sequence.push(Nucleotide::ALL[best_index].as_byte());
    }

    (edit_loss_floor > 0).then_some(ProfileBound {
        optimal_sequence,
        upper_bound,
        edit_loss_floor,
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
    if best.next().is_some() || index >= Nucleotide::ALL.len() {
        return None;
    }
    Some((index, *score))
}

fn canonical_reference_length(reference: &str, modulo_length: Option<usize>) -> usize {
    modulo_length.unwrap_or(reference.len())
}
