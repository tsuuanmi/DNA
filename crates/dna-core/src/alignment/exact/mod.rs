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
use dna_kernel::model::nucleotide::Nucleotide;
use dna_kernel::read_evidence::EvidenceProfile;

pub(crate) use seeded::{Pruned, align_at_or_above, align_pruned};
pub(crate) use upper_bound::{UpperBoundPlacement, classify};

/// Query locus that seeds never contain: its best substitution is not one
/// unique canonical base, or loses less than one edit unit to the runner-up.
const WILDCARD: u8 = b'N';

/// Score bounds shared by both proof tiers.
#[derive(Debug)]
struct ProfileBound {
    /// Unique profile-optimal base per query locus, or `WILDCARD`.
    optimal_sequence: Vec<u8>,
    /// Sum of every locus's best substitution score over all reference symbols.
    upper_bound: i64,
    /// Minimum score loss of one edit unit that breaks an exact seed.
    edit_loss_floor: i64,
    /// Minimum score loss of one deleted reference base.
    deletion_loss_floor: i64,
}

/// Derives the profile upper bound and the seed edit-loss floor.
///
/// Every alignment scores `upper_bound` minus non-negative per-event losses:
/// an aligned locus loses its best score minus the substitution it receives,
/// an inserted locus loses its best score minus one gap extension, and a
/// deletion run loses its opening plus one extension per base. A query
/// interval aligned without gaps and with every non-wildcard locus on its
/// optimal base is an exact seed. Breaking one costs at least
/// `edit_loss_floor`: an inserted locus loses at least the smallest
/// best-minus-extension, a deletion run at least one open plus one extension,
/// and a non-wildcard substitution at least its best-minus-runner-up loss,
/// which is at least the floor by the wildcard rule.
fn profile_bound(
    profiles: &[Option<EvidenceProfile>],
    config: &AlignmentConfig,
) -> Option<ProfileBound> {
    let gap_extension = scaled(config.gap_extension_score);
    let gap_open = scaled(config.gap_open_score);
    if gap_extension >= 0 || gap_open > 0 {
        return None;
    }
    let deletion_run_loss = gap_open.checked_add(gap_extension)?.checked_neg()?;

    let mut loci = Vec::with_capacity(profiles.len());
    let mut insertion_loss = i64::MAX;
    for profile in profiles {
        let scores = substitution_scores(*profile, config);
        let best_score = *scores.iter().max()?;
        if best_score <= gap_extension {
            return None;
        }
        insertion_loss = insertion_loss.min(best_score.checked_sub(gap_extension)?);
        loci.push((scores, best_score));
    }
    let edit_loss_floor = insertion_loss.min(deletion_run_loss);
    if edit_loss_floor <= 0 {
        return None;
    }

    let mut optimal_sequence = Vec::with_capacity(loci.len());
    let mut upper_bound = 0_i64;
    for (scores, best_score) in loci {
        upper_bound = upper_bound.checked_add(best_score)?;
        optimal_sequence.push(canonical_seed_base(scores, edit_loss_floor));
    }

    Some(ProfileBound {
        optimal_sequence,
        upper_bound,
        edit_loss_floor,
        deletion_loss_floor: gap_extension.checked_neg()?,
    })
}

/// The locus's seed base: the canonical base whose substitution beats every
/// other canonical base by at least `edit_loss_floor`, or `WILDCARD`.
///
/// A non-canonical reference symbol is not a competitor here; the seeded tier
/// budgets one broken seed per such symbol instead.
fn canonical_seed_base(scores: [i64; 5], edit_loss_floor: i64) -> u8 {
    let canonical = &scores[..Nucleotide::ALL.len()];
    let Some((best_index, best_canonical)) = canonical
        .iter()
        .enumerate()
        .max_by_key(|(_, score)| **score)
    else {
        return WILDCARD;
    };
    let separated = canonical
        .iter()
        .enumerate()
        .filter(|(index, _)| *index != best_index)
        .all(|(_, score)| {
            best_canonical
                .checked_sub(*score)
                .is_some_and(|loss| loss >= edit_loss_floor)
        });
    if separated {
        Nucleotide::ALL[best_index].as_byte()
    } else {
        WILDCARD
    }
}

/// The per-locus unique canonical optimum, when every locus has one.
fn unique_optimum(
    profiles: &[Option<EvidenceProfile>],
    config: &AlignmentConfig,
) -> Option<Vec<u8>> {
    profiles
        .iter()
        .map(|profile| {
            unique_canonical_maximum(substitution_scores(*profile, config))
                .map(|(index, _)| Nucleotide::ALL[index].as_byte())
        })
        .collect()
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
