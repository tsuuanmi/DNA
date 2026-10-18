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

use crate::scoring::{scaled, substitution_scores};
use crate::{AlignmentConfig, MAX_ALIGNMENT_CELLS};
use dna_kernel::model::nucleotide::Nucleotide;
use dna_kernel::read_evidence::EvidenceProfile;

pub(crate) use seeded::{Pruned, align_at_or_above, align_pruned};
pub(crate) use upper_bound::{UpperBoundPlacement, classify};

/// Query locus that seeds never contain: its best substitution is not one
/// unique canonical base, or loses less than one edit unit to the runner-up.
const WILDCARD: u8 = b'N';

/// The reference a read aligns against, with the facts every proof reads
/// from it, derived once for both orientations.
#[derive(Debug)]
pub(crate) struct WorkingReference<'a> {
    /// Linear reference, or a circular reference written twice.
    sequence: &'a str,
    /// Circle length of a doubled circular reference.
    modulo_length: Option<usize>,
    /// Non-canonical symbols one alignment can cover: an alignment consumes
    /// at most one circle of a circular reference.
    ambiguous_bases: usize,
}

impl<'a> WorkingReference<'a> {
    pub(crate) fn new(sequence: &'a str, modulo_length: Option<usize>) -> Self {
        let ambiguous_bases = sequence
            .bytes()
            .take(modulo_length.unwrap_or(sequence.len()))
            .filter(|base| Nucleotide::from_char(char::from(*base)).is_none())
            .count();
        Self {
            sequence,
            modulo_length,
            ambiguous_bases,
        }
    }

    pub(crate) fn sequence(&self) -> &'a str {
        self.sequence
    }

    pub(crate) fn modulo_length(&self) -> Option<usize> {
        self.modulo_length
    }

    /// Length of one copy of the reference.
    fn canonical_length(&self) -> usize {
        self.modulo_length.unwrap_or(self.sequence.len())
    }
}

/// One orientation's alignment problem with its profile bound, derived once
/// and shared by every proof tier.
#[derive(Debug)]
pub(crate) struct Problem<'a> {
    query: &'a str,
    profiles: &'a [Option<EvidenceProfile>],
    reference: &'a WorkingReference<'a>,
    config: &'a AlignmentConfig,
    bound: ProfileBound,
}

impl<'a> Problem<'a> {
    /// The problem with its profile bound, or `None` when no proof applies
    /// and the caller must use full-reference Gotoh.
    pub(crate) fn new(
        query: &'a str,
        profiles: &'a [Option<EvidenceProfile>],
        reference: &'a WorkingReference<'a>,
        config: &'a AlignmentConfig,
    ) -> Option<Self> {
        if !problem_is_provable(query, profiles, reference, config) {
            return None;
        }
        Some(Self {
            query,
            profiles,
            reference,
            config,
            bound: profile_bound(profiles, config)?,
        })
    }
}

/// Score bounds shared by both proof tiers.
#[derive(Debug)]
struct ProfileBound {
    /// Unique profile-optimal base per query locus, or `WILDCARD`.
    optimal_sequence: Vec<u8>,
    /// The strictly best canonical base of every locus, when every locus has
    /// one that no other reference symbol ties.
    unique_optimum: Option<Vec<u8>>,
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
///
/// The negative gap scores that `problem_is_provable` requires make both
/// loss floors positive.
fn profile_bound(
    profiles: &[Option<EvidenceProfile>],
    config: &AlignmentConfig,
) -> Option<ProfileBound> {
    let gap_extension = scaled(config.gap_extension_score);
    let deletion_run_loss = scaled(config.gap_open_score)
        .checked_add(gap_extension)?
        .checked_neg()?;
    let deletion_loss_floor = gap_extension.checked_neg()?;

    let mut loci = Vec::with_capacity(profiles.len());
    let mut insertion_loss = i64::MAX;
    for profile in profiles {
        let scores = substitution_scores(*profile, config);
        let best_score = scores.into_iter().fold(i64::MIN, i64::max);
        if best_score <= gap_extension {
            return None;
        }
        insertion_loss = insertion_loss.min(best_score.checked_sub(gap_extension)?);
        loci.push((scores, best_score));
    }
    let edit_loss_floor = insertion_loss.min(deletion_run_loss);

    let mut optimal_sequence = Vec::with_capacity(loci.len());
    let mut unique_optimum = Some(Vec::with_capacity(loci.len()));
    let mut upper_bound = 0_i64;
    for (scores, best_score) in loci {
        upper_bound = upper_bound.checked_add(best_score)?;
        optimal_sequence.push(canonical_seed_base(scores, edit_loss_floor));
        unique_optimum = unique_optimum.and_then(|mut sequence| {
            sequence.push(unique_canonical_maximum(scores)?);
            Some(sequence)
        });
    }

    Some(ProfileBound {
        optimal_sequence,
        unique_optimum,
        upper_bound,
        edit_loss_floor,
        deletion_loss_floor,
    })
}

/// The locus's seed base: the canonical base whose substitution beats every
/// other canonical base by at least `edit_loss_floor`, or `WILDCARD`.
///
/// A non-canonical reference symbol is not a competitor here; the seeded tier
/// budgets one broken seed per such symbol instead.
fn canonical_seed_base(scores: [i64; 5], edit_loss_floor: i64) -> u8 {
    let [a, c, g, t, _] = scores;
    let canonical = [a, c, g, t];
    let best = a.max(c).max(g).max(t);
    // Bases within one edit unit of the best, the best itself included.
    let mut leaders = Nucleotide::ALL
        .into_iter()
        .zip(canonical)
        .filter(|(_, score)| {
            best.checked_sub(*score)
                .is_none_or(|loss| loss < edit_loss_floor)
        });
    match (leaders.next(), leaders.next()) {
        (Some((base, _)), None) => base.as_byte(),
        _ => WILDCARD,
    }
}

/// The locus's canonical base whose score strictly beats every other
/// reference symbol, non-canonical ones included.
fn unique_canonical_maximum(scores: [i64; 5]) -> Option<u8> {
    let best_score = scores.into_iter().fold(i64::MIN, i64::max);
    let mut best = scores
        .into_iter()
        .enumerate()
        .filter(|(_, score)| *score == best_score);
    let (index, _) = best.next()?;
    if best.next().is_some() {
        return None;
    }
    Nucleotide::ALL.get(index).map(|base| base.as_byte())
}

fn problem_is_provable(
    query: &str,
    profiles: &[Option<EvidenceProfile>],
    reference: &WorkingReference<'_>,
    config: &AlignmentConfig,
) -> bool {
    if query.is_empty()
        || reference.sequence.is_empty()
        || profiles.len() != query.len()
        || config.gap_open_score >= 0
        || config.gap_extension_score >= 0
    {
        return false;
    }
    if let Some(length) = reference.modulo_length
        && (length == 0 || query.len() > length || length > reference.sequence.len())
    {
        return false;
    }

    let Some(rows) = query.len().checked_add(1) else {
        return false;
    };
    let Some(width) = reference.sequence.len().checked_add(1) else {
        return false;
    };
    rows.checked_mul(width)
        .is_some_and(|cells| cells <= MAX_ALIGNMENT_CELLS)
}
