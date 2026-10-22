//! Reference confirmations (ADR-0080): withdrawn where a second sequence
//! may hide, beside a difference, or in a run across an unobserved junction.

use std::collections::{BTreeMap, BTreeSet};

use crate::decide::Decision;
use crate::observe::{Observation, ReadObservations, SiteValue};
use crate::runs::reference_run;
use dna_kernel::model::nucleotide::is_canonical;

/// Withdraws the reference confirmations a second sequence may hide
/// (ADR-0080): at a position where a read of the confirming strand shows a
/// second peak its modality does not explain, and where a read of the other
/// strand does so without its second base standing out, or with a second
/// base that a read of the confirming strand shows standing out too. A
/// second on one strand alone is not the two-strand signature of a second
/// sequence, nor is a second base both strands show that the modality
/// explains on the other strand: the confirming call's own second is
/// explained, so neither strand shows an unexplained second there.
pub(crate) fn hide_confirmations(observed: &mut [ReadObservations]) {
    let hidden = observed
        .iter()
        .enumerate()
        .flat_map(|(index, read)| {
            read.positions
                .iter()
                .filter(|(_, observation)| observation.confirms)
                .filter(|(position, _)| {
                    // A second base standing out on a read of the confirming
                    // strand.
                    let confirming_strand_shows = |base: &char| {
                        observed.iter().any(|third| {
                            third.forward == read.forward
                                && third.seconds.get(position) == Some(base)
                        })
                    };
                    observed.iter().any(|other| {
                        other.unexplained.contains(position)
                            && (other.forward == read.forward
                                || other
                                    .seconds
                                    .get(position)
                                    .is_none_or(confirming_strand_shows))
                    })
                })
                .map(move |(position, _)| (index, *position))
        })
        .collect::<Vec<_>>();
    for (index, position) in hidden {
        if let Some(observation) = observed[index].positions.get_mut(&position) {
            withdraw(observation);
        }
    }
}

/// Withdraws every reference confirmation at or beside a position where any
/// read observes something other than the reference, cleanly or not: another
/// base, a deletion, or an insertion at a junction beside it. A read can be
/// shifted there, so a reference base beside its difference shows nothing
/// about the length the reference keeps.
pub(crate) fn withdraw_beside_observed_differences(
    observed: &mut [ReadObservations],
    bases: &[u8],
) {
    let mut differing = BTreeSet::new();
    for read in observed.iter() {
        for (&position, observation) in &read.positions {
            let differs = match &observation.value {
                SiteValue::Base(base) => {
                    is_canonical(*base)
                        && bases
                            .get(position)
                            .is_some_and(|&reference| char::from(reference) != *base)
                }
                SiteValue::Deletion => true,
                SiteValue::Insertion(_) => false,
            };
            if differs {
                differing.insert(position);
            }
        }
        for (&position, observation) in &read.junctions {
            if observation.inserts() {
                differing.insert(position);
                differing.insert(position + 1);
            }
        }
    }
    withdraw_confirmations(observed, |position| {
        (position.saturating_sub(1)..=position + 1).any(|at| differing.contains(&at))
    });
}

/// Withdraws the confirmations of a position no read observes cleanly
/// otherwise, when it forms a run of one base with a neighbour across a
/// junction no read observes cleanly: a confirmation shows a base, never a
/// length, and an unobserved junction in a run is where a length may change.
pub(crate) fn withdraw_unobserved_run_junctions(
    observed: &mut [ReadObservations],
    bases: &[u8],
    seeded_run: usize,
) {
    // The length of the reference run of one base holding `position`.
    let run_length = |position: usize| {
        reference_run(bases, position).map_or(0, |(first, last)| last + 1 - first)
    };
    let observed_cleanly = |position: usize| clean_call_at(observed, position);
    let junction_clean = |before: usize| {
        observed.iter().any(|read| {
            read.junctions
                .get(&before)
                .is_some_and(|observation| observation.clean)
        })
    };
    let same_base = |left: usize, right: usize| {
        bases.get(left).is_some()
            && bases.get(left) == bases.get(right)
            && run_length(left) >= seeded_run
    };
    let open = confirmed(observed)
        .into_iter()
        .filter(|&position| !observed_cleanly(position))
        .filter(|&position| {
            let left = position
                .checked_sub(1)
                .is_some_and(|left| same_base(left, position) && !junction_clean(left));
            let right = same_base(position, position + 1) && !junction_clean(position);
            left || right
        })
        .collect::<BTreeSet<_>>();
    withdraw_confirmations(observed, |position| open.contains(&position));
}

/// Withdraws every reference confirmation within or beside a decided site
/// that differs from the reference (ADR-0080): there the position stays as
/// the other reads leave it, so the site keeps the neighbours it had. A site
/// written partly as `N` differs too, so a confirmation beside a partly
/// undecided site is withdrawn as well. Whether any was withdrawn.
pub(crate) fn withdraw_near_differences(
    observed: &mut [ReadObservations],
    decided: &BTreeMap<usize, (usize, Decision)>,
    bases: &[u8],
) -> bool {
    let differing = decided
        .iter()
        .filter(|(first, (last, decision))| {
            let reference = (**first..=*last)
                .map(|position| char::from(bases[position]))
                .collect::<String>();
            decision
                .value
                .as_deref()
                .is_some_and(|value| value != reference)
        })
        .map(|(&first, (last, _))| first.saturating_sub(1)..=last + 1)
        .collect::<Vec<_>>();
    withdraw_confirmations(observed, |position| {
        differing.iter().any(|near| near.contains(&position))
    })
}

/// A reference confirmation withdrawn: the observation is unclean, as its
/// call's verdict leaves it.
pub(crate) fn withdraw(observation: &mut Observation) {
    observation.clean = false;
    observation.confirms = false;
}

/// Withdraws every reference confirmation at a position `near` accepts.
/// Whether any was withdrawn.
pub(crate) fn withdraw_confirmations(
    observed: &mut [ReadObservations],
    near: impl Fn(usize) -> bool,
) -> bool {
    let mut withdrawn = false;
    for read in observed {
        for (&position, observation) in &mut read.positions {
            if observation.confirms && near(position) {
                withdraw(observation);
                withdrawn = true;
            }
        }
    }
    withdrawn
}

/// The positions a reference confirmation observes.
fn confirmed(observed: &[ReadObservations]) -> BTreeSet<usize> {
    observed
        .iter()
        .flat_map(|read| {
            read.positions
                .iter()
                .filter(|(_, observation)| observation.confirms)
                .map(|(position, _)| *position)
        })
        .collect()
}

/// Some read observes `position` cleanly by its own call.
fn clean_call_at(observed: &[ReadObservations], position: usize) -> bool {
    observed.iter().any(|read| {
        read.positions
            .get(&position)
            .is_some_and(Observation::clean_call)
    })
}

/// The positions only reference confirmations observe cleanly.
pub(crate) fn confirmed_only(observed: &[ReadObservations]) -> BTreeSet<usize> {
    confirmed(observed)
        .into_iter()
        .filter(|&position| !clean_call_at(observed, position))
        .collect()
}
