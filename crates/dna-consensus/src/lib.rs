//! Sample consensus: adjudicates between placed reads and assembles one
//! consensus sequence per covered reference interval (PROP-0003, ADR-0073).
//!
//! Each read observes the reference positions it covers and the junctions
//! between them. Wherever a read cleanly differs from the reference, the site
//! is grown to whole reference runs, including the runs it joins, and decided
//! as one stretch by its run structure (ADR-0074): the ordered run bases,
//! then each run's length with how it is known. Reads that describe one
//! haplotype through different alignments therefore agree. Everywhere else
//! each position is decided on its own. Only clean observations (unmasked,
//! trusted, no support veto) decide, under decision rule version 1.

pub mod config;
pub mod model;

mod decide;
pub mod decided;
mod insertion;
pub mod length;
pub mod noise;
mod observe;

use std::collections::{BTreeMap, BTreeSet};

use crate::model::consensus::{
    Consensus, ConsensusRun, ConsensusSegment, ConsensusSite, ConsensusSummary, RunEvidence,
    SiteState,
};
use decide::{Decision, decide};
use dna_call::model::called_read::CalledRead;
use dna_call::read_call::CallConfig;
use dna_kernel::model::nucleotide::is_canonical;
use dna_kernel::model::reference::Reference;
use dna_kernel::plugin::{Contract, PluginDescriptor, PluginFamily};
use observe::{Char, Observation, ReadObservations, SiteValue, Trust, observe};

pub use config::{ConsensusConfig, RawConsensusConfig};
pub use insertion::{InsertionView, insertion_views};

/// Builds the consensus of the sample's admitted reads.
///
/// Reads are expected in the sample's deterministic read order; the reads
/// listed at each site follow it.
#[must_use]
pub fn build(
    reads: &[&CalledRead],
    reference: &Reference,
    config: &CallConfig,
    consensus_config: &ConsensusConfig,
) -> Consensus {
    let bases = reference.sequence.as_bytes();
    let mut observed = reads
        .iter()
        .map(|read| {
            observe(
                read,
                &config.variant_calling,
                bases,
                consensus_config.seeded_run,
            )
        })
        .collect::<Vec<_>>();
    let unreliable_only = observe::unreliable_only(reads, &config.variant_calling);
    // Reference confirmations (ADR-0080): none where a second sequence may
    // hide, none within or beside a site that differs from the reference, and
    // the positions they alone decide bound no neighbour's run when called.
    hide_confirmations(&mut observed);
    withdraw_beside_observed_differences(&mut observed, bases);
    withdraw_unobserved_run_junctions(&mut observed, bases, consensus_config.seeded_run);
    let mut decided = decide_positions(&observed, bases, consensus_config.seeded_run);
    // Deciding again can move a differing site beside another confirmation,
    // or leave a reference decision with dissent to set aside: each pass
    // withdraws at least one clean observation, so this ends.
    while withdraw_near_differences(&mut observed, &decided, bases)
        | set_aside_dissent(
            &mut observed,
            &decided,
            bases,
            config.variant_calling.read_end_margin,
            consensus_config.dissent_edge,
        )
    {
        decided = decide_positions(&observed, bases, consensus_config.seeded_run);
    }
    // An insertion into a short repeat that a read of the other strand reads
    // through without is decided by neither strand (SRS-CONS-015).
    insertion::contest_one_strand(&mut decided, &observed, bases, consensus_config.seeded_run);
    two_strand(&mut decided, reads, reference, config, consensus_config);
    let reference_only = confirmed_only(&observed);

    let mut consensus = Consensus {
        segments: Vec::new(),
        sites: Vec::new(),
        summary: ConsensusSummary::default(),
        reference_only,
        unreliable_only,
    };
    for piece in pieces(&decided, consensus_config.maximum_undecided_gap) {
        assemble(
            &piece,
            &decided,
            bases,
            consensus_config.minimum_resolved_bases,
            &mut consensus,
        );
    }
    consensus
}

/// One run of one read over a reference stretch, as the consensus reads it
/// (ADR-0074): its base, its call count, and whether the read shows its
/// length in phase.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReadRun {
    /// The run's base, on the reference strand.
    pub base: char,
    /// The read's calls of the run.
    pub length: usize,
    /// Every call of the run is clean and counts, and the calls just before
    /// and after it are unmasked calls of another base that resolve the
    /// run's ends: the `in_phase` length evidence of a run.
    pub in_phase: bool,
}

/// Each read's runs over the 0-based reference stretch `first..=last`, by
/// the consensus's own reading of the read (its observations, framed by the
/// calls beside the stretch); `None` for a read that does not cover the
/// stretch or reads an unresolved call in it. Reads come in `reads` order.
#[must_use]
pub fn read_runs_over(
    reads: &[&CalledRead],
    reference: &Reference,
    config: &CallConfig,
    consensus_config: &ConsensusConfig,
    first: usize,
    last: usize,
) -> Vec<Option<Vec<ReadRun>>> {
    let bases = reference.sequence.as_bytes();
    reads
        .iter()
        .map(|read| {
            let observed = observe(
                read,
                &config.variant_calling,
                bases,
                consensus_config.seeded_run,
            );
            let (left, chars, right) = observed.span(first, last, bases.len())?;
            if chars.iter().any(|char| char.base == 'N') {
                return None;
            }
            let runs = read_runs(&observed.id, left, &resolved(&chars), right, true);
            Some(
                runs.runs
                    .iter()
                    .map(|run| ReadRun {
                        base: run.base,
                        length: run.length,
                        in_phase: run.clean
                            && run.counted
                            && run.left == End::InPhase
                            && run.right == End::InPhase,
                    })
                    .collect(),
            )
        })
        .collect()
}

/// The runs at which an estimated length gives way to the two-strand rule.
const fn estimated(evidence: RunEvidence) -> bool {
    matches!(
        evidence,
        RunEvidence::AnchoredEnd | RunEvidence::PhaseLoss | RunEvidence::ReferenceFrame
    )
}

/// Decides by the two-strand rule over the reads' tract decodes (SRS-CONS-013,
/// PROP-0012 increment 4) the runs a stretch's decision left as an estimate
/// (an anchored end, a phase loss, or the reference frame), in every pair of
/// runs of one base the stretch splits by a single-base anchor, where its
/// call keeps the reference's run composition. A run the rule decides takes
/// its length, with `two_strand` evidence. A run a read shows in phase, or
/// whose ends two reads confirm, keeps its length whatever the rule says;
/// the review flags their disagreement. The pattern is the stretch's own
/// structure, never a position: the rule names no region.
fn two_strand(
    decided: &mut BTreeMap<usize, (usize, Decision)>,
    reads: &[&CalledRead],
    reference: &Reference,
    config: &CallConfig,
    consensus_config: &ConsensusConfig,
) {
    let bases = reference.sequence.as_bytes();
    for (&start, (end, decision)) in decided.iter_mut() {
        if !decision.runs.iter().any(|run| estimated(run.evidence)) {
            continue;
        }
        let reference_runs = runs_of(bases[start..=*end].iter().map(|&base| char::from(base)));
        if reference_runs.len() != decision.runs.len()
            || reference_runs
                .iter()
                .zip(&decision.runs)
                .any(|((base, _), run)| *base != run.base)
        {
            continue;
        }
        // Each reference run's first position, 1-based.
        let firsts = reference_runs
            .iter()
            .scan(start + 1, |position, (_, length)| {
                let first = *position;
                *position += length;
                Some(first)
            })
            .collect::<Vec<_>>();
        let mut changed = false;
        for index in 0..reference_runs.len().saturating_sub(2) {
            let (repeat, anchor) = (reference_runs[index], reference_runs[index + 1]);
            if reference_runs[index + 2].0 != repeat.0
                || anchor.1 != 1
                || decision.runs[index + 1].length != 1
                || !(estimated(decision.runs[index].evidence)
                    || estimated(decision.runs[index + 2].evidence))
            {
                continue;
            }
            let window = length::Window {
                first: firsts[index],
                last: firsts[index + 2] + reference_runs[index + 2].1 - 1,
                anchor: firsts[index + 1],
                repeat: repeat.0,
                anchor_base: anchor.0,
            };
            let shown = read_runs_over(
                reads,
                reference,
                config,
                consensus_config,
                window.first - 1,
                window.last - 1,
            );
            let views = length::views(&window, reads, &shown);
            for (run, slot) in [(1, index), (2, index + 2)] {
                if !estimated(decision.runs[slot].evidence) {
                    continue;
                }
                let (length::Outcome::Decided(length), _) =
                    length::decide(&views[run - 1], run, consensus_config.length_rule())
                else {
                    continue;
                };
                decision.runs[slot].length = length;
                decision.runs[slot].evidence = RunEvidence::TwoStrand;
                changed = true;
            }
        }
        if changed {
            decision.value = Some(
                decision
                    .runs
                    .iter()
                    .map(|run| run.base.to_string().repeat(run.length))
                    .collect(),
            );
        }
    }
}

/// Decides every stretch by its run structure and every other covered
/// position by its clean observations, then leaves unplaced the length
/// changes nothing beside them places.
fn decide_positions(
    observed: &[ReadObservations],
    bases: &[u8],
    seeded_run: usize,
) -> BTreeMap<usize, (usize, Decision)> {
    let covered = observed
        .iter()
        .flat_map(|read| read.positions.keys().copied())
        .collect::<BTreeSet<_>>();
    // A stretch is decided without reference confirmations, by position as
    // well: a confirmation shows a base, not a length or a composition.
    let mut strict = observed.to_vec();
    for read in &mut strict {
        for observation in read.positions.values_mut() {
            if observation.confirms {
                withdraw(observation);
            }
        }
    }
    let mut decided: BTreeMap<usize, (usize, Decision)> = BTreeMap::new();
    for (first, last) in stretches(&strict, bases, seeded_run) {
        decided.insert(first, (last, decide_stretch(&strict, first, last, bases)));
    }
    let inside = decided
        .iter()
        .flat_map(|(&first, (last, _))| first..=*last)
        .collect::<BTreeSet<_>>();
    for &position in covered.difference(&inside) {
        let observations = observed
            .iter()
            .filter_map(|read| {
                read.positions.get(&position).map(|observation| {
                    let sequence = match &observation.value {
                        SiteValue::Base(base) => base.to_string(),
                        SiteValue::Deletion | SiteValue::Insertion(_) => String::new(),
                    };
                    (read.id.as_str(), observation.clean.then_some(sequence))
                })
            })
            .collect::<Vec<_>>();
        let expected = char::from(bases[position]).to_string();
        decided.insert(position, (position, decide(&observations, &expected)));
    }
    unplace_unbounded_indels(&mut decided, bases.len());
    decided
}

/// Withdraws the reference confirmations a second sequence may hide
/// (ADR-0080): at a position where a read of the confirming strand shows a
/// second peak its modality does not explain, where a read of the other
/// strand does so without its second base standing out, or with a second
/// base that a read of the confirming strand shows standing out too, and
/// where a read of the other strand shows the confirming call's own second
/// base standing out. A second on one strand alone is not the two-strand
/// signature of a second sequence.
fn hide_confirmations(observed: &mut [ReadObservations]) {
    let hidden = observed
        .iter()
        .enumerate()
        .flat_map(|(index, read)| {
            read.positions
                .iter()
                .filter(|(_, observation)| observation.confirms)
                .filter(|(position, _)| {
                    let second = read.seconds.get(position);
                    // A second base standing out on a read of the confirming
                    // strand.
                    let confirming_strand_shows = |base: &char| {
                        observed.iter().any(|third| {
                            third.forward == read.forward
                                && third.seconds.get(position) == Some(base)
                        })
                    };
                    observed.iter().any(|other| {
                        let other_strand = other.forward != read.forward;
                        (other.unexplained.contains(position)
                            && (!other_strand
                                || other
                                    .seconds
                                    .get(position)
                                    .is_none_or(confirming_strand_shows)))
                            || (other_strand
                                && second.is_some()
                                && other.seconds.get(position) == second)
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
fn withdraw_beside_observed_differences(observed: &mut [ReadObservations], bases: &[u8]) {
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
            if matches!(&observation.value, SiteValue::Insertion(inserted) if !inserted.is_empty())
            {
                differing.insert(position);
                differing.insert(position + 1);
            }
        }
    }
    let near = |position: usize| {
        (position.saturating_sub(1)..=position + 1).any(|at| differing.contains(&at))
    };
    for read in observed.iter_mut() {
        for (position, observation) in &mut read.positions {
            if observation.confirms && near(*position) {
                withdraw(observation);
            }
        }
    }
}

/// Withdraws the confirmations of a position no read observes cleanly
/// otherwise, when it forms a run of one base with a neighbour across a
/// junction no read observes cleanly: a confirmation shows a base, never a
/// length, and an unobserved junction in a run is where a length may change.
fn withdraw_unobserved_run_junctions(
    observed: &mut [ReadObservations],
    bases: &[u8],
    seeded_run: usize,
) {
    // The reference run of one base holding `position`.
    let run_length = |position: usize| {
        let Some(&base) = bases.get(position) else {
            return 0;
        };
        let mut first = position;
        while first > 0 && bases[first - 1] == base {
            first -= 1;
        }
        let mut last = position;
        while last + 1 < bases.len() && bases[last + 1] == base {
            last += 1;
        }
        last + 1 - first
    };
    let observed_cleanly = |position: usize| {
        observed.iter().any(|read| {
            read.positions
                .get(&position)
                .is_some_and(|observation| observation.clean && !observation.confirms)
        })
    };
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
    let confirmed = observed
        .iter()
        .flat_map(|read| {
            read.positions
                .iter()
                .filter(|(_, observation)| observation.confirms)
                .map(|(position, _)| *position)
        })
        .collect::<BTreeSet<_>>();
    let open = confirmed
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
    for read in observed.iter_mut() {
        for (position, observation) in &mut read.positions {
            if observation.confirms && open.contains(position) {
                withdraw(observation);
            }
        }
    }
}

/// Withdraws every reference confirmation within or beside a decided site
/// that differs from the reference (ADR-0080): there the position stays as
/// the other reads leave it, so the site keeps the neighbours it had. A site
/// written partly as `N` differs too, so a confirmation beside a partly
/// undecided site is withdrawn as well. Whether any was withdrawn.
fn withdraw_near_differences(
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
    let mut withdrawn = false;
    for read in observed {
        for (position, observation) in &mut read.positions {
            if observation.confirms && differing.iter().any(|near| near.contains(position)) {
                withdraw(observation);
                withdrawn = true;
            }
        }
    }
    withdrawn
}

/// Sets aside the dissent from a reference decision that a read of the other
/// strand contradicts (SRS-CONS-011). At a site or stretch decided as the
/// reference, the reads opposing it are set aside, their observations there
/// withdrawn, when every one of them is on a strand other than a supporting
/// read's that observes every position cleanly, shows no second base
/// standing out there, and lies far enough from its own ends, and either:
///
/// - every opposing read shows a second base standing out there and no
///   second peak its modality leaves unexplained, and the supporting read
///   lies at least `margin` calls from its ends; or
/// - every opposing read differs only by an insertion or deletion within
///   `dissent_edge` calls of its own ends, and the supporting read lies at
///   least `dissent_edge` calls from its ends.
///
/// At a site or stretch called by a strict majority as other than the
/// reference, the first rule alone applies, and only to reads that show the
/// decided runs in order and differ from the decision in a run's length
/// alone.
///
/// A set-aside read shows nothing there and no longer opposes the decision,
/// which keeps its value. Whether any was set aside.
fn set_aside_dissent(
    observed: &mut [ReadObservations],
    decided: &BTreeMap<usize, (usize, Decision)>,
    bases: &[u8],
    margin: usize,
    dissent_edge: usize,
) -> bool {
    let length = bases.len();
    let mut aside = Vec::new();
    for (&first, (last, decision)) in decided {
        let last = *last;
        if decision.opposing.is_empty() || first > last {
            continue;
        }
        let reference = bases[first..=last]
            .iter()
            .map(|&base| char::from(base))
            .collect::<String>();
        let Some(value) = decision.value.as_deref() else {
            continue;
        };
        // Off the reference, only a called decision's dissent in a run's
        // length alone may be set aside.
        let at_reference = value == reference;
        if !at_reference && decision.state != SiteState::Called {
            continue;
        }
        let index = |name: &String| observed.iter().position(|read| &read.id == name);
        let opposing = decision
            .opposing
            .iter()
            .filter_map(index)
            .collect::<Vec<_>>();
        let supporting = decision
            .supporting
            .iter()
            .filter_map(index)
            .collect::<Vec<_>>();
        if opposing.len() != decision.opposing.len() {
            continue;
        }
        let positions = || first..=last;
        // The junctions inside the site and at its two edges.
        let junctions = || std::iter::once((first + length - 1) % length).chain(positions());
        let anchored = |needed: usize| {
            supporting.iter().any(|&read| {
                let read = &observed[read];
                opposing
                    .iter()
                    .all(|&other| observed[other].forward != read.forward)
                    && positions().all(|position| {
                        read.positions
                            .get(&position)
                            .is_some_and(|observation| observation.clean && !observation.confirms)
                            && !read.seconds.contains_key(&position)
                            && read
                                .edges
                                .get(&position)
                                .is_some_and(|&edge| edge >= needed)
                    })
            })
        };
        let explained = |read: &ReadObservations| {
            positions().any(|position| read.seconds.contains_key(&position))
                && !positions().any(|position| read.unexplained.contains(&position))
        };
        let indel_near_end = |read: &ReadObservations| {
            let mut indels = Vec::new();
            for position in positions() {
                match read
                    .positions
                    .get(&position)
                    .map(|observation| &observation.value)
                {
                    Some(SiteValue::Deletion) => indels.push(position),
                    Some(SiteValue::Base(base)) if *base != char::from(bases[position]) => {
                        return false;
                    }
                    _ => {}
                }
            }
            for junction in junctions() {
                if read.junctions.get(&junction).is_some_and(|observation| {
                    matches!(&observation.value, SiteValue::Insertion(inserted) if !inserted.is_empty())
                }) {
                    indels.push(junction);
                }
            }
            !indels.is_empty()
                && indels.iter().all(|position| {
                    read.edges
                        .get(position)
                        .is_some_and(|&edge| edge < dissent_edge)
                })
        };
        // The read shows the decided runs in order, every call resolved, and
        // differs from the decision only in a run's length.
        let composition = |sequence: &str| {
            runs_of(sequence.chars())
                .into_iter()
                .map(|(base, _)| base)
                .collect::<String>()
        };
        let length_only = |read: &ReadObservations| {
            read.span(first, last, length).is_some_and(|(_, chars, _)| {
                let observed = chars.iter().map(|char| char.base).collect::<String>();
                chars.iter().all(|char| is_canonical(char.base))
                    && observed != value
                    && composition(&observed) == composition(value)
            })
        };
        let artefact = opposing.iter().all(|&read| {
            explained(&observed[read]) && (at_reference || length_only(&observed[read]))
        }) && anchored(margin);
        let edge_indel = at_reference
            && opposing.iter().all(|&read| indel_near_end(&observed[read]))
            && anchored(dissent_edge);
        if artefact || edge_indel {
            aside.extend(opposing.into_iter().map(|read| (read, first, last)));
        }
    }
    // A set-aside read shows nothing at the site: neither a base nor, inside
    // a stretch, a run's length or end.
    let mut withdrawn = false;
    let mut set_aside = |observation: &mut Observation| {
        if observation.clean
            || observation
                .chars
                .iter()
                .any(|char| char.trust > Trust::Nothing)
        {
            withdraw(observation);
            for char in &mut observation.chars {
                char.trust = Trust::Nothing;
            }
            withdrawn = true;
        }
    };
    for (read, first, last) in aside {
        let read = &mut observed[read];
        for position in first..=last {
            if let Some(observation) = read.positions.get_mut(&position) {
                set_aside(observation);
            }
        }
        for junction in std::iter::once((first + length - 1) % length).chain(first..=last) {
            if let Some(observation) = read.junctions.get_mut(&junction) {
                set_aside(observation);
            }
        }
    }
    withdrawn
}

/// A reference confirmation withdrawn: the observation is unclean, as its
/// call's verdict leaves it.
fn withdraw(observation: &mut Observation) {
    observation.clean = false;
    observation.confirms = false;
}

/// The positions only reference confirmations observe cleanly.
fn confirmed_only(observed: &[ReadObservations]) -> BTreeSet<usize> {
    let confirmed = observed
        .iter()
        .flat_map(|read| {
            read.positions
                .iter()
                .filter(|(_, observation)| observation.confirms)
                .map(|(position, _)| *position)
        })
        .collect::<BTreeSet<_>>();
    confirmed
        .into_iter()
        .filter(|position| {
            !observed.iter().any(|read| {
                read.positions
                    .get(position)
                    .is_some_and(|observation| observation.clean && !observation.confirms)
            })
        })
        .collect()
}

/// Leaves undecided every decided site that changes the reference length
/// while a position beside it is undecided or uncovered, unless the site was
/// decided by its run structure with every run's length shown by a read
/// (`in_phase` or `anchored_end`): those reads show the run's bounds
/// themselves. Otherwise the bases beside an insertion or deletion place
/// it: inside a repeat they decide which unit changes, so without them the
/// length change is not placed. Undoing one site can leave another one's
/// neighbour undecided, so this repeats until nothing changes. `length` is
/// the reference length; positions wrap.
fn unplace_unbounded_indels(decided: &mut BTreeMap<usize, (usize, Decision)>, length: usize) {
    if length == 0 {
        return;
    }
    // The decided base at `position`, `None` when it is undecided.
    let base_at = |decided: &BTreeMap<usize, (usize, Decision)>, position: usize| {
        let (&first, (last, decision)) = decided.range(..=position).next_back()?;
        if *last < position {
            return None;
        }
        let value = decision.value.as_deref()?;
        // A site of the reference length writes one base per position;
        // otherwise only its two ends are its own.
        let base = if value.len() == last - first + 1 {
            value.chars().nth(position - first)
        } else if position == first {
            value.chars().next()
        } else {
            value.chars().last()
        }?;
        (base != 'N').then_some(base)
    };
    loop {
        let unplaced = decided
            .iter()
            .filter(|(first, (last, decision))| {
                let (first, last) = (**first, *last);
                let shown = !decision.runs.is_empty()
                    && decision.runs.iter().all(|run| {
                        matches!(
                            run.evidence,
                            RunEvidence::InPhase | RunEvidence::AnchoredEnd
                        )
                    });
                !shown
                    && decision
                        .value
                        .as_ref()
                        .is_some_and(|value| value.len() != last - first + 1)
                    && (base_at(decided, (first + length - 1) % length).is_none()
                        || base_at(decided, (last + 1) % length).is_none())
            })
            .map(|(&first, _)| first)
            .collect::<Vec<_>>();
        if unplaced.is_empty() {
            return;
        }
        for first in unplaced {
            if let Some((_, decision)) = decided.get_mut(&first) {
                decision.value = None;
                decision.state = SiteState::Unresolved;
                decision.runs.clear();
            }
        }
    }
}

/// How a read shows one end of a run (ADR-0071, ADR-0074).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum End {
    /// The call beyond the run is unmasked and reads another base.
    InPhase,
    /// The call beyond the run reads another base and anchors the alignment,
    /// but is masked, or one unresolved call precedes it: it shows at most the
    /// dominant end of a mixture.
    Anchored,
    /// The read loses phase at the run: the call beyond it is masked or
    /// unresolved and shows no end. The run is at least as long as its
    /// in-phase calls.
    PhaseLoss,
    /// The read does not show where the run ends.
    Open,
}

impl End {
    /// The read shows where the run ends.
    fn shown(self) -> bool {
        matches!(self, Self::InPhase | Self::Anchored)
    }
}

/// How `boundary`, the call just beyond a run of `run`, shows that end. A
/// masked, unresolved, or unclean call that reads no other base is where the
/// read loses phase; a clean call of the run's base continues the run.
fn end(boundary: Option<Char>, run: char) -> End {
    match boundary {
        Some(char) if char.resolving() && char.base != run => {
            if char.unmasked() && !char.beside_unresolved {
                End::InPhase
            } else {
                End::Anchored
            }
        }
        Some(char) if !char.clean() => End::PhaseLoss,
        None | Some(_) => End::Open,
    }
}

/// One run of a read's bases over a stretch.
#[derive(Debug, Clone)]
struct RunView {
    base: char,
    length: usize,
    /// Every call of the run is clean: the read shows the run's base.
    clean: bool,
    /// Every call of the run is unmasked: the calls count the run.
    counted: bool,
    left: End,
    right: End,
    /// Reference positions of the calls just before and after the run.
    left_at: Option<usize>,
    right_at: Option<usize>,
    /// Bases of the calls just before and after the run, masked or not.
    left_base: Option<char>,
    right_base: Option<char>,
    /// Reference positions of the run's clean calls.
    clean_at: Vec<usize>,
}

impl RunView {
    /// The read shows the run's length: unmasked calls, both ends shown
    /// (ADR-0071).
    fn bounded(&self) -> bool {
        self.counted && self.left.shown() && self.right.shown()
    }

    /// The run's in-phase length when the read shows one end and loses phase
    /// at the other. Within a length mixture it lies near the dominant
    /// length, so it is an estimate, not a bound.
    fn lower_bound(&self) -> Option<usize> {
        let one_sided = (self.left.shown() && self.right == End::PhaseLoss)
            || (self.left == End::PhaseLoss && self.right.shown());
        (self.counted && one_sided).then_some(self.length)
    }
}

/// One read's runs over a stretch, and the clean parts at its two ends.
struct ReadRuns<'a> {
    /// The read's identity, its input's content SHA-256.
    id: &'a str,
    runs: Vec<RunView>,
    /// No call over the stretch is unresolved (`N`).
    complete: bool,
    /// Runs read cleanly from the stretch's start: whole clean runs, then the
    /// run holding the first unclean call when that run starts cleanly.
    prefix: usize,
    /// The same from the stretch's end, counted backwards.
    suffix: usize,
    /// Reference position of the last clean call of the clean start.
    prefix_last: Option<usize>,
    /// Reference position of the first clean call of the clean end.
    suffix_first: Option<usize>,
}

impl ReadRuns<'_> {
    fn bases(runs: &[RunView]) -> String {
        runs.iter().map(|run| run.base).collect()
    }

    /// Every run is read cleanly, and the read loses phase at neither edge of
    /// the stretch, so it shows the stretch's whole composition.
    fn whole(&self) -> bool {
        self.complete
            && self.runs.iter().all(|run| run.clean)
            && self
                .runs
                .first()
                .is_none_or(|run| run.left != End::PhaseLoss)
            && self
                .runs
                .last()
                .is_none_or(|run| run.right != End::PhaseLoss)
    }

    /// The read's runs mapped onto `composition` by run index: all of them
    /// when they have its bases, otherwise those of its clean start and end
    /// that match its first and last runs.
    fn mapped(&self, composition: &str) -> Vec<(usize, RunView)> {
        let bases = composition.chars().collect::<Vec<_>>();
        if self.complete && Self::bases(&self.runs) == composition {
            return self.runs.iter().cloned().enumerate().collect();
        }
        // A partly clean run at the inner edge of a clean start or end shows
        // its base only: it may hold further runs, so neither its length nor
        // its inner end counts.
        let partial = |run: &RunView, inner_left: bool| {
            let mut run = run.clone();
            if !run.clean {
                run.counted = false;
                if inner_left {
                    run.left = End::Open;
                } else {
                    run.right = End::Open;
                }
            }
            run
        };
        // A clean start or end that holds every run of the composition while
        // the read reads further runs over the stretch places a run end of
        // the composition inside the stretch: it does not show this
        // composition, and maps nothing.
        let covers = |edge: usize| edge == bases.len() && self.runs.len() > edge;
        let mut mapped = Vec::new();
        let prefix = &self.runs[..self.prefix];
        if prefix.len() <= bases.len()
            && !covers(prefix.len())
            && prefix
                .iter()
                .zip(&bases)
                .all(|(run, base)| run.base == *base)
        {
            mapped.extend(prefix.iter().map(|run| partial(run, false)).enumerate());
        }
        let suffix = &self.runs[self.runs.len() - self.suffix..];
        if let Some(offset) = bases.len().checked_sub(suffix.len())
            && !covers(suffix.len())
            && suffix
                .iter()
                .zip(&bases[offset..])
                .all(|(run, base)| run.base == *base)
        {
            for (index, run) in suffix.iter().enumerate() {
                if !mapped.iter().any(|(mapped, _)| *mapped == offset + index) {
                    mapped.push((offset + index, partial(run, true)));
                }
            }
        }
        mapped
    }
}

/// A read's runs over a stretch from its framed observation.
fn read_runs<'a>(
    id: &'a str,
    left: Option<Char>,
    chars: &[Char],
    right: Option<Char>,
    complete: bool,
) -> ReadRuns<'a> {
    let groups = groups_of(chars);
    let runs = groups
        .iter()
        .enumerate()
        .map(|(index, group)| {
            let base = group[0].base;
            let before = index
                .checked_sub(1)
                .map_or(left, |previous| groups[previous].last().copied());
            let after = groups
                .get(index + 1)
                .map_or(right, |next| next.first().copied());
            RunView {
                base,
                length: group.len(),
                clean: group.iter().all(|char| char.clean()),
                counted: group.iter().all(|char| char.unmasked()),
                left: end(before, base),
                right: end(after, base),
                left_at: before.and_then(|char| char.position),
                right_at: after.and_then(|char| char.position),
                left_base: before.map(|char| char.base),
                right_base: after.map(|char| char.base),
                clean_at: group
                    .iter()
                    .filter(|char| char.clean())
                    .filter_map(|char| char.position)
                    .collect(),
            }
        })
        .collect::<Vec<_>>();
    let clean_part = |groups: &mut dyn Iterator<Item = &Vec<Char>>| {
        let mut count = 0;
        for group in groups {
            if group.iter().all(|char| char.clean()) {
                count += 1;
            } else {
                break;
            }
        }
        count
    };
    let mut prefix = clean_part(&mut groups.iter());
    if groups.get(prefix).is_some_and(|group| group[0].clean()) {
        prefix += 1;
    }
    let mut suffix = clean_part(&mut groups.iter().rev());
    if suffix < groups.len()
        && groups[groups.len() - 1 - suffix]
            .last()
            .is_some_and(|char| char.clean())
    {
        suffix += 1;
    }
    let prefix_last = chars
        .iter()
        .take_while(|char| char.clean())
        .filter_map(|char| char.position)
        .last();
    let suffix_first = chars
        .iter()
        .rev()
        .take_while(|char| char.clean())
        .filter_map(|char| char.position)
        .last();
    // A clean start or end anchors the stretch's edge only where the read
    // does not lose phase there.
    if runs.first().is_some_and(|run| run.left == End::PhaseLoss) {
        prefix = 0;
    }
    if runs.last().is_some_and(|run| run.right == End::PhaseLoss) {
        suffix = 0;
    }
    ReadRuns {
        id,
        runs,
        complete,
        prefix: prefix.min(groups.len()),
        suffix: suffix.min(groups.len()),
        prefix_last,
        suffix_first,
    }
}

/// The run bases of one read's clean start joined to another read's clean
/// end, when they overlap: the clean calls of both reach a common reference
/// position, and one way only matches the start's last runs to the end's
/// first runs. `None` otherwise.
fn join(start: &ReadRuns<'_>, finish: &ReadRuns<'_>) -> Option<String> {
    let head = &start.runs[..start.prefix];
    let tail = &finish.runs[finish.runs.len() - finish.suffix..];
    let (Some(last), Some(first)) = (start.prefix_last, finish.suffix_first) else {
        return None;
    };
    if head.is_empty() || tail.is_empty() || last < first {
        return None;
    }
    let mut joined = None;
    for offset in 0..head.len() {
        let overlap = head.len() - offset;
        if overlap > tail.len()
            || head[offset..]
                .iter()
                .zip(tail)
                .any(|(left, right)| left.base != right.base)
        {
            continue;
        }
        if joined.is_some() {
            return None;
        }
        joined = Some(format!(
            "{}{}",
            ReadRuns::bases(&head[..offset]),
            ReadRuns::bases(tail)
        ));
    }
    joined
}

/// Demotes to `PhaseLoss` every anchored run end that a read places
/// elsewhere than where reads showing that end in phase place it: in-phase
/// evidence outweighs a masked call's (ADR-0074).
fn reconcile(mapped: &mut [(&str, Vec<(usize, RunView)>)], runs: usize) {
    type Side = fn(&mut RunView) -> (&mut End, Option<usize>);
    let sides: [Side; 2] = [
        |view| (&mut view.left, view.left_at),
        |view| (&mut view.right, view.right_at),
    ];
    for index in 0..runs {
        for side in sides {
            let mut in_phase = BTreeSet::new();
            for (_, views) in mapped.iter_mut() {
                for (_, view) in views.iter_mut().filter(|(mapped, _)| *mapped == index) {
                    let (end, at) = side(view);
                    if *end == End::InPhase {
                        in_phase.insert(at);
                    }
                }
            }
            let [located] = in_phase.into_iter().collect::<Vec<_>>()[..] else {
                continue;
            };
            for (_, views) in mapped.iter_mut() {
                for (_, view) in views.iter_mut().filter(|(mapped, _)| *mapped == index) {
                    let (end, at) = side(view);
                    if *end == End::Anchored && at != located {
                        *end = End::PhaseLoss;
                    }
                }
            }
        }
    }
}

/// The reference-frame length of run `index`, which no read shows both ends
/// of: the reference positions between the call before it, located by every
/// read that shows its start, and the call after it, located by every read
/// that shows its end, when each of those positions is read cleanly by some
/// read (ADR-0074). `None` when an end is shown by no read, the reads place an
/// end differently, a position between is not read cleanly, or no position is
/// left between.
fn reference_frame(mapped: &[(&str, Vec<(usize, RunView)>)], index: usize) -> Option<usize> {
    let views = mapped
        .iter()
        .flat_map(|(_, views)| views.iter())
        .filter(|(mapped, _)| *mapped == index)
        .map(|(_, view)| view)
        .collect::<Vec<_>>();
    let located = |side: fn(&RunView) -> (End, Option<usize>)| {
        let mut at = views.iter().filter_map(|view| {
            let (end, at) = side(view);
            end.shown().then_some(at)
        });
        let first = at.next()??;
        at.all(|other| other == Some(first)).then_some(first)
    };
    let before = located(|view| (view.left, view.left_at))?;
    let after = located(|view| (view.right, view.right_at))?;
    let read = views
        .iter()
        .flat_map(|view| view.clean_at.iter().copied())
        .collect::<BTreeSet<_>>();
    let span = after.checked_sub(before)?.checked_sub(1)?;
    (span >= 1 && (before + 1..after).all(|position| read.contains(&position))).then_some(span)
}

/// Decides one stretch by its run structure (ADR-0074).
///
/// 1. **Composition**, the ordered run bases: voted from reads that read
///    every run cleanly; otherwise the one composition that joins one read's
///    clean start to another read's clean end; otherwise the reference's, when
///    a read shows its run structure.
/// 2. **Lengths**: each run's length is voted from reads that show both of
///    its ends (ADR-0071), labelled `in_phase` or `anchored_end`.
/// 3. **Reference frame**: when exactly one run has no length, it takes the
///    reference positions between its ends as located by the reads that show
///    them, labelled `reference_frame` (see [`reference_frame`]).
///
/// A stretch left undecided is decided position by position.
fn decide_stretch(
    observed: &[ReadObservations],
    first: usize,
    last: usize,
    bases: &[u8],
) -> Decision {
    let length = bases.len();
    let expected = (first..=last)
        .map(|position| char::from(bases[position]))
        .collect::<String>();
    let reference_runs = runs_of(expected.chars());
    let reference_composition = reference_runs.iter().map(|run| run.0).collect::<String>();
    let reads = observed
        .iter()
        .filter_map(|read| {
            read.span(first, last, length).map(|(left, chars, right)| {
                let complete = chars.iter().all(|char| char.base != 'N');
                read_runs(read.id.as_str(), left, &resolved(&chars), right, complete)
            })
        })
        .collect::<Vec<_>>();
    let names = || reads.iter().map(|read| read.id);
    let fallback = || by_position(observed, first, last, bases, names());
    let whole = reads
        .iter()
        .map(|read| (read.id, read.whole().then(|| ReadRuns::bases(&read.runs))))
        .collect::<Vec<_>>();
    let composition = if whole.iter().any(|(_, bases)| bases.is_some()) {
        let composition = decide(&whole, &reference_composition);
        // A read that reads only an edge of the stretch cleanly, covering the
        // stretch or not, still shows the bases there: when their runs
        // contradict the composition's first or last runs, the read opposes
        // it, with an observation no other read shares.
        let contradicting = composition.value.as_deref().map_or_else(Vec::new, |value| {
            observed
                .iter()
                .filter(|read| {
                    !whole
                        .iter()
                        .any(|(name, bases)| *name == read.id && bases.is_some())
                })
                .filter(|read| {
                    let edge = |from_end| {
                        runs_of(read.clean_edge(first, last, length, from_end).chars())
                            .into_iter()
                            .map(|(base, _)| base)
                            .collect::<String>()
                    };
                    let (prefix, suffix) = (edge(false), edge(true));
                    (!prefix.is_empty() && !value.starts_with(&prefix))
                        || (!suffix.is_empty() && !value.ends_with(&suffix))
                })
                .map(|read| read.id.as_str())
                .collect::<Vec<_>>()
        });
        if contradicting.is_empty() {
            composition
        } else {
            let mut observations = whole
                .iter()
                .filter(|(name, _)| !contradicting.contains(name))
                .cloned()
                .collect::<Vec<_>>();
            observations.extend(
                contradicting
                    .iter()
                    .map(|name| (*name, Some(format!("~{name}")))),
            );
            decide(&observations, &reference_composition)
        }
    } else {
        let mut joined: BTreeMap<String, BTreeSet<&str>> = BTreeMap::new();
        for start in &reads {
            for finish in &reads {
                if let Some(bases) = join(start, finish) {
                    joined
                        .entry(bases)
                        .or_default()
                        .extend([start.id, finish.id]);
                }
            }
        }
        let structured = reads
            .iter()
            .any(|read| read.complete && ReadRuns::bases(&read.runs) == reference_composition);
        match joined.len() {
            1 => {
                let (bases, supporting) = joined.into_iter().next().unwrap_or_default();
                let observations = reads
                    .iter()
                    .map(|read| (read.id, supporting.contains(read.id).then(|| bases.clone())))
                    .collect::<Vec<_>>();
                decide(&observations, &bases)
            }
            0 if structured => Decision {
                value: Some(reference_composition.clone()),
                state: SiteState::Called,
                supporting: Vec::new(),
                opposing: Vec::new(),
                uninformative: Vec::new(),
                runs: Vec::new(),
            },
            _ => return fallback(),
        }
    };
    let Some(runs) = composition.value.clone() else {
        return fallback();
    };
    if runs.is_empty() {
        // Every position deleted: there is no run to measure.
        return composition;
    }
    let at_reference = runs == reference_composition;

    let mut mapped = reads
        .iter()
        .map(|read| (read.id, read.mapped(&runs)))
        .collect::<Vec<_>>();
    let count = runs.chars().count();
    reconcile(&mut mapped, count);
    // Reads that show every run's length vote on the whole sequence first, so
    // that run lengths read through different run boundaries never mix; a tie
    // between them leaves the stretch to the position-by-position decision.
    let lengths = |views: &[(usize, RunView)]| {
        (0..count)
            .map(|index| {
                let view = &views.iter().find(|(mapped, _)| *mapped == index)?.1;
                view.bounded()
                    .then(|| view.base.to_string().repeat(view.length))
            })
            .collect::<Option<String>>()
    };
    let sequences = mapped
        .iter()
        .map(|(name, views)| (*name, lengths(views)))
        .collect::<Vec<_>>();
    if sequences.iter().any(|(_, sequence)| sequence.is_some()) {
        let mut decision = decide(&sequences, &expected);
        let Some(value) = decision.value.as_deref() else {
            return fallback();
        };
        decision.runs = runs_of(value.chars())
            .into_iter()
            .enumerate()
            .map(|(index, (base, length))| {
                let in_phase = mapped.iter().any(|(name, views)| {
                    decision.supporting.iter().any(|read| read == name)
                        && views.iter().any(|(mapped, view)| {
                            *mapped == index
                                && view.left == End::InPhase
                                && view.right == End::InPhase
                        })
                });
                ConsensusRun {
                    base,
                    length,
                    evidence: if in_phase {
                        RunEvidence::InPhase
                    } else {
                        RunEvidence::AnchoredEnd
                    },
                }
            })
            .collect();
        return with_composition(decision, &composition, &names().collect::<Vec<_>>());
    }
    let mut decisions = runs
        .chars()
        .enumerate()
        .map(|(index, base)| {
            let observations = mapped
                .iter()
                .filter_map(|(name, views)| {
                    let view = &views.iter().find(|(mapped, _)| *mapped == index)?.1;
                    Some((
                        *name,
                        view.bounded().then(|| base.to_string().repeat(view.length)),
                    ))
                })
                .collect::<Vec<_>>();
            let reference = if at_reference {
                base.to_string().repeat(reference_runs[index].1)
            } else {
                String::new()
            };
            let decision = decide(&observations, &reference);
            let in_phase = mapped.iter().any(|(name, views)| {
                decision.supporting.iter().any(|read| read == name)
                    && views.iter().any(|(mapped, view)| {
                        *mapped == index && view.left == End::InPhase && view.right == End::InPhase
                    })
            });
            let evidence = if in_phase {
                RunEvidence::InPhase
            } else {
                RunEvidence::AnchoredEnd
            };
            (base, decision, evidence)
        })
        .collect::<Vec<_>>();

    // A run no read shows the length of takes the length two reads count in
    // phase from its opposite ends, every call of the run unmasked, when they
    // agree; otherwise the reference frame between its observed ends; without
    // one, the longest run read in phase before a read loses phase.
    for (index, (base, decision, evidence)) in decisions.iter_mut().enumerate() {
        if decision.value.is_some() || decision.state != SiteState::Unresolved {
            continue;
        }
        let entries = |from_left: bool| {
            mapped
                .iter()
                .filter_map(|(name, views)| {
                    let view = &views.iter().find(|(mapped, _)| *mapped == index)?.1;
                    let (entry, exit, beyond) = if from_left {
                        (view.left, view.right, view.right_base)
                    } else {
                        (view.right, view.left, view.left_base)
                    };
                    let other = beyond.is_some_and(|beyond| beyond != view.base);
                    // The read loses phase at the exit call, which does not
                    // read the run's base: it shows no longer run.
                    (view.counted && entry == End::InPhase && exit == End::PhaseLoss && other)
                        .then_some((*name, view.length))
                })
                .collect::<Vec<_>>()
        };
        let (from_left, from_right) = (entries(true), entries(false));
        let counts = from_left
            .iter()
            .chain(&from_right)
            .map(|(_, length)| *length)
            .collect::<BTreeSet<_>>();
        if !from_left.is_empty() && !from_right.is_empty() && counts.len() == 1 {
            let length = counts.into_iter().next().unwrap_or_default();
            decision.value = Some(base.to_string().repeat(length));
            decision.state = SiteState::Called;
            decision.supporting = from_left
                .iter()
                .chain(&from_right)
                .map(|(name, _)| (*name).to_owned())
                .collect();
            *evidence = RunEvidence::ConfirmedEnd;
            continue;
        }
        let bounds = mapped
            .iter()
            .filter_map(|(name, views)| {
                let view = &views.iter().find(|(mapped, _)| *mapped == index)?.1;
                Some((*name, view.lower_bound()?))
            })
            .collect::<Vec<_>>();
        let bound = bounds.iter().map(|(_, length)| *length).max();
        let (length, source) = match (reference_frame(&mapped, index), bound) {
            (Some(frame), _) => (frame, RunEvidence::ReferenceFrame),
            (None, Some(bound)) => (bound, RunEvidence::PhaseLoss),
            (None, None) => continue,
        };
        decision.value = Some(base.to_string().repeat(length));
        decision.state = SiteState::Called;
        decision.supporting = bounds.iter().map(|(name, _)| (*name).to_owned()).collect();
        *evidence = source;
    }
    if decisions
        .iter()
        .any(|(_, decision, _)| decision.value.is_none())
    {
        if !at_reference {
            return fallback();
        }
        let combined = combine(
            &decisions
                .iter()
                .enumerate()
                .map(|(index, (base, decision, _))| {
                    (*base, reference_runs[index].1, decision.clone())
                })
                .collect::<Vec<_>>(),
            names(),
        );
        return if combined.state == SiteState::Unresolved {
            fallback()
        } else {
            combined
        };
    }
    let mut combined = combine(
        &decisions
            .iter()
            .map(|(base, decision, _)| (*base, 0, decision.clone()))
            .collect::<Vec<_>>(),
        names(),
    );
    combined.runs = decisions
        .iter()
        .map(|(base, decision, evidence)| ConsensusRun {
            base: *base,
            length: decision.value.as_ref().map_or(0, String::len),
            evidence: *evidence,
        })
        .collect();
    with_composition(combined, &composition, &names().collect::<Vec<_>>())
}

/// Adds a stretch's composition decision to its length decision: reads that
/// disagree on the composition oppose the stretch, reads that only showed the
/// composition support it, and a contested composition contests the stretch.
/// Reads are listed in `order`.
fn with_composition(mut decision: Decision, composition: &Decision, order: &[&str]) -> Decision {
    for read in &composition.opposing {
        if !decision.opposing.contains(read) {
            decision.opposing.push(read.clone());
        }
    }
    for read in &composition.supporting {
        if !decision.supporting.contains(read) && !decision.opposing.contains(read) {
            decision.supporting.push(read.clone());
        }
    }
    decision
        .supporting
        .retain(|read| !composition.opposing.contains(read));
    if composition.state == SiteState::Contested {
        decision.state = SiteState::Contested;
    }
    let rank = |read: &String| order.iter().position(|name| name == read);
    decision.supporting.sort_by_key(rank);
    decision.opposing.sort_by_key(rank);
    let Decision {
        supporting,
        opposing,
        uninformative,
        ..
    } = &mut decision;
    uninformative.retain(|read| !supporting.contains(read) && !opposing.contains(read));
    decision
}

/// Decides a stretch position by position, each read without an insertion or
/// deletion inside it contributing the base its alignment places on each
/// position, when it is clean; reads with an indel there are uninformative.
fn by_position<'a>(
    observed: &[ReadObservations],
    first: usize,
    last: usize,
    bases: &[u8],
    covering: impl Iterator<Item = &'a str>,
) -> Decision {
    let length = bases.len();
    // Only clean indels shift a read's bases off their positions.
    let unchanged = |read: &ReadObservations| {
        (first..=last).all(|position| {
            read.positions.get(&position).is_none_or(|observation| {
                observation.value != SiteValue::Deletion || !observation.clean
            })
        }) && std::iter::once((first + length - 1) % length)
            .chain(first..=last)
            .all(|position| {
                read.junctions
                    .get(&position)
                    .is_none_or(|observation| observation.chars.iter().all(|char| !char.clean()))
            })
    };
    let decisions = (first..=last)
        .map(|position| {
            let observations = observed
                .iter()
                .filter_map(|read| {
                    read.positions.get(&position).map(|observation| {
                        let call = match observation.value {
                            SiteValue::Base(base) if observation.clean && unchanged(read) => {
                                Some(base.to_string())
                            }
                            _ => None,
                        };
                        (read.id.as_str(), call)
                    })
                })
                .collect::<Vec<_>>();
            let base = char::from(bases[position]);
            (base, 1, decide(&observations, &base.to_string()))
        })
        .collect::<Vec<_>>();
    combine(&decisions, covering)
}

/// One stretch decision from its runs' decisions: an undecided run is written
/// as `N` over its reference length; the stretch is contested when a run is,
/// and unresolved only when every run is.
fn combine<'a>(
    decisions: &[(char, usize, Decision)],
    covering: impl Iterator<Item = &'a str>,
) -> Decision {
    let value = decisions
        .iter()
        .any(|(_, _, decision)| decision.value.is_some())
        .then(|| {
            decisions
                .iter()
                .map(|(_, reference_length, decision)| {
                    decision
                        .value
                        .clone()
                        .unwrap_or_else(|| "N".repeat(*reference_length))
                })
                .collect::<String>()
        });
    let state = if decisions
        .iter()
        .any(|(_, _, decision)| decision.state == SiteState::Contested)
    {
        SiteState::Contested
    } else if value.is_none() {
        SiteState::Unresolved
    } else {
        SiteState::Called
    };
    let opposing = decisions
        .iter()
        .flat_map(|(_, _, decision)| decision.opposing.iter().cloned())
        .collect::<BTreeSet<_>>();
    let supporting = decisions
        .iter()
        .flat_map(|(_, _, decision)| decision.supporting.iter().cloned())
        .filter(|read| !opposing.contains(read))
        .collect::<BTreeSet<_>>();
    let ordered = |set: &BTreeSet<String>, reads: &[&str]| {
        reads
            .iter()
            .filter(|read| set.contains(**read))
            .map(|read| (*read).to_owned())
            .collect::<Vec<_>>()
    };
    let reads = covering.collect::<Vec<_>>();
    let uninformative = reads
        .iter()
        .filter(|read| !supporting.contains(**read) && !opposing.contains(**read))
        .map(|read| (*read).to_owned())
        .collect();
    Decision {
        value,
        state,
        supporting: ordered(&supporting, &reads),
        opposing: ordered(&opposing, &reads),
        uninformative,
        runs: Vec::new(),
    }
}

/// The observed bases without unresolved calls (`N`). The bases next to a
/// dropped call can neither be clean, nor count towards a run, nor show where
/// a run ends, because the unresolved call may extend either run.
///
/// The exception is a lone unresolved call between two different bases of
/// which exactly one is unmasked: the read loses phase there (ADR-0074). The
/// unmasked base still counts towards its run, the masked one still shows
/// that run's end, and both are marked `beside_unresolved`, so the run's
/// length is known only as `anchored_end`.
fn resolved(chars: &[Char]) -> Vec<Char> {
    let lone_exit = |index: usize| {
        let (Some(before), Some(after)) = (
            index.checked_sub(1).and_then(|before| chars.get(before)),
            chars.get(index + 1),
        ) else {
            return false;
        };
        before.base != 'N'
            && after.base != 'N'
            && before.base != after.base
            && before.unmasked() != after.unmasked()
    };
    let mut kept: Vec<Char> = Vec::with_capacity(chars.len());
    let mut after_unresolved = None;
    for (index, char) in chars.iter().enumerate() {
        let mut char = *char;
        if char.base == 'N' {
            let exit = lone_exit(index);
            if let Some(previous) = kept.last_mut() {
                taint(previous, exit);
            }
            after_unresolved = Some(exit);
            continue;
        }
        if let Some(exit) = after_unresolved.take() {
            taint(&mut char, exit);
        }
        kept.push(char);
    }
    kept
}

/// Marks a base next to a dropped unresolved call: at a phase-loss exit it
/// keeps its count and run end, otherwise it loses both.
fn taint(char: &mut Char, exit: bool) {
    if exit {
        char.trust = char.trust.min(Trust::Unmasked);
        char.beside_unresolved = true;
    } else {
        char.trust = Trust::Nothing;
    }
}

/// Runs of equal bases as `(base, length)`.
fn runs_of(bases: impl Iterator<Item = char>) -> Vec<(char, usize)> {
    let mut runs: Vec<(char, usize)> = Vec::new();
    for base in bases {
        match runs.last_mut() {
            Some((previous, count)) if *previous == base => *count += 1,
            _ => runs.push((base, 1)),
        }
    }
    runs
}

/// Observed bases grouped into runs of equal bases.
fn groups_of(chars: &[Char]) -> Vec<Vec<Char>> {
    let mut groups: Vec<Vec<Char>> = Vec::new();
    for &char in chars {
        match groups.last_mut() {
            Some(group)
                if group
                    .last()
                    .is_some_and(|previous| previous.base == char.base) =>
            {
                group.push(char);
            }
            _ => groups.push(vec![char]),
        }
    }
    groups
}

/// Reference intervals where any read cleanly differs from the reference (an
/// unclean difference alone decides nothing and does not widen a stretch),
/// together with the neighbouring reference runs a substituted base or a
/// deletion joins, and every homopolymer of at least `seeded_run` bases that
/// no read shows cleanly with both ends while some read differs on it: each
/// grown to whole reference runs on both sides, overlapping or touching
/// intervals merged. A stretch never crosses the reference origin. Shorter
/// runs hold few length variants, and their disputes are mostly mixed signal
/// at a heteroplasmic base; `consensus.seeded_run` matches the human-mtDNA
/// profile's conformance minimum run length, which the consensus does not
/// read.
fn stretches(
    observed: &[ReadObservations],
    bases: &[u8],
    seeded_run: usize,
) -> Vec<(usize, usize)> {
    let length = bases.len();
    let mut seeds = BTreeSet::new();
    for read in observed {
        for (&position, observation) in &read.positions {
            if observation.clean
                && observation.value != SiteValue::Base(char::from(bases[position]))
            {
                seeds.insert(position);
                // A substituted base that equals a neighbouring reference
                // base joins that neighbour's run (ADR-0074).
                if let SiteValue::Base(base) = observation.value {
                    let neighbours = [position.checked_sub(1), Some(position + 1)];
                    for neighbour in neighbours.into_iter().flatten() {
                        if neighbour < length && char::from(bases[neighbour]) == base {
                            seeds.insert(neighbour);
                        }
                    }
                }
            }
        }
        // A clean deletion that brings two runs of one base together joins
        // them (ADR-0074).
        let deleted = read
            .positions
            .iter()
            .filter(|(_, observation)| {
                observation.clean && observation.value == SiteValue::Deletion
            })
            .map(|(&position, _)| position)
            .collect::<Vec<_>>();
        for block in deleted.chunk_by(|left, right| left + 1 == *right) {
            if let (Some(&first), Some(&last)) = (block.first(), block.last())
                && let Some(before) = first.checked_sub(1)
                && last + 1 < length
                && bases[before] == bases[last + 1]
            {
                seeds.extend([before, last + 1]);
            }
        }
        for (&position, observation) in &read.junctions {
            if observation.chars.iter().any(|char| char.clean()) && position + 1 < length {
                seeds.insert(position);
                seeds.insert(position + 1);
            }
        }
    }
    // A homopolymer of at least `seeded_run` bases that no read shows whole
    // and cleanly, with both of its ends, while some read differs on it,
    // cleanly or not, is decided by its run structure: decided position by
    // position, it would assert a length no read shows (ADR-0071).
    let mut first = 1;
    while first + 1 < length {
        let mut last = first;
        while last + 1 < length && bases[last + 1] == bases[first] {
            last += 1;
        }
        if last + 1 - first >= seeded_run && last + 1 < length && !seeds.contains(&first) {
            let bounded = observed.iter().any(|read| {
                (first - 1..=last + 1).all(|position| {
                    read.positions
                        .get(&position)
                        .is_some_and(|observation| observation.clean)
                })
            });
            let differs = observed.iter().any(|read| {
                (first - 1..=last + 1).any(|position| {
                    read.positions.get(&position).is_some_and(|observation| {
                        observation.value != SiteValue::Base(char::from(bases[position]))
                    })
                }) || (first - 1..=last).any(|position| {
                    read.junctions.get(&position).is_some_and(|observation| {
                        observation.value != SiteValue::Insertion(String::new())
                    })
                })
            });
            if !bounded && differs {
                seeds.insert(first);
            }
        }
        first = last + 1;
    }
    let run_start = |mut position: usize| {
        while position > 0 && bases[position - 1] == bases[position] {
            position -= 1;
        }
        position
    };
    let run_end = |mut position: usize| {
        while position + 1 < length && bases[position + 1] == bases[position] {
            position += 1;
        }
        position
    };
    let mut merged: Vec<(usize, usize)> = Vec::new();
    for seed in seeds {
        let (first, last) = (run_start(seed), run_end(seed));
        match merged.last_mut() {
            Some((_, previous)) if first <= *previous + 1 => *previous = (*previous).max(last),
            _ => merged.push((first, last)),
        }
    }
    merged
}

/// Splits the decided sites into pieces of contiguous positions, from one
/// cleanly decided site to the next: a gap of at most `maximum_gap` undecided
/// positions stays inside a piece as `N`, while a longer gap or an uncovered
/// position splits it. Each piece lists its sites by first position.
fn pieces(decided: &BTreeMap<usize, (usize, Decision)>, maximum_gap: usize) -> Vec<Vec<usize>> {
    let undecided = |first: &usize| decided[first].1.state == SiteState::Unresolved;
    let close = |current: &mut Vec<usize>, pieces: &mut Vec<Vec<usize>>| {
        while current.last().is_some_and(undecided) {
            current.pop();
        }
        if !current.is_empty() {
            pieces.push(std::mem::take(current));
        }
    };
    let mut pieces = Vec::new();
    let mut current = Vec::new();
    let mut previous_last: Option<usize> = None;
    let mut observed_last: Option<usize> = None;
    for (&first, (last, decision)) in decided {
        let observed = decision.state != SiteState::Unresolved;
        let contiguous = previous_last.is_some_and(|end| end + 1 == first);
        let distant = observed && observed_last.is_some_and(|end| first - end - 1 > maximum_gap);
        if !contiguous || distant {
            close(&mut current, &mut pieces);
            current.clear();
            observed_last = None;
        }
        if observed {
            observed_last = Some(*last);
        }
        if observed_last.is_some() {
            current.push(first);
        }
        previous_last = Some(*last);
    }
    close(&mut current, &mut pieces);
    pieces
}

/// Adds one piece as a segment, with its sites, when its sequence holds at
/// least `minimum` resolved bases.
fn assemble(
    piece: &[usize],
    decided: &BTreeMap<usize, (usize, Decision)>,
    bases: &[u8],
    minimum: usize,
    consensus: &mut Consensus,
) {
    let mut sequence = String::new();
    let mut sites = Vec::new();
    let mut summary = ConsensusSummary::default();
    for first in piece {
        let (last, decision) = &decided[first];
        let span = last - first + 1;
        if let Some(value) = &decision.value {
            sequence.push_str(value);
            summary.unresolved_positions += value.chars().filter(|base| *base == 'N').count();
        } else {
            sequence.extend(std::iter::repeat_n('N', span));
            summary.unresolved_positions += span;
        }
        for run in &decision.runs {
            match run.evidence {
                RunEvidence::ReferenceFrame => summary.reference_frame_runs += 1,
                RunEvidence::PhaseLoss => summary.phase_loss_runs += 1,
                RunEvidence::InPhase
                | RunEvidence::AnchoredEnd
                | RunEvidence::ConfirmedEnd
                | RunEvidence::TwoStrand => {}
            }
        }
        match decision.state {
            SiteState::Called => summary.called_positions += span,
            SiteState::Contested => summary.contested_sites += 1,
            SiteState::Unresolved => {}
        }
        let reference = String::from_utf8_lossy(&bases[*first..=*last]).into_owned();
        // A run length that is an estimate keeps its label even where it
        // equals the reference.
        let estimated = decision.runs.iter().any(|run| {
            matches!(
                run.evidence,
                RunEvidence::PhaseLoss | RunEvidence::ReferenceFrame
            )
        });
        if decision.value.as_deref() != Some(reference.as_str())
            || decision.state != SiteState::Called
            || !decision.opposing.is_empty()
            || estimated
        {
            sites.push(ConsensusSite {
                first_0based: *first,
                last_0based: *last,
                reference,
                call: decision.value.clone(),
                state: decision.state,
                supporting: decision.supporting.clone(),
                opposing: decision.opposing.clone(),
                uninformative: decision.uninformative.clone(),
                runs: decision.runs.clone(),
            });
        }
    }
    // The core cannot place a piece with fewer resolved bases.
    if sequence.chars().filter(|base| *base != 'N').count() < minimum {
        return;
    }
    let (Some(first), Some(last)) = (piece.first(), piece.last()) else {
        return;
    };
    consensus.segments.push(ConsensusSegment {
        start_0based: *first,
        end_0based: decided[last].0,
        sequence,
    });
    consensus.sites.extend(sites);
    consensus.summary.called_positions += summary.called_positions;
    consensus.summary.contested_sites += summary.contested_sites;
    consensus.summary.unresolved_positions += summary.unresolved_positions;
    consensus.summary.reference_frame_runs += summary.reference_frame_runs;
    consensus.summary.phase_loss_runs += summary.phase_loss_runs;
}

/// The sample consensus plugin.
pub const PLUGIN: PluginDescriptor = PluginDescriptor {
    id: "consensus",
    family: PluginFamily::Core,
    version: 22,
    provides: &[Contract::Consensus],
    requires: &[Contract::CalledReads],
    config_sections: &["consensus"],
};

#[cfg(test)]
mod tests {
    use dna_align::{AlignmentConfig, RawAlignmentConfig, align_best};
    use dna_call::RawVariantCallingConfig;
    use dna_call::model::variant::VariantCallingResult;
    use dna_call::read_call::RawCallConfig;
    use dna_kernel::error::Result;
    use dna_kernel::model::reference::ReferenceTopology;
    use dna_kernel::read_evidence::{EvidenceReason, Placement, ReadEvidence};

    use super::*;

    const REFERENCE: &str = "TTGACGTCAGTACGATCGTACCTGAGTACGAGGACT";

    fn reference() -> Reference {
        Reference {
            name: "ref".into(),
            sequence: REFERENCE.into(),
            topology: ReferenceTopology::Linear,
            sequence_sha256: String::new(),
        }
    }

    fn alignment() -> Result<AlignmentConfig> {
        RawAlignmentConfig {
            match_score: 3,
            mismatch_score: -5,
            ambiguous_score: 0,
            gap_open_score: -10,
            gap_extension_score: -4,
            minimum_callable_bases: 1,
            minimum_identity: 0.8,
            max_stretch: 64,
        }
        .validate()
    }

    fn config() -> Result<CallConfig> {
        RawCallConfig {
            alignment: RawAlignmentConfig {
                match_score: 3,
                mismatch_score: -5,
                ambiguous_score: 0,
                gap_open_score: -10,
                gap_extension_score: -4,
                minimum_callable_bases: 1,
                minimum_identity: 0.8,
                max_stretch: 64,
            },
            variant_calling: RawVariantCallingConfig {
                max_indel_length: 50,
                read_end_margin: 0,
            },
        }
        .validate()
    }

    /// A read of `sequence` placed on the reference, with the calls in
    /// `masked` masked as dephased.
    fn placed(name: &str, sequence: &str, masked: &[usize]) -> Result<CalledRead> {
        placed_with(name, sequence, masked, &[])
    }

    /// As [`placed`], with the calls in `unresolved` masked as unresolved.
    fn placed_with(
        name: &str,
        sequence: &str,
        masked: &[usize],
        unresolved: &[usize],
    ) -> Result<CalledRead> {
        let mut evidence = ReadEvidence::clean(sequence);
        let masks = masked
            .iter()
            .map(|&index| {
                (
                    Placement::Anchoring(EvidenceReason::new("dephased_signal")),
                    index,
                )
            })
            .chain(unresolved.iter().map(|&index| {
                (
                    Placement::Unresolved(EvidenceReason::new("unresolved_signal")),
                    index,
                )
            }));
        for (placement, index) in masks {
            let call = evidence.calls()[index].masked(placement);
            evidence = evidence.with_call(index, call);
        }
        let alignment = align_best(&evidence, &reference(), &alignment()?)?;
        Ok(CalledRead {
            input_name: name.into(),
            input_sha256: name.into(),
            reference_sha256: String::new(),
            configuration_sha256: String::new(),
            evidence,
            alignment,
            variants: VariantCallingResult::fixture(Vec::new(), Vec::new()),
        })
    }

    /// The observations of a read `read_sequence` on `reference_sequence`,
    /// both `ACGTACGTAC G CCCCC T AG` + `tail` in shape, with calls 19-25
    /// masked as a segment that starts right after a repeat run, under a
    /// read-end margin of 8.
    fn exempt_observations(
        reference_sequence: &str,
        read_sequence: &str,
    ) -> std::result::Result<ReadObservations, Box<dyn std::error::Error>> {
        let reference = Reference {
            name: "ref".into(),
            sequence: reference_sequence.into(),
            topology: ReferenceTopology::Linear,
            sequence_sha256: String::new(),
        };
        let mut evidence = ReadEvidence::clean(read_sequence);
        for index in 19..26 {
            let mut call = evidence.calls()[index];
            call = call.masked(Placement::Unresolved(EvidenceReason::new(
                "post_homopolymer",
            )));
            evidence = evidence.with_call(index, call);
        }
        // The modality marks the C run (11-15) and its bounds as a run bounded
        // before the post-repeat phase loss at 19 (ADR-0076).
        for index in 10..=16 {
            let mut call = evidence.calls()[index];
            call.post_repeat.call = (11..=15).contains(&index);
            call.post_repeat.bound = index == 10 || index == 16;
            evidence = evidence.with_call(index, call);
        }
        let evidence = evidence.with_signal_edges(8);
        let alignment = align_best(&evidence, &reference, &alignment()?)?;
        let read = CalledRead {
            input_name: "r".into(),
            input_sha256: "r".into(),
            reference_sha256: String::new(),
            configuration_sha256: String::new(),
            evidence,
            alignment,
            variants: VariantCallingResult::fixture(Vec::new(), Vec::new()),
        };
        let settings = RawVariantCallingConfig {
            max_indel_length: 50,
            read_end_margin: 8,
        }
        .validate()?;
        Ok(observe(&read, &settings, reference.sequence.as_bytes(), 4))
    }

    // SRS-CONS-003
    #[test]
    fn an_exempt_run_trusts_its_matches_but_never_a_substituted_call()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        // The read's C run 11-15 covers the reference A at 11 and is bounded
        // by G10 and T16, which match the reference.
        let tail = "AAAAAACGTACGTACGTACG";
        let observed = exempt_observations(
            &format!("ACGTACGTACGACCCCTAG{tail}"),
            &format!("ACGTACGTACGCCCCCTAG{tail}"),
        )?;
        let clean = |position: usize| observed.positions.get(&position).map(|site| site.clean);
        // The substituted call keeps `read_end`; the run's matches are clean;
        // the bound that reads the reference base is clean too, since it can
        // never make a variant.
        assert_eq!(clean(11), Some(false));
        assert!((12..=15).all(|position| clean(position) == Some(true)));
        assert_eq!(clean(16), Some(true));
        Ok(())
    }

    // SRS-CONS-003
    #[test]
    fn an_exempt_run_bound_on_another_base_keeps_its_trust()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        // The read's bound at 16 reads A where the reference has T.
        let tail = "AAAAAACGTACGTACGTACG";
        let observed = exempt_observations(
            &format!("ACGTACGTACGCCCCCTAG{tail}"),
            &format!("ACGTACGTACGCCCCCAAG{tail}"),
        )?;
        let clean = |position: usize| observed.positions.get(&position).map(|site| site.clean);
        assert!((11..=15).all(|position| clean(position) == Some(true)));
        assert_eq!(clean(16), Some(false));
        Ok(())
    }

    // SRS-CONS-003
    #[test]
    fn each_placement_and_verdict_maps_to_one_trust_level()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        use super::observe::{Trust, observe};
        // A read equal to the reference: call `i` observes position `i`.
        let mut read = placed_with("r", REFERENCE, &[5], &[8])?;
        let peak = EvidenceReason::new("peak_below_minimum");
        let mut vetoed = read.evidence.calls()[11];
        vetoed.identity.withheld.insert(peak)?;
        let mut base_only = read.evidence.calls()[14];
        base_only.bound = dna_kernel::read_evidence::Bound::Withheld;
        read.evidence = read.evidence.with_call(11, vetoed).with_call(14, base_only);
        let observed = observe(&read, &config()?.variant_calling, REFERENCE.as_bytes(), 4);
        let char_at = |position: usize| {
            observed
                .positions
                .get(&position)
                .and_then(|observation| observation.chars.first().copied())
        };
        let level = |position: usize| char_at(position).map(|char| char.trust);
        // Scored, trusted, nothing withheld: it shows its base.
        assert_eq!(level(2), Some(Trust::Clean));
        // Scored with a reason against its base: it counts towards its run.
        assert_eq!(level(11), Some(Trust::Unmasked));
        // Masked but anchoring: it shows a run's end.
        assert_eq!(level(5), Some(Trust::Anchoring));
        // Masked as unresolved: nothing.
        assert_eq!(level(8), Some(Trust::Nothing));
        // Vouched for by its base only: it shows its base, but neither counts
        // towards a run nor shows a run's end.
        let only = char_at(14).ok_or("no call at 14")?;
        assert!(only.clean() && !only.unmasked() && !only.resolving());
        Ok(())
    }

    // SRS-CONS-003
    #[test]
    fn each_use_of_a_call_reads_its_own_reasons()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        use super::observe::observe;
        // `REFERENCE` with an A deleted after call 10 and a G inserted after
        // call 20 (0-based reference positions 11 and after 20).
        let mut sequence = REFERENCE.to_owned();
        sequence.insert(21, 'G');
        sequence.remove(11);
        let mut read = placed("r", &sequence, &[])?;
        let mixed = EvidenceReason::new("mixed_supporting_dna");
        let peak = EvidenceReason::new("peak_below_minimum");
        // Mixed signal on the inserted G withholds its base, not the
        // insertion; a weak peak beside the deletion withholds that call's
        // base, not the deletion.
        let mut inserted = read.evidence.calls()[20];
        inserted.identity.withheld.insert(mixed)?;
        let mut flank = read.evidence.calls()[10];
        flank.identity.withheld.insert(peak)?;
        flank.indel.withheld.insert(peak)?;
        read.evidence = read.evidence.with_call(20, inserted).with_call(10, flank);
        let observed = observe(&read, &config()?.variant_calling, REFERENCE.as_bytes(), 4);
        assert_eq!(
            observed.positions.get(&11).map(|deleted| deleted.clean),
            Some(true)
        );
        assert_eq!(
            observed.positions.get(&10).map(|base| base.clean),
            Some(false)
        );
        let junction = observed
            .junctions
            .values()
            .find(|junction| !junction.chars.is_empty());
        assert_eq!(junction.map(|junction| junction.clean), Some(true));
        Ok(())
    }

    // SRS-CONS-009
    #[test]
    fn an_indel_with_an_undecided_neighbour_stays_undecided() {
        let decision = |value: Option<&str>| Decision {
            value: value.map(Into::into),
            state: if value.is_some() {
                SiteState::Called
            } else {
                SiteState::Unresolved
            },
            supporting: Vec::new(),
            opposing: Vec::new(),
            uninformative: Vec::new(),
            runs: Vec::new(),
        };
        // 10 undecided, 11-12 deleted, 13 an A; then 20 a C, 21-22 deleted,
        // 23 a G.
        let mut decided = BTreeMap::from([
            (10, (10, decision(None))),
            (11, (12, decision(Some("")))),
            (13, (13, decision(Some("A")))),
            (20, (20, decision(Some("C")))),
            (21, (22, decision(Some("")))),
            (23, (23, decision(Some("G")))),
        ]);
        // 29 undecided, 30-31 a CC run lengthened to CCC that a read shows
        // in phase: its reads show the run's bounds themselves.
        let mut shown = decision(Some("CCC"));
        shown.runs = vec![run('C', 3, RunEvidence::InPhase)];
        decided.insert(29, (29, decision(None)));
        decided.insert(30, (31, shown));
        unplace_unbounded_indels(&mut decided, 100);
        assert_eq!(decided[&11].1.state, SiteState::Unresolved);
        assert!(decided[&11].1.value.is_none());
        assert_eq!(decided[&21].1.value.as_deref(), Some(""));
        assert_eq!(decided[&30].1.value.as_deref(), Some("CCC"));
    }

    /// A tract decode of a run of `run` entered at call `call_start`, the
    /// anchor `mark` after it: its first run `first` or one longer, and the
    /// total `total` or one longer, each at shares 0.95 and 0.05.
    fn decoded(
        run: char,
        mark: char,
        call_start: usize,
        first: usize,
        total: usize,
    ) -> dna_kernel::read_evidence::TractLengths {
        use dna_kernel::read_evidence::{LengthMarginal, Margin, StructureFit, TractLengths};
        let marginal = |length: usize| LengthMarginal {
            lengths: vec![length, length + 1],
            shares: vec![0.95, 0.05],
            errors: vec![0.01, 0.01],
            dominant: Some(length),
            next: Some(length + 1),
            z: Some(10.0),
            z_adjusted: Some(10.0),
            margin: Margin::Clear,
        };
        TractLengths {
            onset: call_start + total,
            call_start,
            extent: total,
            window: total,
            run,
            mark: Some(mark),
            candidates: vec![
                StructureFit {
                    mark: None,
                    residual: 0.02,
                },
                StructureFit {
                    mark: Some(mark),
                    residual: 0.01,
                },
            ],
            residual: 0.01,
            sigma: 0.01,
            autocorrelation: 0.0,
            first: marginal(first),
            total: marginal(total),
            cleaned: None,
        }
    }

    // SRS-CONS-013
    #[test]
    fn both_strands_decide_a_window_run_left_as_an_estimate()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        // A stretch of C7 T C5 at 11-23 (1-based), the anchor T at 18, its
        // window by its own structure. Both strands decode C8 T C5: the
        // forward read enters run 1 at 11, the reverse read run 2 at 23, and
        // each derives the other run.
        let window_reference = Reference {
            name: "ref".into(),
            sequence: "TTGACGTCAGCCCCCCCTCCCCCGATCGTACCTGAGTACGAGGACT".into(),
            topology: ReferenceTopology::Linear,
            sequence_sha256: String::new(),
        };
        let read = |name: &str, sequence: &str, tract| -> Result<CalledRead> {
            let evidence = ReadEvidence::clean(sequence).with_tracts(vec![tract])?;
            let alignment = align_best(&evidence, &window_reference, &alignment()?)?;
            Ok(CalledRead {
                input_name: name.into(),
                input_sha256: name.into(),
                reference_sha256: String::new(),
                configuration_sha256: String::new(),
                evidence,
                alignment,
                variants: VariantCallingResult::fixture(Vec::new(), Vec::new()),
            })
        };
        let sequence = window_reference.sequence.as_str();
        let reverse_sequence = dna_kernel::model::nucleotide::reverse_complement(sequence);
        let forward = read("f", sequence, decoded('C', 'T', 10, 8, 14))?;
        let reverse = read(
            "r",
            &reverse_sequence,
            decoded('G', 'A', sequence.len() - 23, 5, 14),
        )?;
        let decided_with = |evidence: RunEvidence| {
            let decision = Decision {
                value: Some("CCCCCCCTCCCCC".into()),
                state: SiteState::Called,
                supporting: vec!["f".into(), "r".into()],
                opposing: Vec::new(),
                uninformative: Vec::new(),
                runs: vec![
                    run('C', 7, evidence),
                    run('T', 1, RunEvidence::InPhase),
                    run('C', 5, RunEvidence::InPhase),
                ],
            };
            BTreeMap::from([(10, (22, decision))])
        };
        let settle =
            |evidence: RunEvidence| -> std::result::Result<Decision, Box<dyn std::error::Error>> {
                let mut decided = decided_with(evidence);
                two_strand(
                    &mut decided,
                    &[&forward, &reverse],
                    &window_reference,
                    &config()?,
                    &config::tests::config(),
                );
                Ok(decided.remove(&10).ok_or("missing stretch")?.1)
            };
        // A phase-lost run 1 takes the strands' length; run 2, in phase,
        // keeps its own.
        let decision = settle(RunEvidence::PhaseLoss)?;
        assert_eq!(
            decision.runs,
            [
                run('C', 8, RunEvidence::TwoStrand),
                run('T', 1, RunEvidence::InPhase),
                run('C', 5, RunEvidence::InPhase)
            ]
        );
        assert_eq!(decision.value.as_deref(), Some("CCCCCCCCTCCCCC"));
        // A run whose ends two reads confirm keeps its length.
        let decision = settle(RunEvidence::ConfirmedEnd)?;
        assert_eq!(decision.runs[0], run('C', 7, RunEvidence::ConfirmedEnd));
        assert_eq!(decision.value.as_deref(), Some("CCCCCCCTCCCCC"));
        Ok(())
    }

    fn site_at(consensus: &Consensus, position: usize) -> Option<&ConsensusSite> {
        consensus
            .sites
            .iter()
            .find(|site| (site.first_0based..=site.last_0based).contains(&position))
    }

    // SRS-CONS-004, SRS-CONS-005
    #[test]
    fn agreeing_reads_give_their_haplotype() -> std::result::Result<(), Box<dyn std::error::Error>>
    {
        // A substitution at 10 (T>G) and a GG insertion after 20 in both reads.
        let haplotype = "TTGACGTCAGGACGATCGTACGGCTGAGTACGAGGACT";
        let forward = placed("f", haplotype, &[])?;
        let reverse = placed(
            "r",
            &dna_kernel::model::nucleotide::reverse_complement(haplotype),
            &[],
        )?;
        let consensus = build(
            &[&forward, &reverse],
            &reference(),
            &config()?,
            &config::tests::config(),
        );
        assert_eq!(consensus.segments.len(), 1);
        assert_eq!(consensus.segments[0].sequence, haplotype);
        // The substituted G joins the G before it: the stretch is that run.
        let substitution = site_at(&consensus, 10).ok_or("missing substitution site")?;
        assert_eq!(
            (substitution.first_0based, substitution.last_0based),
            (9, 10)
        );
        assert_eq!(substitution.call.as_deref(), Some("GG"));
        assert_eq!(substitution.state, SiteState::Called);
        assert_eq!(substitution.supporting, ["f", "r"]);
        assert_eq!(
            substitution.runs,
            [ConsensusRun {
                base: 'G',
                length: 2,
                evidence: RunEvidence::InPhase
            }]
        );
        // The insertion is decided over the whole CC run it lengthens.
        let insertion = site_at(&consensus, 20).ok_or("missing insertion site")?;
        assert_eq!((insertion.first_0based, insertion.last_0based), (20, 21));
        assert_eq!(insertion.call.as_deref(), Some("CGGC"));
        Ok(())
    }

    // SRS-CONS-005
    #[test]
    fn a_lone_read_against_a_reference_read_is_contested()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        let variant = placed("v", "TTGACGTCAGGACGATCGTACCTGAGTACGAGGACT", &[])?;
        let plain = placed("p", REFERENCE, &[])?;
        let consensus = build(
            &[&variant, &plain],
            &reference(),
            &config()?,
            &config::tests::config(),
        );
        assert_eq!(consensus.segments[0].sequence, REFERENCE);
        let site = site_at(&consensus, 10).ok_or("missing contested site")?;
        assert_eq!(site.state, SiteState::Contested);
        assert_eq!(
            (site.supporting.as_slice(), site.opposing.as_slice()),
            (["p".to_owned()].as_slice(), ["v".to_owned()].as_slice())
        );
        assert_eq!(consensus.summary.contested_sites, 1);
        Ok(())
    }

    /// `sequence` placed with its call `index` withheld as mixed and, when
    /// `reference_only`, asserted reference-only, its second base `second`
    /// standing out.
    fn confirming(
        name: &str,
        sequence: &str,
        index: usize,
        reference_only: bool,
    ) -> std::result::Result<CalledRead, Box<dyn std::error::Error>> {
        let mut read = placed(name, sequence, &[])?;
        let mut call = read.evidence.calls()[index];
        call.identity
            .withheld
            .insert(EvidenceReason::new("mixed_supporting_dna"))?;
        call.reference_only = reference_only;
        call.second = Some(dna_kernel::read_evidence::SecondBase {
            base: 'C',
            share: 0.3,
            stands_out: true,
        });
        read.evidence = read.evidence.with_call(index, call);
        Ok(read)
    }

    // SRS-CONS-010
    #[test]
    fn a_reference_only_call_confirms_the_reference_and_nothing_else()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        let build_alone =
            |read: &CalledRead| -> std::result::Result<Consensus, Box<dyn std::error::Error>> {
                Ok(build(
                    &[read],
                    &reference(),
                    &config()?,
                    &config::tests::config(),
                ))
            };
        // Withheld, the call at 14 leaves its position undecided.
        let withheld = build_alone(&confirming("r", REFERENCE, 14, false)?)?;
        assert_eq!(&withheld.segments[0].sequence[13..16], "GNT");
        // Reference-only, it confirms the reference there, and the position
        // is marked as decided by that confirmation alone.
        let confirmed = build_alone(&confirming("r", REFERENCE, 14, true)?)?;
        assert_eq!(confirmed.segments[0].sequence, REFERENCE);
        assert_eq!(confirmed.reference_only, BTreeSet::from([14]));
        // A reference-only call reading another base than the reference's
        // decides nothing.
        let other = "TTGACGTCAGTACGTTCGTACCTGAGTACGAGGACT";
        let differing = build_alone(&confirming("r", other, 14, true)?)?;
        assert_eq!(&differing.segments[0].sequence[13..16], "GNT");
        assert!(differing.reference_only.is_empty());
        Ok(())
    }

    // SRS-CONS-010
    #[test]
    fn a_confirmation_beside_a_differing_site_is_withdrawn()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        // The read substitutes G>T at 13, cleanly, and its call at 14 is a
        // reference-only A: the substitution's neighbour stays undecided, as
        // the read's verdict leaves it.
        let substituted = "TTGACGTCAGTACTATCGTACCTGAGTACGAGGACT";
        let read = confirming("r", substituted, 14, true)?;
        let consensus = build(&[&read], &reference(), &config()?, &config::tests::config());
        assert_eq!(&consensus.segments[0].sequence[12..16], "CTNT");
        assert!(consensus.reference_only.is_empty());
        Ok(())
    }

    // SRS-CONS-002
    #[test]
    fn reads_of_one_file_name_count_by_content()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        // Three reads of one strand and one file name: two substitute G>T at
        // 13, one reads the reference. Each counts once, by its content.
        let substituted = "TTGACGTCAGTACTATCGTACCTGAGTACGAGGACT";
        let mut reads = vec![
            placed("sha-1", substituted, &[])?,
            placed("sha-2", substituted, &[])?,
            placed("sha-3", REFERENCE, &[])?,
        ];
        for read in &mut reads {
            read.input_name = "read.ab1".into();
        }
        let consensus = build(
            &reads.iter().collect::<Vec<_>>(),
            &reference(),
            &config()?,
            &config::tests::config(),
        );
        let site = consensus
            .sites
            .iter()
            .find(|site| site.first_0based <= 13 && 13 <= site.last_0based)
            .ok_or("the substitution must be a site")?;
        assert_eq!(
            site.call.as_deref().map(|call| call.contains('T')),
            Some(true)
        );
        assert_eq!(site.supporting, ["sha-1", "sha-2"]);
        assert_eq!(site.opposing, ["sha-3"]);
        Ok(())
    }

    /// `sequence` placed with the calls at `edges` at a signal edge and
    /// reference-only: confirmations the modality's edge does not stop.
    fn at_edges(
        name: &str,
        sequence: &str,
        edges: &[usize],
    ) -> std::result::Result<CalledRead, Box<dyn std::error::Error>> {
        let mut read = placed(name, sequence, &[])?;
        for &index in edges {
            let mut call = read.evidence.calls()[index];
            call.identity.edge = true;
            call.indel.edge = true;
            call.reference_only = true;
            read.evidence = read.evidence.with_call(index, call);
        }
        Ok(read)
    }

    // SRS-CONS-010
    #[test]
    fn a_reference_only_call_at_an_edge_confirms_clear_of_the_read_ends()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        // At an edge, the call at 14 would leave its position undecided; as
        // a reference-only call it confirms the reference there.
        let read = at_edges("r", REFERENCE, &[14])?;
        let consensus = build(&[&read], &reference(), &config()?, &config::tests::config());
        assert!(consensus.reference_only.contains(&14));
        assert_eq!(consensus.segments[0].sequence, REFERENCE);
        Ok(())
    }

    // SRS-CONS-010
    #[test]
    fn a_vouched_call_within_the_end_margin_confirms_the_reference()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        use super::observe::observe;
        // A read equal to the reference, its ends unvouched: calls 0-2 lie
        // within a margin of three calls of its start, and call `i` observes
        // position `i`.
        let mut read = placed("r", REFERENCE, &[])?;
        let mut edge = read.evidence.calls()[0];
        edge.identity.edge = true;
        let mut second = read.evidence.calls()[1];
        second.second = Some(dna_kernel::read_evidence::SecondBase {
            base: 'A',
            share: 0.3,
            stands_out: true,
        });
        read.evidence = read.evidence.with_call(0, edge).with_call(1, second);
        let settings = RawVariantCallingConfig {
            max_indel_length: 50,
            read_end_margin: 3,
        }
        .validate()?;
        let observed = observe(&read, &settings, REFERENCE.as_bytes(), 4);
        let site = |position: usize| {
            observed
                .positions
                .get(&position)
                .map(|site| (site.clean, site.confirms))
        };
        // The modality's edge verdict, or a second base standing out, leaves
        // the reference unconfirmed; the call at 2, vouched for its base,
        // confirms it. Beyond the margin, calls are clean on their own.
        assert_eq!(site(0), Some((false, false)));
        assert_eq!(site(1), Some((false, false)));
        assert_eq!(site(2), Some((true, true)));
        assert_eq!(site(3), Some((true, false)));
        Ok(())
    }

    // SRS-CONS-010
    #[test]
    fn a_confirmation_beside_any_observed_difference_is_withdrawn()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        // A second read reads T at 16 (C>T), uncleanly: a shifted read can
        // show that, so the confirmation at 15 beside it is withdrawn, while
        // the one at 14 stands.
        let confirming_read = at_edges("r", REFERENCE, &[14, 15])?;
        let mut shifted = placed("s", "TTGACGTCAGTACGATTGTACCTGAGTACGAGGACT", &[])?;
        for index in [14, 15, 16] {
            let mut call = shifted.evidence.calls()[index];
            call.identity.edge = true;
            shifted.evidence = shifted.evidence.with_call(index, call);
        }
        let consensus = build(
            &[&confirming_read, &shifted],
            &reference(),
            &config()?,
            &config::tests::config(),
        );
        assert!(consensus.reference_only.contains(&14));
        assert!(!consensus.reference_only.contains(&15));
        Ok(())
    }

    // SRS-CONS-010
    #[test]
    fn a_stretch_is_decided_without_confirmations() {
        // T12C merges the tract into one stretch (7-16). A second read reads
        // the tract uncleanly but confirms the reference at 9: the stretch's
        // decision is the same with and without that confirmation. This
        // guards the outcome only: the stretch here does not fall back to its
        // position-by-position decision, where the exclusion acts.
        let merged = "GTCAAAACCCCCCCCCCATG";
        let observed = |confirm: bool| {
            let clean = rows("f", 0, TRACT, merged, &"c".repeat(TRACT.len()));
            let mut other = rows("g", 0, TRACT, TRACT, &"v".repeat(TRACT.len()));
            if confirm && let Some(observation) = other.positions.get_mut(&9) {
                observation.confirms = true;
                observation.clean = true;
            }
            vec![clean, other]
        };
        let decided = |confirm: bool| {
            decide_positions(&observed(confirm), TRACT.as_bytes(), 4)
                .get(&7)
                .map(|(last, decision)| (*last, decision.value.clone(), decision.state))
        };
        assert!(decided(false).is_some());
        assert_eq!(decided(true), decided(false));
    }

    // SRS-CONS-010
    #[test]
    fn a_confirmation_never_closes_an_unobserved_junction_in_a_long_run() {
        // The C run at 7-11 of the tract, five long. A confirmation at 9 with
        // no read observing the junction after 8 cleanly is withdrawn where
        // such a run is long enough to be a stretch, and stands otherwise.
        let confirmations = |seeded_run: usize| {
            let mut read = rows("r", 0, TRACT, TRACT, &"c".repeat(TRACT.len()));
            if let Some(observation) = read.positions.get_mut(&9) {
                observation.confirms = true;
            }
            if let Some(junction) = read.junctions.get_mut(&8) {
                junction.clean = false;
            }
            let mut observed = vec![read];
            withdraw_unobserved_run_junctions(&mut observed, TRACT.as_bytes(), seeded_run);
            observed[0].positions[&9].confirms
        };
        assert!(!confirmations(4));
        assert!(confirmations(6));
    }

    // SRS-CONS-010
    #[test]
    fn a_second_sequence_any_read_shows_withdraws_a_confirmation() {
        // Two reads at position 14: the first confirms the reference with a
        // second C; the second does not observe the position cleanly.
        let observations = |first_forward: bool, second: (bool, Option<char>, bool)| {
            let mut first = rows("a", 0, REFERENCE, REFERENCE, &"c".repeat(REFERENCE.len()));
            first.forward = first_forward;
            if let Some(observation) = first.positions.get_mut(&14) {
                observation.confirms = true;
            }
            first.seconds.insert(14, 'C');
            let mut other = rows("b", 0, REFERENCE, REFERENCE, &"c".repeat(REFERENCE.len()));
            other.forward = second.0;
            if let Some(base) = second.1 {
                other.seconds.insert(14, base);
            }
            if second.2 {
                other.unexplained.insert(14);
            }
            let mut observed = vec![first, other];
            hide_confirmations(&mut observed);
            observed[0].positions[&14].confirms
        };
        // Nothing on the other read: the confirmation stands.
        assert!(observations(true, (false, None, false)));
        // The other strand shows the same second base standing out.
        assert!(!observations(true, (false, Some('C'), false)));
        // The same strand doing so is not a second sequence by itself.
        assert!(observations(true, (true, Some('C'), false)));
        // A read of the confirming strand with a second peak its modality
        // does not explain.
        assert!(!observations(true, (true, None, true)));
        assert!(!observations(true, (true, Some('A'), true)));
    }

    // SRS-CONS-010
    #[test]
    fn an_unexplained_second_on_the_other_strand_alone_keeps_a_confirmation() {
        // Read `a` (forward) confirms the reference at 14 with a second C,
        // read `b` (reverse) shows an unexplained second peak there, and read
        // `c` (forward) optionally shows a second base standing out.
        let observations = |second: Option<char>, third: Option<char>| {
            let read = |name: &str, forward: bool| {
                let mut read = rows(name, 0, REFERENCE, REFERENCE, &"c".repeat(REFERENCE.len()));
                read.forward = forward;
                read
            };
            let mut first = read("a", true);
            if let Some(observation) = first.positions.get_mut(&14) {
                observation.confirms = true;
            }
            first.seconds.insert(14, 'C');
            let mut other = read("b", false);
            other.unexplained.insert(14);
            if let Some(base) = second {
                other.seconds.insert(14, base);
            }
            let mut same = read("c", true);
            if let Some(base) = third {
                same.seconds.insert(14, base);
            }
            let mut observed = vec![first, other, same];
            hide_confirmations(&mut observed);
            observed[0].positions[&14].confirms
        };
        // A second base the confirming strand does not show: one strand
        // alone, so the confirmation stands.
        assert!(observations(Some('A'), None));
        assert!(observations(Some('A'), Some('G')));
        // The confirming strand shows the same second base standing out.
        assert!(!observations(Some('A'), Some('A')));
        // The confirming call's own second base: the other-strand guard.
        assert!(!observations(Some('C'), None));
        // No second base standing out to compare: still withdrawn.
        assert!(!observations(None, None));
    }

    // SRS-CONS-012
    #[test]
    fn positions_only_unreliable_spans_observe_are_unreliable_only()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        let reason = EvidenceReason::new("poor_start_signal");
        let mut poor = placed("p", REFERENCE, &[])?;
        poor.evidence =
            poor.evidence
                .with_unreliable(vec![dna_kernel::read_evidence::UnreliableSpan {
                    calls: 0..5,
                    reason,
                }])?;
        // Alone, the read observes its span's positions only unreliably.
        let alone = build(&[&poor], &reference(), &config()?, &config::tests::config());
        assert_eq!(
            alone.unreliable_only.into_iter().collect::<Vec<_>>(),
            (0..5)
                .map(|position| (position, reason))
                .collect::<Vec<_>>()
        );
        // A read that observes 3 and 4 cleanly leaves 0..=2.
        let late = placed("l", &REFERENCE[3..], &[])?;
        let both = build(
            &[&poor, &late],
            &reference(),
            &config()?,
            &config::tests::config(),
        );
        assert_eq!(
            both.unreliable_only.keys().copied().collect::<Vec<_>>(),
            [0, 1, 2]
        );
        Ok(())
    }

    // SRS-CONS-003, SRS-CONS-006
    #[test]
    fn masked_calls_never_decide() -> std::result::Result<(), Box<dyn std::error::Error>> {
        // The variant read's call at 10 is masked: the reference read decides.
        let variant = placed("v", "TTGACGTCAGGACGATCGTACCTGAGTACGAGGACT", &[10])?;
        let plain = placed("p", REFERENCE, &[])?;
        let consensus = build(
            &[&variant, &plain],
            &reference(),
            &config()?,
            &config::tests::config(),
        );
        assert_eq!(consensus.segments[0].sequence, REFERENCE);
        assert!(site_at(&consensus, 10).is_none());
        // Alone, a masked call leaves the position unresolved.
        let alone = build(
            &[&variant],
            &reference(),
            &config()?,
            &config::tests::config(),
        );
        assert_eq!(&alone.segments[0].sequence[8..12], "AGNA");
        assert_eq!(alone.summary.unresolved_positions, 1);
        Ok(())
    }

    // SRS-CONS-003
    #[test]
    fn no_masked_call_places_a_substitution_or_an_indel()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        // A substitution (A14T), an insertion (a G after 14), and a deletion
        // (of 14), each read only through masked calls, anchoring or
        // unresolved: a masked call shows at most a run's end (ADR-0071), so
        // none is decided, beside a reference read or alone.
        let substituted = "TTGACGTCAGTACGTTCGTACCTGAGTACGAGGACT";
        let inserted = "TTGACGTCAGTACGAGTCGTACCTGAGTACGAGGACT";
        let deleted = "TTGACGTCAGTACGTCGTACCTGAGTACGAGGACT";
        let cases: [(&str, &[usize]); 3] = [
            (substituted, &[14]),
            (inserted, &[15]),
            (deleted, &[13, 14]),
        ];
        let plain = placed("p", REFERENCE, &[])?;
        for (sequence, masked) in cases {
            for read in [
                placed("v", sequence, masked)?,
                placed_with("v", sequence, &[], masked)?,
            ] {
                let both = build(
                    &[&read, &plain],
                    &reference(),
                    &config()?,
                    &config::tests::config(),
                );
                assert_eq!(both.segments[0].sequence, REFERENCE);
                let alone = build(&[&read], &reference(), &config()?, &config::tests::config());
                // Every segment reads its reference interval, or `N`.
                for segment in &alone.segments {
                    let expected = &REFERENCE[segment.start_0based..=segment.end_0based];
                    assert_eq!(segment.sequence.len(), expected.len(), "{sequence}");
                    assert!(
                        segment
                            .sequence
                            .chars()
                            .zip(expected.chars())
                            .all(|(call, reference)| call == reference || call == 'N'),
                        "{sequence} gave {}",
                        segment.sequence
                    );
                }
            }
        }
        Ok(())
    }

    // SRS-CONS-003
    #[test]
    fn an_insertion_into_a_run_the_read_does_not_bound_never_decides()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        // A C inserted into the CC run at 20-21; the A before the run is call 19.
        let lengthened = "TTGACGTCAGTACGATCGTACCCTGAGTACGAGGACT";
        let bounded = build(
            &[&placed("b", lengthened, &[])?],
            &reference(),
            &config()?,
            &config::tests::config(),
        );
        assert_eq!(bounded.segments[0].sequence, lengthened);
        // Unresolved, the A no longer shows where the run starts, so the read
        // cannot tell the run's length and the reference length stays.
        let unbounded = build(
            &[&placed_with("u", lengthened, &[], &[19])?],
            &reference(),
            &config()?,
            &config::tests::config(),
        );
        let mut expected = REFERENCE.to_owned();
        expected.replace_range(19..20, "N");
        assert_eq!(unbounded.segments[0].sequence, expected);
        assert!(site_at(&unbounded, 20).is_none());
        Ok(())
    }
    /// A read observed over `reference` from 0-based `first`, given as two
    /// alignment rows (`-` for an inserted or deleted base) and one kind per
    /// column: `c` clean, `m` masked but anchoring, `v` unmasked but vetoed;
    /// an `N` in the read row is an unresolved call.
    fn rows(
        name: &str,
        first: usize,
        reference: &str,
        read: &str,
        kinds: &str,
    ) -> ReadObservations {
        use super::observe::Observation;
        let char_of = |base: char, kind: char, position: Option<usize>| {
            let canonical = base != 'N';
            Char {
                base,
                trust: match kind {
                    _ if !canonical => Trust::Nothing,
                    'c' => Trust::Clean,
                    'v' => Trust::Unmasked,
                    _ => Trust::Anchoring,
                },
                position,
                beside_unresolved: false,
                base_only: false,
            }
        };
        let mut positions = BTreeMap::new();
        let mut junctions = BTreeMap::new();
        let mut position = first;
        let mut previous: Option<usize> = None;
        let mut inserted: Vec<Char> = Vec::new();
        for ((reference, read), kind) in reference.chars().zip(read.chars()).zip(kinds.chars()) {
            if reference == '-' {
                inserted.push(char_of(read, kind, None));
                continue;
            }
            if let Some(previous) = previous {
                junctions.insert(
                    previous,
                    Observation {
                        value: SiteValue::Insertion(
                            inserted.iter().map(|char| char.base).collect(),
                        ),
                        confirms: false,
                        clean: inserted.iter().all(|char| char.clean()),
                        chars: std::mem::take(&mut inserted),
                    },
                );
            }
            let observation = if read == '-' {
                Observation {
                    value: SiteValue::Deletion,
                    confirms: false,
                    clean: true,
                    chars: Vec::new(),
                }
            } else {
                let char = char_of(read, kind, Some(position));
                Observation {
                    value: SiteValue::Base(read),
                    confirms: false,
                    clean: char.clean(),
                    chars: vec![char],
                }
            };
            positions.insert(position, observation);
            previous = Some(position);
            position += 1;
        }
        ReadObservations {
            forward: true,
            seconds: BTreeMap::new(),
            unexplained: BTreeSet::new(),
            edges: BTreeMap::new(),
            id: name.into(),
            positions,
            junctions,
        }
    }

    /// A forward read of the reference and a reverse read that differs at 4,
    /// both `edge` calls from their ends everywhere.
    fn dissent(reverse: &str, kinds: &str, edge: usize) -> Vec<ReadObservations> {
        const SHORT: &str = "ACGTACGTAC";
        let mut forward = rows("f", 0, SHORT, SHORT, "cccccccccc");
        let mut opposing = rows("r", 0, SHORT, reverse, kinds);
        opposing.forward = false;
        for read in [&mut forward, &mut opposing] {
            read.edges = (0..SHORT.len()).map(|position| (position, edge)).collect();
        }
        vec![forward, opposing]
    }

    /// The decision of the site or stretch holding `position`.
    fn decided_at(observed: &[ReadObservations], position: usize) -> Decision {
        let decided = decide_positions(observed, b"ACGTACGTAC", 4);
        decided
            .range(..=position)
            .next_back()
            .map_or_else(|| decide(&[], "A"), |(_, (_, decision))| decision.clone())
    }

    /// The decision of the stretch holding `position` once one-strand
    /// insertions are contested (SRS-CONS-015).
    fn contested_at(observed: &[ReadObservations], bases: &[u8], position: usize) -> Decision {
        let mut decided = decide_positions(observed, bases, 4);
        insertion::contest_one_strand(&mut decided, observed, bases, 4);
        decided
            .range(..=position)
            .next_back()
            .map_or_else(|| decide(&[], "A"), |(_, (_, decision))| decision.clone())
    }

    // SRS-CONS-015
    #[test]
    fn a_short_insertion_the_other_strand_reads_through_is_contested() {
        const SHORT: &str = "ACGATTCAGC";
        let reads = |aligned: &str, reverse: &str, kinds: &str| {
            let forward = rows("f", 0, "ACGATT-CAGC", "ACGATTTCAGC", "ccccccccccc");
            let mut other = rows("r", 0, aligned, reverse, kinds);
            other.forward = false;
            vec![forward, other]
        };
        // The reverse read reads TT C masked but anchoring: not clean, and
        // showing no run's length, but it reads the site through without the
        // inserted T.
        let through = reads(SHORT, SHORT, "cccaaaaccc");
        let before = decided_at_with(&through, SHORT.as_bytes(), 5);
        assert_eq!(
            (before.value.as_deref(), before.state),
            (Some("TTTC"), SiteState::Called)
        );
        let after = contested_at(&through, SHORT.as_bytes(), 5);
        assert_eq!(
            (after.value.as_deref(), after.state),
            (Some("TTC"), SiteState::Contested)
        );
        assert_eq!(
            (after.supporting, after.opposing),
            (vec!["r".to_owned()], vec!["f".to_owned()])
        );
        assert!(after.runs.is_empty());
        // A reverse read that shows the insertion, even masked, confirms it.
        let shows = reads("ACGATT-CAGC", "ACGATTTCAGC", "cccaaaaaccc");
        assert_eq!(
            contested_at(&shows, SHORT.as_bytes(), 5).value.as_deref(),
            Some("TTTC")
        );
        // An unresolved call there gives the reverse strand no view.
        let unresolved = reads(SHORT, "ACGATNCAGC", "cccaaaaccc");
        assert_eq!(
            contested_at(&unresolved, SHORT.as_bytes(), 5)
                .value
                .as_deref(),
            Some("TTTC")
        );
        // A read of the same strand reading it through does not contest it.
        let mut same_strand = reads(SHORT, SHORT, "cccaaaaccc");
        same_strand[1].forward = true;
        assert_eq!(
            contested_at(&same_strand, SHORT.as_bytes(), 5)
                .value
                .as_deref(),
            Some("TTTC")
        );
    }

    // SRS-CONS-015
    #[test]
    fn a_read_views_an_insertion_site_by_what_it_shows() {
        const SITE: &str = "ACGATTCAGC";
        let view = |first: usize, aligned: &str, read: &str, kinds: &str| {
            let observed = rows("r", first, aligned, read, kinds);
            insertion::view(&observed, 4, 6, "TTTC", SITE.len())
        };
        assert_eq!(
            view(0, "ACGATT-CAGC", "ACGATTTCAGC", "aaaaaaaaaaa"),
            InsertionView::Shows
        );
        assert_eq!(
            view(0, SITE, SITE, "cccaaaaccc"),
            InsertionView::ReadsThrough
        );
        // An unresolved call, or one outside the informative interval, shows
        // no base.
        assert_eq!(
            view(0, SITE, "ACGATNCAGC", "cccaaaaccc"),
            InsertionView::Unresolved
        );
        // Reaching the site but not the position after it.
        assert_eq!(
            view(0, "ACGATT", "ACGATT", "cccccc"),
            InsertionView::ReadEnd
        );
        // Ending before the site.
        assert_eq!(view(0, "ACG", "ACG", "ccc"), InsertionView::NotCovered);
    }

    // SRS-CONS-015
    #[test]
    fn an_insertion_into_a_long_run_is_left_to_the_run_rules() {
        const LONG: &str = "ACGTTTTTCAGC";
        let forward = rows("f", 0, "ACGTTTTT-CAGC", "ACGTTTTTTCAGC", "ccccccccccccc");
        let mut other = rows("r", 0, LONG, LONG, "cccaaaaaaccc");
        other.forward = false;
        let observed = vec![forward, other];
        assert_eq!(
            contested_at(&observed, LONG.as_bytes(), 5).value.as_deref(),
            Some("TTTTTTC")
        );
    }

    /// The decision of the site or stretch holding `position` over `bases`.
    fn decided_at_with(observed: &[ReadObservations], bases: &[u8], position: usize) -> Decision {
        decide_positions(observed, bases, 4)
            .range(..=position)
            .next_back()
            .map_or_else(|| decide(&[], "A"), |(_, (_, decision))| decision.clone())
    }

    // SRS-CONS-011
    #[test]
    fn an_explained_second_on_the_other_strand_is_set_aside() {
        let mut observed = dissent("ACGTTCGTAC", "cccccccccc", 50);
        assert_eq!(decided_at(&observed, 4).opposing, ["r"]);
        // Its second base stands out, and the modality explains its peak.
        observed[1].seconds.insert(4, 'A');
        let decided = decide_positions(&observed, b"ACGTACGTAC", 4);
        let mut unexplained = observed.clone();
        unexplained[1].unexplained.insert(4);
        assert!(!set_aside_dissent(
            &mut unexplained,
            &decided,
            b"ACGTACGTAC",
            12,
            40
        ));
        let mut same_strand = observed.clone();
        same_strand[1].forward = true;
        assert!(!set_aside_dissent(
            &mut same_strand,
            &decided,
            b"ACGTACGTAC",
            12,
            40
        ));
        let mut near_end = observed.clone();
        near_end[0].edges.insert(4, 11);
        assert!(!set_aside_dissent(
            &mut near_end,
            &decided,
            b"ACGTACGTAC",
            12,
            40
        ));
        assert!(set_aside_dissent(
            &mut observed,
            &decided,
            b"ACGTACGTAC",
            12,
            40
        ));
        // Set aside, the read shows nothing at the stretch its substitution
        // made, and the position is decided alone.
        let after = decided_at(&observed, 4);
        assert_eq!(after.value.as_deref(), Some("A"));
        assert_eq!(after.state, SiteState::Called);
        assert!(after.opposing.is_empty());
    }

    // SRS-CONS-011
    #[test]
    fn a_run_length_dissent_from_a_called_variant_is_set_aside() {
        const RUN: &str = "ACATGGGGTACA";
        let reads = |opposing: &str| {
            let forward = rows("f", 0, "ACATGGGG-TACA", "ACATGGGGGTACA", "ccccccccccccc");
            let mut reverse = rows("r", 0, "ACATGGGG-TACA", "ACATGGGGGTACA", "ccccccccccccc");
            reverse.forward = false;
            let mut other = rows("o", 0, RUN, opposing, "cccccccccccc");
            // Its second base stands out, and the modality explains its peak.
            other.seconds.insert(8, 'G');
            let mut observed = vec![forward, reverse, other];
            for read in &mut observed {
                read.edges = (0..RUN.len()).map(|position| (position, 50)).collect();
            }
            observed
        };
        let set_aside = |observed: &mut Vec<ReadObservations>| {
            let decided = decide_positions(observed, RUN.as_bytes(), 4);
            set_aside_dissent(observed, &decided, RUN.as_bytes(), 12, 40)
        };
        let decided_at = |observed: &[ReadObservations]| {
            decide_positions(observed, RUN.as_bytes(), 4)
                .range(..=5)
                .next_back()
                .map(|(_, (_, decision))| decision.clone())
        };
        let mut length = reads(RUN);
        let before = decided_at(&length).map(|decision| (decision.value, decision.opposing));
        assert_eq!(
            before,
            Some((Some("GGGGGT".to_owned()), vec!["o".to_owned()]))
        );
        assert!(set_aside(&mut length));
        let after = decided_at(&length).map(|decision| (decision.value, decision.opposing));
        assert_eq!(after, Some((Some("GGGGGT".to_owned()), Vec::new())));
        // Without a second base standing out, the dissent stands.
        let mut plain = reads(RUN);
        plain[2].seconds.clear();
        assert!(!set_aside(&mut plain));
        // A read showing other runs differs in more than a length.
        let mut substituted = reads("ACATGGAGTACA");
        assert!(!set_aside(&mut substituted));
    }

    // SRS-CONS-011
    #[test]
    fn an_indel_near_its_read_end_is_set_aside() {
        let near = dissent("ACGT-CGTAC", "cccccccccc", 50);
        let decided = decide_positions(&near, b"ACGTACGTAC", 4);
        assert_eq!(decided_at(&near, 4).opposing, ["r"]);
        // Far from its ends, the deletion stands.
        assert!(!set_aside_dissent(
            &mut near.clone(),
            &decided,
            b"ACGTACGTAC",
            12,
            40
        ));
        let mut close = near.clone();
        close[1].edges.insert(4, 5);
        // The supporting read must itself lie that far from its ends.
        let mut short = close.clone();
        short[0].edges.insert(4, 39);
        assert!(!set_aside_dissent(
            &mut short,
            &decided,
            b"ACGTACGTAC",
            12,
            40
        ));
        assert!(set_aside_dissent(
            &mut close,
            &decided,
            b"ACGTACGTAC",
            12,
            40
        ));
        assert!(decided_at(&close, 4).opposing.is_empty());
        // A substitution near the end is no indel.
        let mut substituted = dissent("ACGTTCGTAC", "cccccccccc", 50);
        substituted[1].edges.insert(4, 5);
        let decided = decide_positions(&substituted, b"ACGTACGTAC", 4);
        assert!(!set_aside_dissent(
            &mut substituted,
            &decided,
            b"ACGTACGTAC",
            12,
            40
        ));
    }

    fn run(base: char, length: usize, evidence: RunEvidence) -> ConsensusRun {
        ConsensusRun {
            base,
            length,
            evidence,
        }
    }

    /// An A-tract, then the C-tract `CCCCCTCCCC` (0-based 7-16).
    const TRACT: &str = "GTCAAAACCCCCTCCCCATG";

    // SRS-CONS-004, SRS-CONS-008
    #[test]
    fn a_run_two_reads_enter_from_its_ends_takes_their_agreed_count() {
        // T12C merges the tract. The forward read enters it in phase and
        // loses phase at its end, on an unresolved call; the reverse read the
        // other way round. Both read every call of the run unmasked and count
        // ten.
        let forward = rows(
            "f",
            0,
            TRACT,
            "GTCAAAACCCCCCCCCCNTG",
            "cccccccccccccccccmmm",
        );
        let reverse = rows(
            "r",
            0,
            TRACT,
            "GTCAAANCCCCCCCCCCATG",
            "mmmmmmmccccccccccccc",
        );
        let decision = decide_stretch(&[forward.clone(), reverse.clone()], 7, 16, TRACT.as_bytes());
        assert_eq!(decision.value.as_deref(), Some("CCCCCCCCCC"));
        assert_eq!(decision.runs, [run('C', 10, RunEvidence::ConfirmedEnd)]);
        assert_eq!(decision.supporting, ["f", "r"]);
        // When the reverse read loses phase one call earlier, the counts
        // differ: the run takes the reference frame between its observed ends.
        let early = rows(
            "r",
            0,
            TRACT,
            "GTCAAANCCCCCCCCCCATG",
            "mmmmmmmmcccccccccccc",
        );
        let decision = decide_stretch(&[forward, early], 7, 16, TRACT.as_bytes());
        assert_eq!(decision.value.as_deref(), Some("CCCCCCCCCC"));
        assert_eq!(decision.runs, [run('C', 10, RunEvidence::ReferenceFrame)]);
        // When the forward read's exit call reads the run's base, masked, the
        // run goes on past its count: the count bounds the run from below.
        let going_on = rows(
            "f",
            0,
            TRACT,
            "GTCAAAACCCCCCCCCCCTG",
            "cccccccccccccccccmmm",
        );
        let decision = decide_stretch(&[going_on, reverse], 7, 16, TRACT.as_bytes());
        assert_eq!(decision.runs, [run('C', 10, RunEvidence::ReferenceFrame)]);
    }

    // SRS-CONS-004
    #[test]
    fn a_clean_end_holding_every_run_does_not_stand_for_a_read_with_more_runs() {
        // The tract CCCCC T CCCC (7-16). One read shows its run structure
        // with its C run masked. Another reads C at 7 and T at 11 unclean,
        // a C inserted after 11, then T and CCCC cleanly: its clean end, C T
        // CCCC, has the composition's runs, but the read reads two more runs
        // before them, so it shows neither where the first C run starts nor
        // its length.
        let structured = rows("s", 0, TRACT, TRACT, "cccccccmmmmmcccccccc");
        let longer = rows(
            "l",
            0,
            "GTCAAAACCCCC-TCCCCATG",
            "GTCAAAACCCCTCTCCCCATG",
            "cccccccvcccvccccccccc",
        );
        let decision = decide_stretch(&[structured, longer], 7, 16, TRACT.as_bytes());
        assert_ne!(decision.value.as_deref(), Some("CTCCCC"));
        assert!(
            decision.runs.first().is_none_or(|run| run.length != 1),
            "{decision:?}"
        );
    }

    // SRS-CONS-004
    #[test]
    fn a_read_reading_only_the_exit_cleanly_opposes_a_run_that_swallows_it() {
        // The forward read reads the whole stretch 3-17 as `AAAA` and eleven
        // C, its exit beyond 17 anchored: the tract swallows the A at 17. The
        // reverse read reads only the exit cleanly, `TCCCCA` at 12-17, in
        // phase: it contradicts the composition's last runs.
        let forward = rows(
            "f",
            0,
            TRACT,
            "GTCAAAACCCCCCCCCCCTG",
            "ccccccccccccccccccmm",
        );
        let reverse = rows("r", 0, TRACT, TRACT, "mmmmmmmmmmmmcccccccc");
        let decision = decide_stretch(&[forward.clone(), reverse], 3, 17, TRACT.as_bytes());
        assert_ne!(decision.value.as_deref(), Some("AAAACCCCCCCCCCC"));
        assert_ne!(decision.state, SiteState::Called);
        // Alone, the forward read decides the stretch.
        let alone = decide_stretch(&[forward], 3, 17, TRACT.as_bytes());
        assert_eq!(alone.value.as_deref(), Some("AAAACCCCCCCCCCC"));
    }

    // SRS-CONS-004, SRS-CONS-008
    #[test]
    fn a_read_losing_phase_in_a_run_does_not_decide_what_follows_it() {
        // The forward read reads an eleventh C in phase over the A at 17 before
        // it loses phase; the reverse read reads that A in phase.
        let forward = rows(
            "f",
            0,
            TRACT,
            "GTCAAAACCCCCCCCCCCCG",
            "ccccccccccccccccccmm",
        );
        let reverse = rows(
            "r",
            0,
            TRACT,
            "GTCAAACCCCCCCCCCCATG",
            "mmmmmmmccccccccccccc",
        );
        let decision = decide_stretch(&[forward, reverse], 7, 17, TRACT.as_bytes());
        assert_eq!(decision.value.as_deref(), Some("CCCCCCCCCCA"));
        assert_eq!(
            decision.runs,
            [
                run('C', 10, RunEvidence::ReferenceFrame),
                run('A', 1, RunEvidence::InPhase)
            ]
        );
    }

    /// A C7 run, T, and a C5 run (0-based 6-18), as in a poly-C window.
    const WINDOW: &str = "GTCAAACCCCCCCTCCCCCGCTT";
    const WINDOW_ROW: &str = "GTCAAACCCCCCC--TCCCCC-GCTT";

    // SRS-CONS-004, SRS-CONS-008
    #[test]
    fn a_run_read_until_phase_loss_takes_its_in_phase_length() {
        // A C9 T C6 haplotype; the only read loses phase where the C9 run
        // starts, so no read locates that end.
        let reverse = rows(
            "r",
            0,
            WINDOW_ROW,
            "GTCAACCCCCCCCCCTCCCCCCGCTT",
            "mmmmmmcccccccccccccccccccc",
        );
        let decision = decide_stretch(&[reverse], 6, 18, WINDOW.as_bytes());
        assert_eq!(decision.value.as_deref(), Some("CCCCCCCCCTCCCCCC"));
        assert_eq!(
            decision.runs,
            [
                run('C', 9, RunEvidence::PhaseLoss),
                run('T', 1, RunEvidence::InPhase),
                run('C', 6, RunEvidence::InPhase)
            ]
        );
    }

    // SRS-CONS-004, SRS-CONS-008
    #[test]
    fn a_lone_unresolved_call_at_a_run_end_shows_an_anchored_end() {
        // The same read, its phase lost through one unresolved call before a
        // masked A: the run's length is known only as the dominant one.
        let reverse = rows(
            "r",
            0,
            WINDOW_ROW,
            "GTCAANCCCCCCCCCTCCCCCCGCTT",
            "mmmmmmcccccccccccccccccccc",
        );
        let decision = decide_stretch(&[reverse], 6, 18, WINDOW.as_bytes());
        assert_eq!(decision.value.as_deref(), Some("CCCCCCCCCTCCCCCC"));
        assert_eq!(decision.runs[0], run('C', 9, RunEvidence::AnchoredEnd));
    }

    // SRS-CONS-004
    #[test]
    fn a_disputed_run_no_read_bounds_seeds_a_stretch() {
        // A C inserted into the C7 run (0-based 6-12) through masked calls: the
        // read shows neither run end cleanly, yet it differs on the run.
        let inserting = || {
            rows(
                "i",
                0,
                "GTCAAACCCCCCC-TCCCCCGCTT",
                "GTCAAACCCCCCCCTCCCCCGCTT",
                "ccccccmmmmmmmmmccccccccc",
            )
        };
        assert_eq!(stretches(&[inserting()], WINDOW.as_bytes(), 4), [(6, 12)]);
        // A read that shows both ends of the run cleanly bounds it.
        let bounding = rows("b", 0, WINDOW, WINDOW, &"c".repeat(WINDOW.len()));
        assert!(stretches(&[inserting(), bounding], WINDOW.as_bytes(), 4).is_empty());
    }

    // SRS-CONS-007
    #[test]
    fn an_estimated_length_equal_to_the_reference_keeps_its_site() {
        let estimated = |evidence: RunEvidence| {
            let decision = Decision {
                value: Some(TRACT[7..=16].into()),
                state: SiteState::Called,
                supporting: vec!["f".into()],
                opposing: Vec::new(),
                uninformative: Vec::new(),
                runs: vec![
                    run('C', 5, evidence),
                    run('T', 1, RunEvidence::InPhase),
                    run('C', 4, RunEvidence::InPhase),
                ],
            };
            let decided = BTreeMap::from([(7, (16, decision))]);
            let mut consensus = Consensus {
                segments: Vec::new(),
                sites: Vec::new(),
                summary: ConsensusSummary::default(),
                reference_only: BTreeSet::new(),
                unreliable_only: BTreeMap::new(),
            };
            assemble(&[7], &decided, TRACT.as_bytes(), 1, &mut consensus);
            consensus.sites.len()
        };
        assert_eq!(estimated(RunEvidence::ReferenceFrame), 1);
        assert_eq!(estimated(RunEvidence::InPhase), 0);
    }

    // SRS-CONS-004
    #[test]
    fn reads_that_place_a_run_boundary_differently_are_decided_by_position() {
        // One read deletes a T, the other substitutes it: their run lengths
        // never mix, and the read without an indel decides each position.
        let reference = "GTCAGCCCTTGACT";
        let deleting = rows("d", 0, reference, "GTCAGCCC-TGACT", "cccccccccccccc");
        let substituting = rows("s", 0, reference, "GTCAGCCCCTGACT", "cccccccccccccc");
        let decision = decide_stretch(&[deleting, substituting], 5, 9, reference.as_bytes());
        assert_eq!(decision.value.as_deref(), Some("CCCCT"));
        assert!(decision.runs.is_empty());
    }
}
