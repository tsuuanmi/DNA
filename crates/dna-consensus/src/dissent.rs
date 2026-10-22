//! Dissent set aside (SRS-CONS-011): the reads opposing a decision that a
//! clean read of the other strand contradicts.

use std::collections::BTreeMap;

use crate::confirm::withdraw;
use crate::decide::Decision;
use crate::model::consensus::SiteState;
use crate::observe::{Observation, ReadObservations, SiteValue, Trust, junctions};
use crate::runs::runs_of;
use dna_kernel::model::nucleotide::is_canonical;

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
pub(crate) fn set_aside_dissent(
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
        // Every opposing read is one of the sample's reads; should one not
        // be, nothing is set aside.
        debug_assert_eq!(
            opposing.len(),
            decision.opposing.len(),
            "an opposing read is not among the observed reads"
        );
        if opposing.len() != decision.opposing.len() {
            continue;
        }
        let positions = || first..=last;
        let anchored = |needed: usize| {
            supporting.iter().any(|&read| {
                let read = &observed[read];
                opposing
                    .iter()
                    .all(|&other| observed[other].forward != read.forward)
                    && positions().all(|position| {
                        read.positions
                            .get(&position)
                            .is_some_and(Observation::clean_call)
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
            for junction in junctions(first, last, length) {
                if read
                    .junctions
                    .get(&junction)
                    .is_some_and(Observation::inserts)
                {
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
        for junction in junctions(first, last, length) {
            if let Some(observation) = read.junctions.get_mut(&junction) {
                set_aside(observation);
            }
        }
    }
    withdrawn
}
