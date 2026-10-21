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

mod confirm;
mod decide;
pub mod decided;
mod dissent;
mod insertion;
pub mod length;
pub mod noise;
mod observe;
mod segments;
mod stretch;
mod unplace;

use std::collections::{BTreeMap, BTreeSet};

use crate::model::consensus::{Consensus, ConsensusSummary};
use confirm::{
    confirmed_only, hide_confirmations, withdraw_beside_observed_differences,
    withdraw_confirmations, withdraw_near_differences, withdraw_unobserved_run_junctions,
};
use decide::{Decision, decide};
use dissent::set_aside_dissent;
use dna_call::model::called_read::CalledRead;
use dna_call::read_call::CallConfig;
use dna_kernel::model::reference::Reference;
use dna_kernel::plugin::{Contract, PluginDescriptor, PluginFamily};
use length::two_strand;
use observe::{ReadObservations, SiteValue, observe};
use segments::{assemble, pieces};
use stretch::{decide_stretch, read_runs, resolved, stretches};
use unplace::unplace_unbounded_indels;

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
                        in_phase: run.clean && run.counted && run.ends_in_phase(),
                    })
                    .collect(),
            )
        })
        .collect()
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
    withdraw_confirmations(&mut strict, |_| true);
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
mod tests;
