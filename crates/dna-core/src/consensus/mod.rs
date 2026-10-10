//! Sample consensus: adjudicates between placed reads and assembles one
//! consensus sequence per covered reference interval (PROP-0003, ADR-0073).
//!
//! Each read observes the reference positions it covers and the junctions
//! between them. Wherever a read cleanly differs from the reference, the site
//! is grown to whole reference runs and decided as one stretch (as a whole,
//! run by run, or position by position), so reads that describe one
//! haplotype through different alignments agree. Everywhere else each
//! position is decided on its own. Only clean observations (unmasked,
//! trusted, no support veto) decide, under decision rule version 1.

mod decide;
mod observe;

use std::collections::{BTreeMap, BTreeSet};

use crate::model::called_read::CalledRead;
use crate::model::consensus::{
    Consensus, ConsensusSegment, ConsensusSite, ConsensusSummary, SiteState,
};
use crate::read_call::CoreConfig;
use decide::{Decision, decide};
use dna_kernel::model::reference::Reference;
use dna_kernel::plugin::{Contract, PluginDescriptor, PluginFamily};
use observe::{Char, ReadObservations, SiteValue, observe};

/// Builds the consensus of the sample's admitted reads.
///
/// Reads are expected in the sample's deterministic read order; the reads
/// listed at each site follow it.
#[must_use]
pub fn build(reads: &[&CalledRead], reference: &Reference, config: &CoreConfig) -> Consensus {
    let length = reference.len();
    let bases = reference.sequence.as_bytes();
    let observed = reads
        .iter()
        .map(|read| observe(read, &config.variant_calling, length))
        .collect::<Vec<_>>();
    let covered = observed
        .iter()
        .flat_map(|read| read.positions.keys().copied())
        .collect::<BTreeSet<_>>();

    let mut decided: BTreeMap<usize, (usize, Decision)> = BTreeMap::new();
    for (first, last) in stretches(&observed, bases) {
        decided.insert(first, (last, decide_stretch(&observed, first, last, bases)));
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
                    (read.name.as_str(), observation.clean.then_some(sequence))
                })
            })
            .collect::<Vec<_>>();
        let expected = char::from(bases[position]).to_string();
        decided.insert(position, (position, decide(&observations, &expected)));
    }

    let mut consensus = Consensus {
        segments: Vec::new(),
        sites: Vec::new(),
        summary: ConsensusSummary::default(),
    };
    let minimum = config.alignment.minimum_callable_bases;
    for piece in pieces(&decided, minimum) {
        assemble(&piece, &decided, bases, minimum, &mut consensus);
    }
    consensus
}

/// One read over a stretch: its name, the base before the stretch, its
/// resolved bases inside, the base after it, and whether no base inside was
/// unresolved.
type Span<'a> = (&'a str, Option<Char>, Vec<Char>, Option<Char>, bool);

/// Decides one stretch.
///
/// When a read observes the whole stretch cleanly with another run structure
/// than the reference (a substitution inside a run, a new run), the stretch
/// is decided as a whole from such whole observations. Otherwise every
/// reference run is decided on its own: a read with the reference run
/// structure contributes a run's length when the run's calls are clean and
/// the bases on both sides show where it ends (ADR-0071), so one read can
/// decide one run and another read the next.
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
    let expected_runs = runs_of(expected.chars());
    let spans: Vec<Span<'_>> = observed
        .iter()
        .filter_map(|read| {
            read.span(first, last, length).map(|(left, chars, right)| {
                let complete = chars.iter().all(|char| char.base != 'N');
                (read.name.as_str(), left, resolved(chars), right, complete)
            })
        })
        .collect();
    let resolving = |frame: &Option<Char>| frame.is_some_and(|char| char.resolving);
    let structured = |chars: &[Char]| {
        let observed = runs_of(chars.iter().map(|char| char.base));
        observed.len() == expected_runs.len()
            && observed
                .iter()
                .zip(&expected_runs)
                .all(|(read, wanted)| read.0 == wanted.0)
    };
    let whole = spans
        .iter()
        .map(|(name, left, chars, right, complete)| {
            let clean = *complete
                && chars.iter().all(|char| char.clean)
                && resolving(left)
                && resolving(right);
            (
                *name,
                clean.then(|| chars.iter().map(|char| char.base).collect::<String>()),
            )
        })
        .collect::<Vec<_>>();
    if whole
        .iter()
        .zip(&spans)
        .any(|((_, sequence), (_, _, chars, _, _))| sequence.is_some() && !structured(chars))
    {
        return decide(&whole, &expected);
    }
    let names = || spans.iter().map(|(name, ..)| *name);
    if !spans.iter().any(|(_, _, chars, _, _)| structured(chars)) {
        return by_position(observed, first, last, bases, names());
    }
    let decisions = expected_runs
        .iter()
        .enumerate()
        .map(|(index, &(base, reference_length))| {
            let observations = spans
                .iter()
                .map(|(name, left, chars, right, _)| {
                    if !structured(chars) {
                        return (*name, None);
                    }
                    let groups = groups_of(chars);
                    let run = &groups[index];
                    let left = if index == 0 {
                        *left
                    } else {
                        groups[index - 1].last().copied()
                    };
                    let right = groups
                        .get(index + 1)
                        .map_or(*right, |next| next.first().copied());
                    let clean =
                        run.iter().all(|char| char.clean) && resolving(&left) && resolving(&right);
                    (*name, clean.then(|| base.to_string().repeat(run.len())))
                })
                .collect::<Vec<_>>();
            (
                base,
                reference_length,
                decide(&observations, &base.to_string().repeat(reference_length)),
            )
        })
        .collect::<Vec<_>>();
    let combined = combine(&decisions, names());
    if combined.state == SiteState::Unresolved {
        return by_position(observed, first, last, bases, names());
    }
    combined
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
                    .is_none_or(|observation| observation.chars.iter().all(|char| !char.clean))
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
                        (read.name.as_str(), call)
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
    }
}

/// The observed bases without unresolved calls (`N`). The bases next to a
/// dropped call can neither be clean nor show where a run ends, because the
/// unresolved call may extend either run.
fn resolved(chars: Vec<Char>) -> Vec<Char> {
    let mut kept: Vec<Char> = Vec::with_capacity(chars.len());
    let mut after_unresolved = false;
    for mut char in chars {
        if char.base == 'N' {
            if let Some(previous) = kept.last_mut() {
                previous.clean = false;
                previous.resolving = false;
            }
            after_unresolved = true;
            continue;
        }
        if after_unresolved {
            char.clean = false;
            char.resolving = false;
            after_unresolved = false;
        }
        kept.push(char);
    }
    kept
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
/// unclean difference decides nothing and must not widen a stretch): each grown
/// to whole reference runs on both sides, overlapping or touching intervals
/// merged. A stretch never crosses the reference origin.
fn stretches(observed: &[ReadObservations], bases: &[u8]) -> Vec<(usize, usize)> {
    let length = bases.len();
    let mut seeds = BTreeSet::new();
    for read in observed {
        for (&position, observation) in &read.positions {
            if observation.clean
                && observation.value != SiteValue::Base(char::from(bases[position]))
            {
                seeds.insert(position);
            }
        }
        for (&position, observation) in &read.junctions {
            if observation.chars.iter().any(|char| char.clean) && position + 1 < length {
                seeds.insert(position);
                seeds.insert(position + 1);
            }
        }
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
/// cleanly decided site to the next: a gap of at most `minimum` undecided
/// positions stays inside a piece as `N`, while a longer gap or an uncovered
/// position splits it. Each piece lists its sites by first position.
fn pieces(decided: &BTreeMap<usize, (usize, Decision)>, minimum: usize) -> Vec<Vec<usize>> {
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
        let distant = observed && observed_last.is_some_and(|end| first - end - 1 > minimum);
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
        match decision.state {
            SiteState::Called => summary.called_positions += span,
            SiteState::Contested => summary.contested_sites += 1,
            SiteState::Unresolved => {}
        }
        let reference = String::from_utf8_lossy(&bases[*first..=*last]).into_owned();
        if decision.value.as_deref() != Some(reference.as_str())
            || decision.state != SiteState::Called
            || !decision.opposing.is_empty()
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
}

/// The sample consensus plugin.
pub const PLUGIN: PluginDescriptor = PluginDescriptor {
    id: "consensus",
    family: PluginFamily::Core,
    version: 1,
    provides: &[Contract::Consensus],
    requires: &[Contract::CalledVariants],
    config_sections: &[],
};

#[cfg(test)]
mod tests {
    use crate::alignment::{AlignmentConfig, align_best};
    use crate::model::variant::VariantCallingResult;
    use dna_kernel::error::Result;
    use dna_kernel::model::reference::ReferenceTopology;
    use dna_kernel::read_evidence::{CallMask, EvidenceReason, MaskedAlignment, ReadEvidence};

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

    fn alignment() -> AlignmentConfig {
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

    fn config() -> CoreConfig {
        CoreConfig {
            alignment: alignment(),
            variant_calling: crate::variant_calling::VariantCallingConfig {
                max_indel_length: 50,
                read_end_margin: 0,
            },
            sample_reconciliation: crate::sample::SampleReconciliationConfig {
                minimum_comparable_bases: 1,
                minimum_overlap_agreement: 0.5,
            },
        }
    }

    /// A read of `sequence` placed on the reference, with the calls in
    /// `masked` masked as dephased.
    fn placed(name: &str, sequence: &str, masked: &[usize]) -> Result<CalledRead> {
        let mut evidence = ReadEvidence::clean(sequence);
        for &index in masked {
            let mut call = evidence.calls()[index];
            call.mask = Some(CallMask {
                alignment: MaskedAlignment::Anchoring,
                reason: EvidenceReason::new("dephased_signal"),
            });
            evidence = evidence.with_call(index, call);
        }
        let alignment = align_best(&evidence, &reference(), &alignment())?;
        Ok(CalledRead {
            input_name: name.into(),
            input_sha256: name.into(),
            reference_sha256: String::new(),
            configuration_sha256: String::new(),
            evidence,
            alignment,
            variants: VariantCallingResult {
                reported: Vec::new(),
                observed: Vec::new(),
                excluded: Vec::new(),
            },
        })
    }

    fn site_at(consensus: &Consensus, position: usize) -> Option<&ConsensusSite> {
        consensus
            .sites
            .iter()
            .find(|site| (site.first_0based..=site.last_0based).contains(&position))
    }

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
        let consensus = build(&[&forward, &reverse], &reference(), &config());
        assert_eq!(consensus.segments.len(), 1);
        assert_eq!(consensus.segments[0].sequence, haplotype);
        let substitution = site_at(&consensus, 10).ok_or("missing substitution site")?;
        assert_eq!(substitution.call.as_deref(), Some("G"));
        assert_eq!(substitution.state, SiteState::Called);
        assert_eq!(substitution.supporting, ["f", "r"]);
        // The insertion is decided over the whole CC run it lengthens.
        let insertion = site_at(&consensus, 20).ok_or("missing insertion site")?;
        assert_eq!((insertion.first_0based, insertion.last_0based), (20, 21));
        assert_eq!(insertion.call.as_deref(), Some("CGGC"));
        Ok(())
    }

    #[test]
    fn a_lone_read_against_a_reference_read_is_contested()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        let variant = placed("v", "TTGACGTCAGGACGATCGTACCTGAGTACGAGGACT", &[])?;
        let plain = placed("p", REFERENCE, &[])?;
        let consensus = build(&[&variant, &plain], &reference(), &config());
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

    #[test]
    fn masked_calls_never_decide() -> std::result::Result<(), Box<dyn std::error::Error>> {
        // The variant read's call at 10 is masked: the reference read decides.
        let variant = placed("v", "TTGACGTCAGGACGATCGTACCTGAGTACGAGGACT", &[10])?;
        let plain = placed("p", REFERENCE, &[])?;
        let consensus = build(&[&variant, &plain], &reference(), &config());
        assert_eq!(consensus.segments[0].sequence, REFERENCE);
        assert!(site_at(&consensus, 10).is_none());
        // Alone, a masked call leaves the position unresolved.
        let alone = build(&[&variant], &reference(), &config());
        assert_eq!(&alone.segments[0].sequence[8..12], "AGNA");
        assert_eq!(alone.summary.unresolved_positions, 1);
        Ok(())
    }
}
