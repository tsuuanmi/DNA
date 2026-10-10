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

mod decide;
mod observe;

use std::collections::{BTreeMap, BTreeSet};

use crate::model::called_read::CalledRead;
use crate::model::consensus::{
    Consensus, ConsensusRun, ConsensusSegment, ConsensusSite, ConsensusSummary, RunEvidence,
    SiteState,
};
use crate::read_call::CoreConfig;
use decide::{Decision, decide};
use dna_kernel::model::reference::Reference;
use dna_kernel::plugin::{Contract, PluginDescriptor, PluginFamily};
use observe::{Char, ReadObservations, SiteValue, Trust, observe};

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
    name: &'a str,
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
        let mut mapped = Vec::new();
        let prefix = &self.runs[..self.prefix];
        if prefix.len() <= bases.len()
            && prefix
                .iter()
                .zip(&bases)
                .all(|(run, base)| run.base == *base)
        {
            mapped.extend(prefix.iter().map(|run| partial(run, false)).enumerate());
        }
        let suffix = &self.runs[self.runs.len() - self.suffix..];
        if let Some(offset) = bases.len().checked_sub(suffix.len())
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
    name: &'a str,
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
        name,
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
                read_runs(read.name.as_str(), left, &resolved(&chars), right, complete)
            })
        })
        .collect::<Vec<_>>();
    let names = || reads.iter().map(|read| read.name);
    let fallback = || by_position(observed, first, last, bases, names());
    let whole = reads
        .iter()
        .map(|read| (read.name, read.whole().then(|| ReadRuns::bases(&read.runs))))
        .collect::<Vec<_>>();
    let composition = if whole.iter().any(|(_, bases)| bases.is_some()) {
        decide(&whole, &reference_composition)
    } else {
        let mut joined: BTreeMap<String, BTreeSet<&str>> = BTreeMap::new();
        for start in &reads {
            for finish in &reads {
                if let Some(bases) = join(start, finish) {
                    joined
                        .entry(bases)
                        .or_default()
                        .extend([start.name, finish.name]);
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
                    .map(|read| {
                        (
                            read.name,
                            supporting.contains(read.name).then(|| bases.clone()),
                        )
                    })
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
        .map(|read| (read.name, read.mapped(&runs)))
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

    // A run no read shows the length of takes the reference frame between its
    // observed ends; without one, the longest run read in phase before a read
    // loses phase.
    for (index, (base, decision, evidence)) in decisions.iter_mut().enumerate() {
        if decision.value.is_some() || decision.state != SiteState::Unresolved {
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
/// unclean difference decides nothing and must not widen a stretch), together
/// with the neighbouring reference runs a substituted base or a deletion joins:
/// each grown
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
        for run in &decision.runs {
            match run.evidence {
                RunEvidence::ReferenceFrame => summary.reference_frame_runs += 1,
                RunEvidence::PhaseLoss => summary.phase_loss_runs += 1,
                RunEvidence::InPhase | RunEvidence::AnchoredEnd => {}
            }
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
    version: 2,
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
            .map(|&index| (index, MaskedAlignment::Anchoring, "dephased_signal"))
            .chain(
                unresolved
                    .iter()
                    .map(|&index| (index, MaskedAlignment::Unresolved, "unresolved_signal")),
            );
        for (index, alignment, reason) in masks {
            let mut call = evidence.calls()[index];
            call.mask = Some(CallMask {
                alignment,
                reason: EvidenceReason::new(reason),
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

    #[test]
    fn an_insertion_into_a_run_the_read_does_not_bound_never_decides()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        // A C inserted into the CC run at 20-21; the A before the run is call 19.
        let lengthened = "TTGACGTCAGTACGATCGTACCCTGAGTACGAGGACT";
        let bounded = build(&[&placed("b", lengthened, &[])?], &reference(), &config());
        assert_eq!(bounded.segments[0].sequence, lengthened);
        // Unresolved, the A no longer shows where the run starts, so the read
        // cannot tell the run's length and the reference length stays.
        let unbounded = build(
            &[&placed_with("u", lengthened, &[], &[19])?],
            &reference(),
            &config(),
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
                        clean: inserted.iter().all(|char| char.clean()),
                        chars: std::mem::take(&mut inserted),
                    },
                );
            }
            let observation = if read == '-' {
                Observation {
                    value: SiteValue::Deletion,
                    clean: true,
                    chars: Vec::new(),
                }
            } else {
                let char = char_of(read, kind, Some(position));
                Observation {
                    value: SiteValue::Base(read),
                    clean: char.clean(),
                    chars: vec![char],
                }
            };
            positions.insert(position, observation);
            previous = Some(position);
            position += 1;
        }
        ReadObservations {
            name: name.into(),
            positions,
            junctions,
        }
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

    #[test]
    fn a_run_whose_ends_two_reads_show_takes_the_reference_frame() {
        // T12C merges the tract. The forward read shows its start and loses
        // phase at its end; the reverse read the other way round.
        let forward = rows(
            "f",
            0,
            TRACT,
            "GTCAAAACCCCCCCCCCCTG",
            "cccccccccccccccccmmm",
        );
        let reverse = rows(
            "r",
            0,
            TRACT,
            "GTCAAACCCCCCCCCCCATG",
            "mmmmmmmccccccccccccc",
        );
        let decision = decide_stretch(&[forward, reverse], 7, 16, TRACT.as_bytes());
        assert_eq!(decision.value.as_deref(), Some("CCCCCCCCCC"));
        assert_eq!(decision.runs, [run('C', 10, RunEvidence::ReferenceFrame)]);
        assert_eq!(decision.supporting, ["f", "r"]);
    }

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
