//! Stretches decided by their run structure (ADR-0074): where a stretch is
//! seeded, each read's runs over it, and its composition and length votes.

use std::collections::{BTreeMap, BTreeSet};

use crate::decide::{Decision, decide};
use crate::model::consensus::{ConsensusRun, RunEvidence, SiteState};
use crate::observe::{Char, ReadObservations, SiteValue, Trust, junctions};

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
pub(crate) struct RunView {
    pub(crate) base: char,
    pub(crate) length: usize,
    /// Every call of the run is clean: the read shows the run's base.
    pub(crate) clean: bool,
    /// Every call of the run is unmasked: the calls count the run.
    pub(crate) counted: bool,
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
    /// The read shows both of the run's ends in phase.
    pub(crate) fn ends_in_phase(&self) -> bool {
        self.left == End::InPhase && self.right == End::InPhase
    }

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
pub(crate) struct ReadRuns<'a> {
    /// The read's identity, its input's content SHA-256.
    id: &'a str,
    pub(crate) runs: Vec<RunView>,
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
pub(crate) fn read_runs<'a>(
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
pub(crate) fn decide_stretch(
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
            .map(|(index, (base, length))| ConsensusRun {
                base,
                length,
                evidence: bounded_evidence(&mapped, &decision.supporting, index),
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
            let evidence = bounded_evidence(&mapped, &decision.supporting, index);
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

/// The evidence of run `index`'s length voted from reads that show both of
/// its ends: `in_phase` when a `supporting` read shows both in phase,
/// `anchored_end` otherwise.
fn bounded_evidence(
    mapped: &[(&str, Vec<(usize, RunView)>)],
    supporting: &[String],
    index: usize,
) -> RunEvidence {
    let in_phase = mapped.iter().any(|(name, views)| {
        supporting.iter().any(|read| read == name)
            && views
                .iter()
                .any(|(mapped, view)| *mapped == index && view.ends_in_phase())
    });
    if in_phase {
        RunEvidence::InPhase
    } else {
        RunEvidence::AnchoredEnd
    }
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
        }) && junctions(first, last, length).all(|position| {
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
pub(crate) fn resolved(chars: &[Char]) -> Vec<Char> {
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
pub(crate) fn runs_of(bases: impl Iterator<Item = char>) -> Vec<(char, usize)> {
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
pub(crate) fn stretches(
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
