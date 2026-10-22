//! What one placed read observes at each reference position and junction.

use std::collections::{BTreeMap, BTreeSet};

use dna_align::model::alignment::Orientation;
use dna_call::Eligibility;
use dna_call::VariantCallingConfig;
use dna_call::model::called_read::CalledRead;
use dna_call::model::variant::VariantKind;
use dna_kernel::model::nucleotide::is_canonical;
use dna_kernel::read_evidence::{CallEvidence, EvidenceReason, Placement};

/// What one read contributes at one reference position or junction.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) enum SiteValue {
    /// A base at a reference position.
    Base(char),
    /// A deleted reference position.
    Deletion,
    /// The bases inserted at a junction; empty when the read spans the
    /// junction without an insertion.
    Insertion(String),
}

/// What one observed call can show, from least to most; each level implies
/// the ones below it. The level follows from the call's placement and its
/// verdict on its base (`ReadEvidence`):
///
/// | Placement | Verdict on its base | Level |
/// |---|---|---|
/// | unresolved, or outside the informative interval, or not canonical | any | `Nothing` |
/// | anchoring | any | `Anchoring` |
/// | scored | at an edge no bounded-run exemption (ADR-0076) lifts, or a modality reason | `Unmasked` |
/// | scored | clear of the edges (or exempt) and no modality reason | `Clean` |
///
/// A scored call whose bound the modality withholds is vouched for by its
/// base only: whatever its level, it neither counts towards a run nor shows
/// a run's end (`Char::base_only`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(super) enum Trust {
    /// Nothing: a non-canonical or uninformative call, or one beside an
    /// unresolved call.
    Nothing,
    /// Canonical and informative, masked but anchoring: it shows where a
    /// neighbouring run ends (ADR-0071).
    Anchoring,
    /// Canonical, informative, and unmasked: a call of the read's in-phase
    /// signal, which counts towards a run's length (ADR-0071).
    Unmasked,
    /// Unmasked, trusted (`read_end`), and no modality reason against its
    /// base: it shows its base.
    Clean,
}

/// One observed base: a call aligned to a position or inserted at a junction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct Char {
    pub(super) base: char,
    pub(super) trust: Trust,
    /// The reference position the base is aligned to; `None` for an inserted
    /// base.
    pub(super) position: Option<usize>,
    /// One unresolved call separates the base from a run where the read
    /// loses phase: the base shows that run's end only as `anchored_end`, the
    /// run being as long as its in-phase calls or one longer (ADR-0074).
    pub(super) beside_unresolved: bool,
    /// The modality vouches for the call's base only: it neither counts
    /// towards a run nor shows a run's end.
    pub(super) base_only: bool,
}

impl Char {
    /// The call shows its base.
    pub(super) fn clean(self) -> bool {
        self.trust == Trust::Clean
    }

    /// The call shows where a neighbouring run ends.
    pub(super) fn resolving(self) -> bool {
        self.trust >= Trust::Anchoring && !self.base_only
    }

    /// The call counts towards its run.
    pub(super) fn unmasked(self) -> bool {
        self.trust >= Trust::Unmasked && !self.base_only
    }
}

/// One read's value at one site, and whether it may decide the site.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct Observation {
    pub(super) value: SiteValue,
    /// The calls behind the value are unmasked, trusted (`read_end`), and
    /// raise no support veto.
    pub(super) clean: bool,
    /// Clean only because a reference-only call confirms the reference here
    /// (ADR-0080).
    pub(super) confirms: bool,
    /// The observed bases: one for a base, none for a deletion, the inserted
    /// ones at a junction.
    pub(super) chars: Vec<Char>,
}

impl Observation {
    /// Clean by the read's own calls, not only as a reference confirmation.
    pub(super) const fn clean_call(&self) -> bool {
        self.clean && !self.confirms
    }

    /// A junction where the read inserts bases.
    pub(super) fn inserts(&self) -> bool {
        matches!(&self.value, SiteValue::Insertion(inserted) if !inserted.is_empty())
    }
}

/// The junctions inside the 0-based stretch `first..=last` and at its two
/// edges, each keyed by the position before it, on a reference of `length`
/// positions whose origin wraps.
pub(super) fn junctions(first: usize, last: usize, length: usize) -> impl Iterator<Item = usize> {
    std::iter::once((first + length - 1) % length).chain(first..=last)
}

/// Every observation of one read, keyed by 0-based reference position; a
/// junction is keyed by the position before it.
#[derive(Debug, Clone)]
pub(super) struct ReadObservations {
    /// The read's identity, its input's content SHA-256: never its file
    /// name, which two inputs can share.
    pub(super) id: String,
    /// The read is placed on the reference strand.
    pub(super) forward: bool,
    pub(super) positions: BTreeMap<usize, Observation>,
    pub(super) junctions: BTreeMap<usize, Observation>,
    /// The second base, on the reference strand, of each position whose call
    /// shows one standing out (ADR-0080).
    pub(super) seconds: BTreeMap<usize, char>,
    /// Positions whose call shows a second peak its modality does not explain
    /// (ADR-0080).
    pub(super) unexplained: BTreeSet<usize>,
    /// Per covered position, the calls between it and the nearer end of the
    /// read's informative interval, through the nearer of the calls on both
    /// sides of a deleted position or a junction after it (SRS-CONS-011).
    pub(super) edges: BTreeMap<usize, usize>,
}

/// Walks the read's alignment columns in reference order.
///
/// A base is clean when its call is clean and canonical. A deletion is clean
/// when the nearest calls on both sides are clean. A junction between two
/// adjacent covered positions is clean when every inserted call and the
/// nearest calls on both sides are clean; it carries the inserted bases, or
/// none. An insertion or deletion that changes the length of a run the read
/// does not bound is never clean (`run_boundary`, ADR-0071): the read does
/// not show how long that run is. A reference-only call whose base is the
/// reference's makes its position's base observation clean as a
/// confirmation of the reference (ADR-0080), and nothing else: its run role
/// stays its call's.
pub(super) fn observe(
    read: &CalledRead,
    config: &VariantCallingConfig,
    bases: &[u8],
    seeded_run: usize,
) -> ReadObservations {
    let reference_length = bases.len();
    let calls = read.evidence.calls();
    let eligibility = Eligibility::new(&read.evidence, config);
    let orientation = read.alignment.orientation();
    // Calls aligned as a match, and calls aligned as a match or inserted:
    // the only bounds, and the only run calls, a bounded-run exemption may
    // trust (ADR-0076); a substituted call keeps `read_end`.
    let mut matches = vec![false; calls.len()];
    let mut length_only = vec![false; calls.len()];
    for column in read.alignment.columns() {
        if let Some(call) = column.original_call_index_0based()
            && call < calls.len()
        {
            let inserted = column.reference_index_0based().is_none();
            matches[call] = !inserted && column.query_base() == column.reference_base();
            length_only[call] = inserted || matches[call];
        }
    }
    // Scored and trusted, whatever the modality's reasons against its base.
    let trusted_call = |index: usize| {
        calls.get(index).is_some_and(|call| {
            call.placement == Placement::Scored
                && (eligibility.trusted(index)
                    || (eligibility.run_exempt(index) && length_only[index])
                    || (eligibility.bound_exempt(index) && matches[index]))
        })
    };
    // Each use reads the modality's reasons for it, as variant calling does:
    // a base its verdict on the base, an inserted base its verdict on an
    // indel, and a call beside an insertion or deletion none.
    let clean_call = |index: usize| {
        trusted_call(index)
            && calls
                .get(index)
                .is_some_and(|call| call.identity.withheld.is_empty())
    };
    let inserted_call = |index: usize| {
        trusted_call(index)
            && calls
                .get(index)
                .is_some_and(|call| call.indel.withheld.is_empty() && !call.base_only())
    };
    // A call that may bound an insertion or deletion: trusted, and not one
    // the modality vouches for by its base only.
    let placing_call = |index: usize| {
        trusted_call(index) && calls.get(index).is_some_and(|call| !call.base_only())
    };
    let interval = read.evidence.informative();
    // A base aligned to `position`, or inserted when it is `None`.
    let observed_char = |base: char, call: usize, position: Option<usize>| Char {
        base,
        trust: match calls.get(call) {
            _ if !is_canonical(base) || !interval.contains(&call) => Trust::Nothing,
            Some(evidence) => match evidence.placement {
                Placement::Scored if position.is_some() && clean_call(call) => Trust::Clean,
                Placement::Scored if position.is_none() && inserted_call(call) => Trust::Clean,
                Placement::Scored => Trust::Unmasked,
                Placement::Anchoring(_) => Trust::Anchoring,
                Placement::Unresolved(_) => Trust::Nothing,
            },
            None => Trust::Nothing,
        },
        position,
        beside_unresolved: false,
        base_only: calls.get(call).is_some_and(CallEvidence::base_only),
    };
    let columns = read.alignment.columns();
    // A call vouched for by its base only that may sit one position off.
    let reference_at = |position: usize| bases.get(position).map(|&base| char::from(base));
    let unanchored = |call: usize, base: char, position: usize| {
        let adjacent = [
            (position + reference_length - 1) % reference_length,
            (position + 1) % reference_length,
        ]
        .map(reference_at);
        Some(base) != reference_at(position) && eligibility.unanchored(call, base, adjacent)
    };
    // A call vouched for by its base only shows a position of a reference
    // homopolymer of at least `seeded_run` bases only where the read bounds
    // that run with calls trusted for more than their base: otherwise the
    // position would carry a run length no read shows.
    let placed_at = columns
        .iter()
        .filter_map(|column| {
            Some((
                column.reference_index_0based()?,
                (column.original_call_index_0based()?, column.query_base()),
            ))
        })
        .collect::<BTreeMap<_, _>>();
    let unbounded_run = |call: usize, position: usize| {
        if !calls.get(call).is_some_and(CallEvidence::base_only) {
            return false;
        }
        let Some(&base) = bases.get(position) else {
            return false;
        };
        let mut first = position;
        while first > 0 && bases[first - 1] == base {
            first -= 1;
        }
        let mut last = position;
        while last + 1 < reference_length && bases[last + 1] == base {
            last += 1;
        }
        let shows = |bound: Option<usize>| {
            bound
                .and_then(|bound| placed_at.get(&bound))
                .is_some_and(|&(call, read)| placing_call(call) && read != char::from(base))
        };
        last + 1 - first >= seeded_run && !(shows(first.checked_sub(1)) && shows(Some(last + 1)))
    };
    // Nearest call at or before, and at or after, every column.
    let mut before = Vec::with_capacity(columns.len());
    let mut last = None;
    for column in columns {
        last = column.original_call_index_0based().or(last);
        before.push(last);
    }
    let mut after = vec![None; columns.len()];
    let mut next = None;
    for (offset, column) in columns.iter().enumerate().rev() {
        next = column.original_call_index_0based().or(next);
        after[offset] = next;
    }
    let bounded = |left: Option<usize>, right: Option<usize>| {
        left.is_some_and(placing_call) && right.is_some_and(placing_call)
    };

    let mut positions = BTreeMap::new();
    let mut junctions = BTreeMap::new();
    let mut seconds = BTreeMap::new();
    let mut unexplained = BTreeSet::new();
    let mut edges = BTreeMap::new();
    // Calls between `call` and the nearer end of the informative interval.
    let edge = |call: usize| {
        if interval.contains(&call) {
            (call - interval.start).min(interval.end - 1 - call)
        } else {
            0
        }
    };
    let mut previous: Option<(usize, usize)> = None; // (reference position, column)
    let mut inserted: Vec<(char, usize)> = Vec::new();
    // Per deleted position: the nearest calls on both sides, and the base.
    let mut deletions = BTreeMap::new();
    for (offset, column) in columns.iter().enumerate() {
        let Some(position) = column.reference_index_0based() else {
            if let Some(call) = column.original_call_index_0based() {
                inserted.push((column.query_base(), call));
            }
            continue;
        };
        if let Some((previous_position, previous_column)) = previous
            && (previous_position + 1) % reference_length == position
        {
            let bases = inserted.iter().map(|(base, _)| *base).collect::<String>();
            let evidence = inserted.iter().map(|(_, call)| *call).collect::<Vec<_>>();
            let unbounded =
                eligibility.unbounded_run(VariantKind::Ins, &bases, orientation, &evidence);
            // As in variant calling: an insertion into a repeat the read
            // bounds is judged on its placements along the repeat.
            let repeat = !unbounded && eligibility.repeat_unvetoed(&evidence);
            let inserted_clean = |call: usize| {
                inserted_call(call)
                    || (repeat
                        && trusted_call(call)
                        && calls.get(call).is_some_and(|call| !call.base_only()))
            };
            let clean = !unbounded
                && inserted
                    .iter()
                    .all(|(base, call)| is_canonical(*base) && inserted_clean(*call))
                && bounded(before[previous_column], after[offset]);
            junctions.insert(
                previous_position,
                Observation {
                    value: SiteValue::Insertion(bases),
                    confirms: false,
                    clean,
                    chars: inserted
                        .iter()
                        .map(|(base, call)| {
                            let mut char = observed_char(*base, *call, None);
                            if unbounded {
                                char.trust = char.trust.min(Trust::Unmasked);
                            } else if repeat && char.trust == Trust::Unmasked && trusted_call(*call)
                            {
                                char.trust = Trust::Clean;
                            }
                            char
                        })
                        .collect(),
                },
            );
        }
        let observation = if column.query_base() == '-' {
            let flanks = (
                offset.checked_sub(1).and_then(|left| before[left]),
                after.get(offset + 1).copied().flatten(),
            );
            deletions.insert(position, (flanks, column.reference_base()));
            Observation {
                value: SiteValue::Deletion,
                confirms: false,
                clean: bounded(flanks.0, flanks.1),
                chars: Vec::new(),
            }
        } else {
            // A call that may sit one position off, or that would show a run
            // the read does not bound, shows nothing here.
            let shifted = column.original_call_index_0based().is_some_and(|call| {
                unanchored(call, column.query_base(), position) || unbounded_run(call, position)
            });
            let call = column
                .original_call_index_0based()
                .and_then(|call| calls.get(call).map(|evidence| (call, evidence)));
            if let Some((_, evidence)) = call {
                if evidence.second_unexplained {
                    unexplained.insert(position);
                }
                if let Some(second) = evidence.second.filter(|second| second.stands_out) {
                    seconds.insert(position, orientation.reference_base(second.base));
                }
            }
            let clean = is_canonical(column.query_base())
                && !shifted
                && column.original_call_index_0based().is_some_and(clean_call);
            // A call reading the reference's base confirms the reference at
            // its position (ADR-0080): a reference-only call, whose modality's
            // own edge verdict does not stand against it while the read's
            // unvouched ends still do; or a call within the margin of an
            // unvouched end that its modality vouches for its base, clear of
            // its edges and withheld for nothing, with no second peak standing
            // out or left unexplained. The margin keeps an indel or a run's
            // end from resting there; a confirmation shows neither.
            let confirms = !clean
                && is_canonical(column.query_base())
                && !shifted
                && reference_at(position) == Some(column.query_base())
                && call.is_some_and(|(call, evidence)| {
                    evidence.placement == Placement::Scored
                        && ((evidence.reference_only
                            && (trusted_call(call) || eligibility.clear_of_ends(call)))
                            || (!eligibility.clear_of_ends(call)
                                && !evidence.identity.edge
                                && evidence.identity.withheld.is_empty()
                                && !evidence.second_unexplained
                                && !evidence.second.is_some_and(|second| second.stands_out)))
                });
            Observation {
                value: SiteValue::Base(column.query_base()),
                confirms,
                clean: clean || confirms,
                chars: column
                    .original_call_index_0based()
                    .map(|call| {
                        let mut char = observed_char(column.query_base(), call, Some(position));
                        if shifted {
                            char.trust = Trust::Nothing;
                        }
                        char
                    })
                    .into_iter()
                    .collect(),
            }
        };
        positions.insert(position, observation);
        let nearest = [before[offset], after[offset]]
            .into_iter()
            .flatten()
            .map(edge)
            .min()
            .unwrap_or(0);
        edges.insert(position, nearest);
        previous = Some((position, offset));
        inserted.clear();
    }
    // Indels beyond the caller's length cap are never clean evidence; neither
    // is a deletion that changes an unbounded run's length.
    let cap = config.max_indel_length();
    for observation in junctions.values_mut() {
        if observation.chars.len() > cap {
            observation.clean = false;
            for char in &mut observation.chars {
                char.trust = char.trust.min(Trust::Unmasked);
            }
        }
    }
    let deleted = deletions.keys().copied().collect::<Vec<_>>();
    for run in deleted.chunk_by(|left, right| left + 1 == *right) {
        let edited = run
            .iter()
            .filter_map(|position| deletions.get(position).map(|(_, base)| *base))
            .collect::<String>();
        let flanking = [
            run.first()
                .and_then(|position| deletions.get(position))
                .and_then(|((left, _), _)| *left),
            run.last()
                .and_then(|position| deletions.get(position))
                .and_then(|((_, right), _)| *right),
        ]
        .into_iter()
        .flatten()
        .collect::<Vec<_>>();
        if run.len() > cap
            || eligibility.unbounded_run(VariantKind::Del, &edited, orientation, &flanking)
        {
            for position in run {
                if let Some(observation) = positions.get_mut(position) {
                    observation.clean = false;
                }
            }
        }
    }
    ReadObservations {
        id: read.input_sha256.clone(),
        forward: orientation == Orientation::Forward,
        positions,
        junctions,
        seconds,
        unexplained,
        edges,
    }
}

impl ReadObservations {
    /// The read's observed bases over `first..=last`, insertions at both
    /// edges and inside included, framed by the bases at the two neighbouring
    /// positions (`None` in the frame for a deleted neighbour); `None` when
    /// the read does not cover the stretch and both its neighbours.
    pub(super) fn span(
        &self,
        first: usize,
        last: usize,
        length: usize,
    ) -> Option<(Option<Char>, Vec<Char>, Option<Char>)> {
        let before = (first + length - 1) % length;
        let after = (last + 1) % length;
        // The base beside the stretch; when it is a lone unresolved call, the
        // masked base beyond it shows the end where the read loses phase.
        let frame = |position: usize, beyond: usize, junction: usize| -> Option<Option<Char>> {
            let char = self.positions.get(&position)?.chars.first().copied();
            if char.is_some_and(|char| char.base == 'N')
                && self
                    .junctions
                    .get(&junction)
                    .is_some_and(|junction| junction.chars.is_empty())
                && let Some(outer) = self
                    .positions
                    .get(&beyond)
                    .and_then(|observation| observation.chars.first().copied())
                && outer.trust == Trust::Anchoring
            {
                return Some(Some(Char {
                    beside_unresolved: true,
                    ..outer
                }));
            }
            Some(char)
        };
        let outward = (before + length - 1) % length;
        let (left, right) = (
            frame(before, outward, outward)?,
            frame(after, (after + 1) % length, after)?,
        );
        let mut chars = self.junctions.get(&before)?.chars.clone();
        let mut position = first;
        loop {
            let observation = self.positions.get(&position)?;
            if observation.value == SiteValue::Deletion && !observation.clean {
                // An untrusted deletion hides the base like an unresolved call.
                chars.push(Char {
                    base: 'N',
                    trust: Trust::Nothing,
                    position: Some(position),
                    beside_unresolved: false,
                    base_only: false,
                });
            }
            chars.extend(observation.chars.iter().copied());
            chars.extend(self.junctions.get(&position)?.chars.iter().copied());
            if position == last {
                break;
            }
            position = (position + 1) % length;
        }
        Some((left, chars, right))
    }

    /// The bases the read shows cleanly from one edge of `first..=last`
    /// inwards, in reference order: from `first` when `from_end` is false,
    /// from `last` otherwise, one clean canonical call per position, until a
    /// position it does not read cleanly or an insertion.
    pub(super) fn clean_edge(
        &self,
        first: usize,
        last: usize,
        length: usize,
        from_end: bool,
    ) -> String {
        let (start, stop) = if from_end {
            (last, first)
        } else {
            (first, last)
        };
        let mut bases = Vec::new();
        let mut position = start;
        while let Some(observation) = self.positions.get(&position) {
            let [char] = observation.chars.as_slice() else {
                break;
            };
            if !observation.clean || !char.clean() || !is_canonical(char.base) {
                break;
            }
            bases.push(char.base);
            if position == stop {
                break;
            }
            let (next, junction) = if from_end {
                let next = (position + length - 1) % length;
                (next, next)
            } else {
                ((position + 1) % length, position)
            };
            if self
                .junctions
                .get(&junction)
                .is_some_and(|junction| !junction.chars.is_empty())
            {
                break;
            }
            position = next;
        }
        if from_end {
            bases.reverse();
        }
        bases.into_iter().collect()
    }
}

/// The reference positions (0-based) reads observe only inside their
/// modality's unreliable spans (SRS-CONS-012): a read's call there lies in
/// such a span, and no read shows the position soundly outside its own
/// spans, by a call aligned there that is scored, trusted, and withheld for
/// nothing on its base. Each with the reason of the first read's span there.
/// A deleted position lies in a span when the calls on both sides of it do.
pub(super) fn unreliable_only(
    reads: &[&CalledRead],
    config: &VariantCallingConfig,
) -> BTreeMap<usize, EvidenceReason> {
    let mut unreliable: BTreeMap<usize, EvidenceReason> = BTreeMap::new();
    let mut sound = BTreeSet::new();
    for read in reads {
        let eligibility = Eligibility::new(&read.evidence, config);
        let spans = read.evidence.unreliable();
        let reason = |call: usize| {
            spans
                .iter()
                .find(|span| span.calls.contains(&call))
                .map(|span| span.reason)
        };
        let columns = read.alignment.columns();
        // The reason of the nearest call at or after each column.
        let mut after = vec![None; columns.len()];
        let mut next = None;
        for (offset, column) in columns.iter().enumerate().rev() {
            if let Some(call) = column.original_call_index_0based() {
                next = Some(reason(call));
            }
            after[offset] = next.flatten();
        }
        let mut before = None;
        for (offset, column) in columns.iter().enumerate() {
            let here = match column.original_call_index_0based() {
                Some(call) => {
                    before = reason(call);
                    before
                }
                None => before.and(after[offset]),
            };
            let Some(position) = column.reference_index_0based() else {
                continue;
            };
            match (here, column.original_call_index_0based()) {
                (Some(reason), _) => {
                    unreliable.entry(position).or_insert(reason);
                }
                (None, Some(call)) if eligibility.shows_base(call) => {
                    sound.insert(position);
                }
                _ => {}
            }
        }
    }
    unreliable.retain(|position, _| !sound.contains(position));
    unreliable
}
