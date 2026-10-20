//! Reference confirmations from a read's decoded downstream sequences
//! (SRS-CONS-010, ADR-0080): after a homopolymer tract, the read is a mixture
//! of molecules whose runs differ in length but share one downstream sequence
//! `D`, which the modality decodes from the read's own signal (SRS-CALL-025).
//! Where `D` reads the reference with its flanks, outside every repeat, the
//! read confirms the reference at that position, even though its own call
//! there is masked: the base comes from the decode, a separate assertion, so
//! a masked call still shows no base of its own.
//!
//! The tract's reference end is found from the reference's own structure,
//! the run base with at most one single-base interruption, and must agree
//! with the alignment of the tract's last call; no position enters the rule.

use std::collections::{BTreeMap, BTreeSet};

use dna_align::model::alignment::{AlignmentColumn, Orientation};
use dna_call::model::called_read::CalledRead;
use dna_kernel::read_evidence::{Margin, Placement, TractLengths};

use crate::config::ConsensusConfig;
use crate::observe::{ReadObservations, SiteValue};

/// What one read's decoded downstream sequences say about reference
/// positions (0-based).
#[derive(Debug, Default, PartialEq, Eq)]
pub(super) struct Verdicts {
    /// Positions where a decode reads the reference, with its flanks,
    /// outside every repeat.
    pub(super) confirms: BTreeSet<usize>,
    /// Positions where a decode reads another base than the reference.
    pub(super) differs: BTreeSet<usize>,
}

const fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        other => other,
    }
}

/// The reference position (0-based) the read's call `call` is aligned to,
/// or, on an inserted column, that of the nearest aligned call before it in
/// the read's order.
fn aligned(read: &CalledRead, call: usize) -> Option<usize> {
    let columns = read.alignment.columns();
    let index = columns
        .iter()
        .position(|column| column.original_call_index_0based() == Some(call))?;
    match read.alignment.orientation() {
        Orientation::Forward => columns[..=index]
            .iter()
            .rev()
            .find_map(AlignmentColumn::reference_index_0based),
        Orientation::Reverse => columns[index..]
            .iter()
            .find_map(AlignmentColumn::reference_index_0based),
    }
}

/// The reference positions one step on in the read's order.
fn step(position: usize, forward: bool, length: usize) -> Option<usize> {
    if forward {
        Some(position + 1).filter(|&next| next < length)
    } else {
        position.checked_sub(1)
    }
}

/// The tract's reference interval `(first, last)` (0-based, inclusive) and
/// its end in the read's order, from `start`, where its first call is
/// aligned: the reference run of the base `run` there, or one position on
/// when the read's first run base substitutes the reference's, back over the
/// run base, and on in the read's order over the run base through at most
/// one single-base interruption followed by the run base.
fn tract_span(bases: &[u8], start: usize, run: u8, forward: bool) -> Option<(usize, usize, usize)> {
    let length = bases.len();
    let start = if bases.get(start) == Some(&run) {
        start
    } else {
        step(start, forward, length).filter(|&next| bases[next] == run)?
    };
    let mut back = start;
    while let Some(previous) = step(back, !forward, length).filter(|&at| bases[at] == run) {
        back = previous;
    }
    let mut end = start;
    let mut interrupted = false;
    loop {
        match step(end, forward, length) {
            Some(next) if bases[next] == run => end = next,
            Some(next)
                if !interrupted
                    && step(next, forward, length).is_some_and(|after| bases[after] == run) =>
            {
                interrupted = true;
                end = next;
            }
            _ => break,
        }
    }
    Some((back.min(end), back.max(end), end))
}

/// Whether reference position `position` lies in or beside a tandem repeat:
/// a homopolymer of at least `seeded_run` bases, or a period-2 or period-3
/// repeat of at least `copies` copies, other than one inside `tract`.
fn in_repeat(
    bases: &[u8],
    position: usize,
    tract: (usize, usize),
    seeded_run: usize,
    copies: usize,
) -> bool {
    [(1, seeded_run), (2, copies), (3, copies)]
        .into_iter()
        .any(|(period, copies)| {
            let need = period * (copies - 1);
            let from = position.saturating_sub(1 + need + period);
            (from..=position + 1).any(|first| {
                let last = first + need + period - 1;
                last < bases.len()
                    && (first..first + need).all(|at| bases[at] == bases[at + period])
                    && first <= position + 1
                    && position <= last + 1
                    && !(tract.0 <= first && last <= tract.1)
            })
        })
}

/// The share of each base, on the read strand, that the decoded mixture
/// puts at the call `offset` calls after the dominant molecule's tract end:
/// each molecule whose total length the decode holds reads its own
/// downstream base there, or its run while it is still in its tract. A
/// molecule that reads past the decoded span adds to no base.
fn mixture_shares(tract: &TractLengths, offset: usize) -> BTreeMap<u8, f64> {
    let downstream = tract.downstream.as_bytes();
    let mut shares = BTreeMap::new();
    let Some(dominant) = tract.total.dominant else {
        return shares;
    };
    for (&length, &share) in tract.total.lengths.iter().zip(&tract.total.shares) {
        let base = match (offset + dominant).checked_sub(length) {
            Some(at) => downstream.get(at).copied(),
            None => Some(tract.run.to_ascii_uppercase() as u8),
        };
        if let Some(base) = base.filter(|_| share > 0.0) {
            *shares.entry(base).or_insert(0.0) += share;
        }
    }
    shares
}

/// What one tract's decoded downstream sequence says, into `verdicts`.
fn tract_verdicts(
    read: &CalledRead,
    tract: &TractLengths,
    bases: &[u8],
    config: &ConsensusConfig,
    verdicts: &mut Verdicts,
) {
    let downstream = tract.downstream.as_bytes();
    if downstream.is_empty() || tract.first.margin == Margin::Undefined {
        return;
    }
    let forward = read.alignment.orientation() == Orientation::Forward;
    let run = tract.run.to_ascii_uppercase() as u8;
    let run = if forward { run } else { complement(run) };
    // The tract's end comes from the reference's own structure; the decoded
    // sequence must then read the reference right after it, over the flank,
    // which anchors the end where a dephased call's alignment cannot.
    let Some(start) = aligned(read, tract.call_start) else {
        return;
    };
    let Some((first, last, end)) = tract_span(bases, start, run, forward) else {
        return;
    };
    let decoded = downstream.len();
    let reference = |offset: usize| -> Option<(usize, u8)> {
        let position = if forward {
            end.checked_add(1 + offset)
                .filter(|&position| position < bases.len())?
        } else {
            end.checked_sub(1 + offset)?
        };
        let base = bases[position];
        Some((position, if forward { base } else { complement(base) }))
    };
    let flank = config.downstream_flank;
    let anchored = (0..=flank).all(|offset| {
        downstream
            .get(offset)
            .is_some_and(|&base| reference(offset).is_some_and(|(_, wanted)| base == wanted))
    });
    if !anchored {
        return;
    }
    let calls = read.evidence.calls();
    let truncated = tract.call_start + tract.window >= calls.len();
    // Each reference position's aligned call.
    let placed = read
        .alignment
        .columns()
        .iter()
        .filter_map(|column| {
            Some((
                column.reference_index_0based()?,
                column.original_call_index_0based()?,
            ))
        })
        .collect::<BTreeMap<_, _>>();
    // The read's own call there is the decoded mixture, as its decode
    // explains it: with a second peak, its second base's share of its top
    // two lies within `downstream_second_z` of the decode's residual sigma
    // of the share the mixture puts on that base against its first, whichever
    // base a near-equal mixture puts first; without one, its base is the
    // mixture's at the position. A call the mixture does not explain may show
    // a second sequence.
    let explained = |position: usize, offset: usize, expected: u8| {
        let Some(call) = placed.get(&position).and_then(|&call| calls.get(call)) else {
            return false;
        };
        let first = call.base.to_ascii_uppercase() as u8;
        let Some(second) = call.second else {
            return first == expected && !call.second_unexplained;
        };
        let shares = mixture_shares(tract, offset);
        let share = |base: u8| shares.get(&base).copied().unwrap_or(0.0);
        let (top, other) = (share(first), share(second.base.to_ascii_uppercase() as u8));
        if top + other <= 0.0 {
            return false;
        }
        (second.share - other / (top + other)).abs() <= config.downstream_second_z * tract.sigma
    };
    for (offset, &base) in downstream.iter().enumerate() {
        let Some((position, expected)) = reference(offset) else {
            break;
        };
        if base != expected {
            verdicts.differs.insert(position);
            continue;
        }
        let flanked = offset + flank < decoded
            && (offset.saturating_sub(flank)..=offset + flank)
                .all(|at| reference(at).is_some_and(|(_, wanted)| downstream[at] == wanted));
        if flanked
            && !truncated
            && explained(position, offset, expected)
            && !in_repeat(
                bases,
                position,
                (first, last),
                config.seeded_run,
                config.downstream_repeat_copies,
            )
        {
            verdicts.confirms.insert(position);
        }
    }
}

/// The reference positions (0-based) where the read's call is in phase
/// (scored) and shows a second peak its modality leaves unexplained: where
/// the read may show a second sequence. A masked call's second peak is its
/// dephasing and shows none.
pub(super) fn unexplained_in_phase(read: &CalledRead) -> BTreeSet<usize> {
    let calls = read.evidence.calls();
    read.alignment
        .columns()
        .iter()
        .filter_map(|column| {
            let call = calls.get(column.original_call_index_0based()?)?;
            (call.placement == Placement::Scored && call.second_unexplained)
                .then_some(column.reference_index_0based())
                .flatten()
        })
        .collect()
}

/// What every tract decode of `read` says about the reference positions
/// after its tracts.
pub(super) fn verdicts(read: &CalledRead, bases: &[u8], config: &ConsensusConfig) -> Verdicts {
    let mut verdicts = Verdicts::default();
    for tract in read.evidence.tracts() {
        tract_verdicts(read, tract, bases, config, &mut verdicts);
    }
    verdicts
}

/// Makes the read's base observation at each of `positions` a clean
/// confirmation of the reference, where it is a base the read does not show
/// cleanly itself: a decode never overrides what a read reads cleanly, and
/// the observation keeps its call's run role.
pub(super) fn confirm<'a>(
    read: &mut ReadObservations,
    positions: impl IntoIterator<Item = &'a usize>,
    bases: &[u8],
) {
    for &position in positions {
        if let Some(observation) = read.positions.get_mut(&position)
            && matches!(observation.value, SiteValue::Base(_))
            && !observation.clean
            && let Some(&reference) = bases.get(position)
        {
            let own = std::mem::replace(
                &mut observation.value,
                SiteValue::Base(char::from(reference)),
            );
            observation.clean = true;
            observation.confirms = true;
            read.decoded.insert(position, own);
        }
    }
}

/// Gives every read's withdrawn decoded confirmation back the value its own
/// call observed.
pub(super) fn restore(observed: &mut [ReadObservations]) {
    for read in observed {
        let withdrawn = read
            .decoded
            .keys()
            .copied()
            .filter(|position| {
                !read
                    .positions
                    .get(position)
                    .is_some_and(|observation| observation.confirms)
            })
            .collect::<Vec<_>>();
        for position in withdrawn {
            if let Some(own) = read.decoded.remove(&position)
                && let Some(observation) = read.positions.get_mut(&position)
            {
                observation.value = own;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use dna_kernel::read_evidence::{LengthMarginal, StructureFit};

    use super::*;

    // SRS-CONS-010
    #[test]
    fn walks_the_tract_through_one_single_base_interruption() {
        // 0-based: AAAA 0-3, CCCCC 4-8, T 9, CCCC 10-13, then ATGCC.
        let bases = b"AAAACCCCCTCCCCATGCC";
        // Forward from the run's first C: through the T to the last C.
        assert_eq!(tract_span(bases, 4, b'C', true), Some((4, 13, 13)));
        // From inside the run, back to its first C.
        assert_eq!(tract_span(bases, 6, b'C', true), Some((4, 13, 13)));
        // A first call on the A before the run, its base substituted.
        assert_eq!(tract_span(bases, 3, b'C', true), Some((4, 13, 13)));
        // Reverse from the last C: down through the T to the first C.
        assert_eq!(tract_span(bases, 13, b'C', false), Some((4, 13, 4)));
        // Two positions off the run base place nothing.
        assert_eq!(tract_span(bases, 2, b'C', true), None);
        // Two interruptions end the walk at the second.
        let twice = b"CCCTCCCACCCG";
        assert_eq!(tract_span(twice, 0, b'C', true), Some((0, 6, 6)));
    }

    // SRS-CONS-010
    #[test]
    fn finds_repeats_beside_a_position_but_not_the_tract_itself() {
        //            0123456789012345678
        let bases = b"GTACACACAGTTTTGCATG";
        // Inside and beside the AC repeat (4 copies, 2-8).
        assert!(in_repeat(bases, 5, (0, 0), 4, 3));
        assert!(in_repeat(bases, 9, (0, 0), 4, 3));
        // Beside the TTTT homopolymer (10-13).
        assert!(in_repeat(bases, 14, (0, 0), 4, 3));
        // Not when the homopolymer is the tract.
        assert!(!in_repeat(bases, 14, (10, 13), 4, 3));
        // Clear of both.
        assert!(!in_repeat(bases, 17, (0, 0), 4, 3));
        // Fewer copies than asked are no repeat.
        assert!(!in_repeat(bases, 5, (0, 0), 4, 5));
    }

    fn marginal(lengths: Vec<usize>, shares: Vec<f64>) -> LengthMarginal {
        let dominant = lengths
            .iter()
            .zip(&shares)
            .max_by(|left, right| left.1.total_cmp(right.1))
            .map(|(&length, _)| length);
        LengthMarginal {
            errors: vec![0.01; lengths.len()],
            dominant,
            next: None,
            z: None,
            z_adjusted: None,
            margin: Margin::Clear,
            lengths,
            shares,
        }
    }

    // SRS-CONS-010
    #[test]
    fn spreads_the_mixture_over_the_downstream_bases() {
        let tract = TractLengths {
            onset: 20,
            call_start: 0,
            extent: 9,
            window: 14,
            run: 'C',
            mark: None,
            candidates: vec![StructureFit {
                mark: None,
                residual: 0.01,
            }],
            residual: 0.01,
            sigma: 0.02,
            autocorrelation: 0.0,
            first: marginal(vec![8, 9], vec![0.25, 0.75]),
            total: marginal(vec![8, 9], vec![0.25, 0.75]),
            downstream: "ATGC".into(),
            cleaned: None,
        };
        // One call after the dominant (9) run's end: D[1] for the C9
        // molecules, D[2] for the C8 ones, one call further on.
        let shares = mixture_shares(&tract, 1);
        assert_eq!(shares.get(&b'T'), Some(&0.75));
        assert_eq!(shares.get(&b'G'), Some(&0.25));
        // At the first call, the C8 molecules read D[1], the C9 ones D[0].
        let shares = mixture_shares(&tract, 0);
        assert_eq!(shares.get(&b'A'), Some(&0.75));
        assert_eq!(shares.get(&b'T'), Some(&0.25));
        // A molecule reading past the decoded span adds to no base.
        let shares = mixture_shares(&tract, 3);
        assert_eq!(shares.get(&b'C'), Some(&0.75));
        assert_eq!(shares.values().sum::<f64>(), 0.75);
    }
}
