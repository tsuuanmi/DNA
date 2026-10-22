//! The read's reference-only calls (SRS-CALL-022, ADR-0080).

use super::{MIXED_SUPPORTING_DNA, WEAK_SIGNAL, quality_vetoes_only};
use crate::model::basecalls::BaseCalls;
use crate::model::callability::ReadCallability;
use dna_kernel::read_evidence::{CallEvidence, Placement};

/// Marks the read's other reference-only calls (SRS-CALL-022, ADR-0080), once
/// its signal edges are known, anywhere in the read. A scored single-channel
/// call without an unexplained second peak is reference-only when it is
/// withheld only by the peak floor or relative quality, or lies at a signal
/// edge every cause of which is benign: the untrimmed start of the read, or a
/// weak-signal mask clear of every repeat run by more than the margin (beside
/// a repeat run a weak segment is the dephased tail). Its second base may
/// stand out when its cause is a neighbour's tail or a dye blob: the call's
/// own signal explains it. The untrimmed 3' end is not benign: there the
/// read's signal decays; nor is the read's physical end. A mixed
/// reference-only call keeps the mark at a signal edge only when the edge is
/// benign by the same test.
pub(super) fn reference_only_calls(
    evidence: &mut [CallEvidence],
    calls: &BaseCalls,
    callability: &ReadCallability,
    (start, end): (usize, usize),
    margin: usize,
) {
    let count = evidence.len();
    let informative = |evidence: &[CallEvidence], at: usize| {
        (start..end).contains(&at) && !matches!(evidence[at].placement, Placement::Unresolved(_))
    };
    // A weak segment is benign only clear of every repeat run: beside one it
    // is the read's dephased tail, not a weak start.
    let clear_of_repeats = |at: usize| {
        callability.segment_at(at).is_some_and(|segment| {
            callability.repeats.iter().all(|repeat| {
                repeat.call_end_0based_exclusive + margin <= segment.call_start_0based
                    || segment.call_end_0based_exclusive + margin <= repeat.call_start_0based
            })
        })
    };
    let benign = |evidence: &[CallEvidence], at: usize| {
        at < start
            || (evidence[at].placement == Placement::Unresolved(WEAK_SIGNAL)
                && clear_of_repeats(at))
    };
    // Whether every cause of the call's edge is benign: the nearest
    // uninformative call within the margin on each side, or else the read's
    // physical end, which never is.
    let benign_edge = |evidence: &[CallEvidence], index: usize| {
        let left = (index.saturating_sub(margin)..index)
            .rev()
            .find(|&at| !informative(evidence, at));
        let right =
            (index + 1..(index + margin + 1).min(count)).find(|&at| !informative(evidence, at));
        let left_ok = match left {
            Some(at) => benign(evidence, at),
            None => index >= margin,
        };
        let right_ok = match right {
            Some(at) => benign(evidence, at),
            None => index + margin < count,
        };
        left_ok && right_ok
    };
    for index in 0..count {
        let call = evidence[index];
        if call.identity.withheld.contains(MIXED_SUPPORTING_DNA) {
            if call.reference_only && call.identity.edge && !benign_edge(evidence, index) {
                evidence[index].reference_only = false;
            }
            continue;
        }
        let single = calls
            .calls
            .get(index)
            .is_some_and(|base| base.qualifying_channels.len() == 1);
        // A second base may stand out: unless its peak is unexplained, it is
        // a neighbour's tail or a dye blob.
        if call.placement != Placement::Scored
            || !single
            || !quality_vetoes_only(&call.identity.withheld)
            || call.second_unexplained
            || call.identity.trusted()
            || (call.identity.edge && !benign_edge(evidence, index))
        {
            continue;
        }
        evidence[index].reference_only = true;
    }
}
