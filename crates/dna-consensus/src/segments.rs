//! Segment assembly: the decided sites split into pieces, each written as a
//! consensus segment with its sites and counts.

use std::collections::BTreeMap;

use crate::decide::Decision;
use crate::model::consensus::{
    Consensus, ConsensusSegment, ConsensusSite, ConsensusSummary, RunEvidence, SiteState,
};

/// Splits the decided sites into pieces of contiguous positions, from one
/// cleanly decided site to the next: a gap of at most `maximum_gap` undecided
/// positions stays inside a piece as `N`, while a longer gap or an uncovered
/// position splits it. Each piece lists its sites by first position.
pub(crate) fn pieces(
    decided: &BTreeMap<usize, (usize, Decision)>,
    maximum_gap: usize,
) -> Vec<Vec<usize>> {
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
pub(crate) fn assemble(
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
        let estimated = decision.runs.iter().any(|run| run.evidence.is_estimate());
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
