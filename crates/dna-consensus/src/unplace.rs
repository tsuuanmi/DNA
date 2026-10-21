//! Length changes left unplaced (SRS-CONS-009): a decided insertion or
//! deletion that nothing beside it places is left undecided.

use std::collections::BTreeMap;

use crate::decide::Decision;
use crate::model::consensus::SiteState;

/// Leaves undecided every decided site that changes the reference length
/// while a position beside it is undecided or uncovered, unless the site was
/// decided by its run structure with every run's length shown by a read
/// (`in_phase` or `anchored_end`): those reads show the run's bounds
/// themselves. Otherwise the bases beside an insertion or deletion place
/// it: inside a repeat they decide which unit changes, so without them the
/// length change is not placed. Undoing one site can leave another one's
/// neighbour undecided, so this repeats until nothing changes. `length` is
/// the reference length; positions wrap.
pub(crate) fn unplace_unbounded_indels(
    decided: &mut BTreeMap<usize, (usize, Decision)>,
    length: usize,
) {
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
                    && decision.runs.iter().all(|run| run.evidence.shown_by_read());
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
