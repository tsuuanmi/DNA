//! Insertions both strands must show (SRS-CONS-015, ADR-0089).
//!
//! An insertion that adds to a short repeat is a composition question that
//! either strand can answer: where both strands cover its site, it is decided
//! only when no read of the strand without clean support reads the site
//! through without it. A long repeat's length is left to the run-structure
//! rules, since a read of the other strand loses phase there and often reads
//! through it.

use std::collections::BTreeMap;

use dna_call::model::called_read::CalledRead;
use dna_call::read_call::CallConfig;
use dna_kernel::model::reference::Reference;

use crate::config::ConsensusConfig;
use crate::decide::Decision;
use crate::model::consensus::SiteState;
use crate::observe::{ReadObservations, Trust, observe};

/// How one read shows a site whose decision inserts bases.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum InsertionView {
    /// The read covers no position of the site.
    NotCovered,
    /// The read covers part of the site, or not both positions beside it: it
    /// ends there.
    ReadEnd,
    /// The read covers the site and both positions beside it, but a call
    /// there shows no base: unresolved, outside the read's informative
    /// interval, or a deleted neighbour.
    Unresolved,
    /// Every call over the site and beside it shows a base, in phase or
    /// masked but anchoring, and together they read another sequence than
    /// the decision: the read reads the site through without the insertion.
    ReadsThrough,
    /// Every call shows a base, and together they read the decision.
    Shows,
}

impl InsertionView {
    /// The view's `snake_case` label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::NotCovered => "not_covered",
            Self::ReadEnd => "read_end",
            Self::Unresolved => "unresolved",
            Self::ReadsThrough => "reads_through",
            Self::Shows => "shows",
        }
    }
}

/// How `read` shows the 0-based site `first..=last` decided as `call`.
pub(crate) fn view(
    read: &ReadObservations,
    first: usize,
    last: usize,
    call: &str,
    length: usize,
) -> InsertionView {
    if !(first..=last).any(|position| read.positions.contains_key(&position)) {
        return InsertionView::NotCovered;
    }
    let Some((left, chars, right)) = read.span(first, last, length) else {
        return InsertionView::ReadEnd;
    };
    let shows_base = |trust: Trust, base: char| trust >= Trust::Anchoring && base != 'N';
    let resolved = [left, right]
        .iter()
        .all(|frame| frame.is_some_and(|char| shows_base(char.trust, char.base)))
        && chars.iter().all(|char| shows_base(char.trust, char.base));
    if !resolved {
        return InsertionView::Unresolved;
    }
    if chars.iter().map(|char| char.base).eq(call.chars()) {
        InsertionView::Shows
    } else {
        InsertionView::ReadsThrough
    }
}

/// Whether `call` decides the 0-based reference site `first..=last` as its
/// reference with one block of bases inserted, where the site holds no
/// reference homopolymer of `seeded_run` bases or more and the block adds to
/// a tandem repeat of its unit (the shortest one it repeats) of fewer than
/// `seeded_run` reference copies: on such a short repeat, a read that reads
/// the site in phase or anchored shows whether the block is there. A long
/// run or repeat is a length question a read of the other strand often
/// reads out of phase; its run rules decide it (SRS-CONS-004).
pub(crate) fn short_insertion(
    bases: &[u8],
    first: usize,
    last: usize,
    call: &str,
    seeded_run: usize,
) -> bool {
    let (Some(reference), call) = (bases.get(first..=last), call.as_bytes()) else {
        return false;
    };
    if call.len() <= reference.len()
        || reference
            .chunk_by(|left, right| left == right)
            .any(|run| run.len() >= seeded_run)
    {
        return false;
    }
    let prefix = reference
        .iter()
        .zip(call)
        .take_while(|(left, right)| left == right)
        .count();
    let suffix = reference[prefix..]
        .iter()
        .rev()
        .zip(call[prefix..].iter().rev())
        .take_while(|(left, right)| left == right)
        .count();
    if prefix + suffix != reference.len() {
        return false;
    }
    let inserted = &call[prefix..call.len() - suffix];
    let unit = (1..=inserted.len())
        .find(|&period| {
            inserted.len() % period == 0
                && inserted
                    .iter()
                    .enumerate()
                    .all(|(index, base)| *base == inserted[index % period])
        })
        .unwrap_or(inserted.len());
    // The repeat around the block, in the haplotype the call decides.
    let reach = seeded_run * unit + unit;
    let start = first.saturating_sub(reach);
    let end = (last + 1 + reach).min(bases.len());
    let haplotype = [&bases[start..first], call, &bases[last + 1..end]].concat();
    let mut from = first - start + prefix;
    let mut to = from + inserted.len();
    while from > 0 && haplotype[from - 1] == haplotype[from - 1 + unit] {
        from -= 1;
    }
    while to < haplotype.len() && haplotype[to] == haplotype[to - unit] {
        to += 1;
    }
    (to - from - inserted.len()) / unit < seeded_run
}

/// Decides as the reference, `contested`, every site called with a short
/// insertion (`short_insertion`) whose clean support lies on one strand,
/// where a read of the other strand reads the site through without it and
/// none of that strand shows it (SRS-CONS-015). The reads that read it
/// through support the reference, and the insertion's reads oppose it.
pub(crate) fn contest_one_strand(
    decided: &mut BTreeMap<usize, (usize, Decision)>,
    observed: &[ReadObservations],
    bases: &[u8],
    seeded_run: usize,
) {
    let forward = observed
        .iter()
        .map(|read| (read.id.as_str(), read.forward))
        .collect::<BTreeMap<_, _>>();
    for (&first, (last, decision)) in decided.iter_mut() {
        let Some(call) = decision.value.as_deref() else {
            continue;
        };
        if decision.state != SiteState::Called
            || !short_insertion(bases, first, *last, call, seeded_run)
        {
            continue;
        }
        let strands = decision
            .supporting
            .iter()
            .filter_map(|name| forward.get(name.as_str()).copied())
            .collect::<Vec<_>>();
        let Some(&strand) = strands.first() else {
            continue;
        };
        if strands.iter().any(|&other| other != strand) {
            continue;
        }
        let views = observed
            .iter()
            .filter(|read| read.forward != strand)
            .map(|read| (read, view(read, first, *last, call, bases.len())))
            .collect::<Vec<_>>();
        if views.iter().any(|(_, view)| *view == InsertionView::Shows) {
            continue;
        }
        let through = views
            .iter()
            .filter(|(_, view)| *view == InsertionView::ReadsThrough)
            .map(|(read, _)| read.id.clone())
            .collect::<Vec<_>>();
        if through.is_empty() {
            continue;
        }
        let reference = String::from_utf8_lossy(&bases[first..=*last]).into_owned();
        let mut opposing = std::mem::take(&mut decision.supporting);
        opposing.extend(
            decision
                .opposing
                .drain(..)
                .filter(|name| !through.contains(name)),
        );
        decision
            .uninformative
            .retain(|name| !through.contains(name));
        opposing.sort_by_key(|name| observed.iter().position(|read| &read.id == name));
        decision.value = Some(reference);
        decision.state = SiteState::Contested;
        decision.supporting = through;
        decision.opposing = opposing;
        decision.runs.clear();
    }
}

/// Each read's view of the 0-based reference site `first..=last` decided as
/// `call`, by the consensus's own reading of the read (SRS-CONS-015); reads
/// come in `reads` order.
#[must_use]
pub fn insertion_views(
    reads: &[&CalledRead],
    reference: &Reference,
    config: &CallConfig,
    consensus_config: &ConsensusConfig,
    first: usize,
    last: usize,
    call: &str,
) -> Vec<InsertionView> {
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
            view(&observed, first, last, call, bases.len())
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    // SRS-CONS-015
    #[test]
    fn an_insertion_into_a_short_repeat_is_short() {
        // TTC -> TTTC: two reference copies of T.
        let bases = b"GACTCTGTTCTTTCAG";
        assert!(short_insertion(bases, 7, 9, "TTTC", 4));
        // A base that repeats nothing: one copy beside it.
        assert!(short_insertion(bases, 1, 2, "AAC", 4));
        // No insertion, or one that substitutes too.
        assert!(!short_insertion(bases, 7, 9, "TTC", 4));
        assert!(!short_insertion(bases, 7, 9, "AATC", 4));
    }

    // SRS-CONS-015
    #[test]
    fn an_insertion_into_a_long_repeat_is_not_short() {
        // AAAAAA -> AAAAAAA: six copies.
        let bases = b"GCAAAAAATTTC";
        assert!(!short_insertion(bases, 2, 7, "AAAAAAA", 4));
        // An AC unit into five AC copies, placed at either end of the repeat.
        let repeat = b"GCAGCACACACACACCGCT";
        assert!(!short_insertion(repeat, 14, 15, "CACC", 4));
        assert!(!short_insertion(repeat, 4, 5, "CACA", 4));
        // A unit of three copies stays short.
        let short = b"GCACACACCGCT";
        assert!(short_insertion(short, 7, 8, "CACC", 4));
        // A base that splits a run of five C: the run's own question.
        let run = b"GATCCCCCGAT";
        assert!(!short_insertion(run, 3, 7, "CCCTCC", 4));
    }
}
