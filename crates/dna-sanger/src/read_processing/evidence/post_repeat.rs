//! The bounded runs before a post-repeat phase loss (ADR-0076).

use crate::read_processing::SangerEvidenceConfig;
use dna_kernel::model::nucleotide::is_canonical;
use dna_kernel::read_evidence::{Bound, CallEvidence, Placement};

/// Marks every run of one base the read bounds right before it loses phase
/// after a repeat run (ADR-0076): a run of at least two canonical calls,
/// scored and counting towards runs, bounded on both sides by calls that
/// show a run's end and read another base, where on each side either more
/// than `margin` informative calls separate the bound from the nearest
/// uninformative call, or at least `post_repeat_minimum_calls` do and that
/// call is masked after a repeat run (`after_repeat`), while some call of the
/// run or its bounds is at a signal edge or within `margin` of the read's
/// ends. The core trusts such a run's length; it takes the run from the read's
/// own calls, and the alignment-dependent half of the rule stays there.
pub(super) fn mark_post_repeat_runs(
    calls: &mut [CallEvidence],
    informative: &std::ops::Range<usize>,
    after_repeat: &[bool],
    config: &SangerEvidenceConfig,
) {
    let margin = config.post_repeat_margin;
    let count = calls.len();
    let carries = |index: usize| {
        informative.contains(&index) && !matches!(calls[index].placement, Placement::Unresolved(_))
    };
    let mut before = vec![0_usize; count];
    let mut run = 0;
    for (index, slot) in before.iter_mut().enumerate() {
        run = if carries(index) { run + 1 } else { 0 };
        *slot = run;
    }
    let mut after = vec![0_usize; count];
    run = 0;
    for index in (0..count).rev() {
        run = if carries(index) { run + 1 } else { 0 };
        after[index] = run;
    }
    let trusted =
        |index: usize| !calls[index].identity.edge && index >= margin && index + margin < count;
    let resolves = |index: usize, base: char| {
        informative.contains(&index)
            && calls[index].base != base
            && calls[index].base != 'N'
            && calls[index].bound != Bound::Withheld
    };
    let side = |run: usize, edge: Option<usize>| {
        run > margin
            || (run >= config.post_repeat_minimum_calls
                && edge.is_some_and(|edge| after_repeat.get(edge).copied().unwrap_or(false)))
    };
    let mut marked = Vec::new();
    let mut start = 0;
    while start < count {
        let base = calls[start].base;
        let mut end = start;
        while end + 1 < count && calls[end + 1].base == base {
            end += 1;
        }
        if end > start
            && is_canonical(base)
            && let Some(left) = start.checked_sub(1)
            && end + 1 < count
        {
            let right = end + 1;
            let bounded = resolves(left, base)
                && resolves(right, base)
                && calls[left..=right].iter().all(|call| {
                    call.placement == Placement::Scored && call.bound != Bound::Withheld
                });
            let exempt = bounded
                && side(before[left], left.checked_sub(before[left]))
                && side(
                    after[right],
                    Some(right + after[right]).filter(|index| *index < count),
                )
                && (left..=right).any(|index| !trusted(index));
            if exempt {
                marked.push((start, end));
            }
        }
        start = end + 1;
    }
    for (start, end) in marked {
        for call in &mut calls[start..=end] {
            call.post_repeat.call = true;
        }
        calls[start - 1].post_repeat.bound = true;
        calls[end + 1].post_repeat.bound = true;
    }
}
