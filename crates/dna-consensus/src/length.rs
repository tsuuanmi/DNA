//! The two-strand run-length rule (PROP-0012 increment 4, "rule A",
//! SRS-CONS-013): it decides the length of a run of one base beside a
//! single-base anchor in a length mixture, from the reads' tract decodes on
//! both strands. The review confirms the runs it decides.
//!
//! Each read gives a view of a run:
//! - from its tract decode, the run it enters first (`direct`), or the next
//!   one (`derived`, the total marginal deconvolved by the direct one);
//! - or, on a strand with no decoded view, the run it reads in phase
//!   (`in_phase`).
//!
//! A view's `p` is its top length's share of its top two adjacent lengths,
//! with a standard error from the decode's own shares, errors, and margins.
//! A run is decided at length `n` when both strands view it, every view puts
//! `n` on top, and every view clears one half by more than the strands'
//! spread `d` or its strand's replicate noise `τ`, whichever is larger, by
//! `k` standard errors. Otherwise a single adjacent pair `n`/`n+1` that every
//! view holds is near-equal, and anything else is a disagreement.

use std::collections::{BTreeMap, BTreeSet};

use dna_align::model::alignment::{AlignmentColumn, Orientation};
use dna_call::model::called_read::CalledRead;
use dna_kernel::nnls::nnls;
use dna_kernel::read_evidence::{LengthMarginal, TractLengths};

use crate::ReadRun;
use crate::decide::Decision;
use crate::model::consensus::RunEvidence;
#[cfg(doc)]
use crate::read_runs_over;
use crate::runs::{anchored_pair, runs_of};

/// Lengths a derived run may take: the shifts the deconvolution fits.
const SHIFTS: std::ops::RangeInclusive<usize> = 1..=15;
/// Most columns the deconvolution's NNLS lets enter (a computation cap).
const DECONVOLUTION_ITERATIONS: usize = 200;
/// Derived shares at or below this, after normalising, are dropped.
const SHARE_FLOOR: f64 = 1e-9;

/// The rule's constants: the confidence level `k` and each strand's
/// replicate excess noise per run.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RuleParams {
    /// Standard errors a view must clear its offset by.
    pub confidence: f64,
    /// Replicate excess noise of a forward view of run 1 and run 2.
    pub forward: [f64; 2],
    /// The same for a reverse view.
    pub reverse: [f64; 2],
}

/// Where a view of a run comes from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Source {
    /// The run the read's tract enters first.
    Direct,
    /// The next run, deconvolved from the read's total.
    Derived,
    /// A run the read reads in phase.
    InPhase,
}

/// One read's view of one run.
#[derive(Debug, Clone, PartialEq)]
pub struct View {
    /// Where the view comes from.
    pub source: Source,
    /// The read is placed on the reference strand.
    pub forward: bool,
    /// The length with the largest share.
    pub top: usize,
    /// The length with the next largest share, if any.
    pub next: Option<usize>,
    /// The top length's share of its top two.
    pub p: f64,
    /// Its standard error; `None` when the decode's margin is undefined.
    pub se: Option<f64>,
}

/// A marginal's top two lengths and the statistics of their split.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Pair {
    /// The length with the largest share.
    pub top: usize,
    /// The length with the next largest share, if any.
    pub next: Option<usize>,
    /// The top's share of the two.
    pub p: f64,
    /// Its standard error; `None` when the margin is undefined.
    pub se: Option<f64>,
}

/// The positive lengths of a marginal by share descending, then length
/// ascending, with their shares and errors.
fn ordered(marginal: &LengthMarginal) -> Vec<(usize, f64, f64)> {
    let mut positive = marginal
        .lengths
        .iter()
        .zip(&marginal.shares)
        .zip(&marginal.errors)
        .filter(|((_, share), _)| **share > 0.0)
        .map(|((&length, &share), &error)| (length, share, error))
        .collect::<Vec<_>>();
    positive.sort_by(|left, right| right.1.total_cmp(&left.1).then(left.0.cmp(&right.0)));
    positive
}

/// The top pair of a decoded marginal: `p = a / (a + b)` and its standard
/// error, from the shares' errors, their difference's margin `z` (for their
/// covariance), and `z / z_adjusted` (the read's residual autocorrelation).
/// `None` when no length holds signal.
#[must_use]
pub fn pair(marginal: &LengthMarginal) -> Option<Pair> {
    let positive = ordered(marginal);
    let &(top, a, error_a) = positive.first()?;
    let Some(&(next, b, error_b)) = positive.get(1) else {
        return Some(Pair {
            top,
            next: None,
            p: 1.0,
            se: Some(0.0),
        });
    };
    let p = a / (a + b);
    let (Some(z), Some(adjusted)) = (marginal.z, marginal.z_adjusted) else {
        return Some(Pair {
            top,
            next: Some(next),
            p,
            se: None,
        });
    };
    if !(z.is_finite() && adjusted.is_finite()) || z <= 0.0 || adjusted <= 0.0 {
        return Some(Pair {
            top,
            next: Some(next),
            p,
            se: None,
        });
    }
    let (var_a, var_b) = (error_a * error_a, error_b * error_b);
    let var_difference = ((a - b) / z).powi(2);
    let covariance = (var_a + var_b - var_difference) / 2.0;
    let inflation = z / adjusted;
    let var_p = (b * b * var_a + a * a * var_b - 2.0 * a * b * covariance) / (a + b).powi(4);
    Some(Pair {
        top,
        next: Some(next),
        p,
        se: Some(var_p.max(0.0).sqrt() * inflation),
    })
}

/// The derived run's lengths and shares: the total marginal fitted, by NNLS,
/// as the first marginal shifted by each derived length plus the mark,
/// normalised; shares at or below the floor dropped. `None` when the fit
/// holds nothing.
#[must_use]
pub fn deconvolve(first: &LengthMarginal, total: &LengthMarginal) -> Option<Vec<(usize, f64)>> {
    let positive = |marginal: &LengthMarginal| {
        marginal
            .lengths
            .iter()
            .zip(&marginal.shares)
            .filter(|(_, share)| **share > 0.0)
            .map(|(&length, &share)| (length, share))
            .collect::<BTreeMap<_, _>>()
    };
    let first = positive(first);
    let total = positive(total);
    let rows = total
        .keys()
        .copied()
        .chain(
            first
                .keys()
                .flat_map(|&length| SHIFTS.map(move |shift| length + shift + 1)),
        )
        .collect::<BTreeSet<_>>();
    let design = SHIFTS
        .map(|shift| {
            rows.iter()
                .map(|&row| {
                    row.checked_sub(shift + 1)
                        .and_then(|length| first.get(&length))
                        .copied()
                        .unwrap_or(0.0)
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let observed = rows
        .iter()
        .map(|row| total.get(row).copied().unwrap_or(0.0))
        .collect::<Vec<_>>();
    let shares = nnls(&design, &observed, DECONVOLUTION_ITERATIONS);
    let sum = shares.iter().sum::<f64>();
    if sum <= 0.0 {
        return None;
    }
    let derived = SHIFTS
        .zip(shares)
        .map(|(shift, share)| (shift, share / sum))
        .filter(|(_, share)| *share > SHARE_FLOOR)
        .collect::<Vec<_>>();
    (!derived.is_empty()).then_some(derived)
}

/// The view a decode gives of the run it enters first.
#[must_use]
pub fn direct(forward: bool, first: &LengthMarginal) -> Option<View> {
    let pair = pair(first)?;
    Some(View {
        source: Source::Direct,
        forward,
        top: pair.top,
        next: pair.next,
        p: pair.p,
        se: pair.se,
    })
}

/// The view a decode gives of the next run, deconvolved; its standard error
/// is the total pair's.
#[must_use]
pub fn derived(forward: bool, first: &LengthMarginal, total: &LengthMarginal) -> Option<View> {
    let mut lengths = deconvolve(first, total)?;
    lengths.sort_by(|left, right| right.1.total_cmp(&left.1).then(left.0.cmp(&right.0)));
    let (top, a) = lengths[0];
    let next = lengths.get(1).copied();
    let total = pair(total)?;
    Some(View {
        source: Source::Derived,
        forward,
        top,
        next: next.map(|(length, _)| length),
        p: a / (a + next.map_or(0.0, |(_, b)| b)),
        se: if total.next.is_none() {
            Some(0.0)
        } else {
            total.se
        },
    })
}

/// The view a read in phase through a run gives.
#[must_use]
pub const fn in_phase(forward: bool, length: usize) -> View {
    View {
        source: Source::InPhase,
        forward,
        top: length,
        next: None,
        p: 1.0,
        se: Some(0.0),
    }
}

/// What the rule makes of a run.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Outcome {
    /// Both strands agree on this length beyond their spread and noise.
    Decided(usize),
    /// Every view holds this length and the next as its top two.
    NearEqual(usize),
    /// The strands, or their reads, view the run otherwise.
    StrandsDiffer,
    /// Only one strand views the run.
    OneStrand,
}

/// The rule's outcome for run `run` (1 or 2) from its `views`, and the
/// strands' spread `d` when the views share a top.
#[must_use]
pub fn decide(views: &[View], run: usize, params: &RuleParams) -> (Outcome, Option<f64>) {
    let forward = views.iter().filter(|view| view.forward).collect::<Vec<_>>();
    let reverse = views
        .iter()
        .filter(|view| !view.forward)
        .collect::<Vec<_>>();
    if forward.is_empty() || reverse.is_empty() {
        return (Outcome::OneStrand, None);
    }
    let pairs = |strand: &[&View]| {
        strand
            .iter()
            .filter_map(|view| {
                view.next
                    .map(|next| (view.top.min(next), view.top.max(next)))
            })
            .collect::<BTreeSet<_>>()
    };
    let (forward_pairs, reverse_pairs) = (pairs(&forward), pairs(&reverse));
    let union = forward_pairs
        .union(&reverse_pairs)
        .copied()
        .collect::<Vec<_>>();
    let flag = || match union.as_slice() {
        [(low, high)] if *high == low + 1 => Outcome::NearEqual(*low),
        _ => Outcome::StrandsDiffer,
    };
    let tops = views.iter().map(|view| view.top).collect::<BTreeSet<_>>();
    let Some(&top) = tops.iter().next().filter(|_| tops.len() == 1) else {
        let outcome = if !forward_pairs.is_empty() && forward_pairs == reverse_pairs {
            flag()
        } else {
            Outcome::StrandsDiffer
        };
        return (outcome, None);
    };
    let spread = forward
        .iter()
        .flat_map(|left| reverse.iter().map(move |right| (left.p - right.p).abs()))
        .fold(0.0, f64::max);
    let index = run.clamp(1, 2) - 1;
    let passes = views.iter().all(|view| {
        let noise = if view.forward {
            params.forward[index]
        } else {
            params.reverse[index]
        };
        view.se
            .is_some_and(|se| (view.p - 0.5) - spread.max(noise) >= params.confidence * se)
    });
    let outcome = if passes {
        Outcome::Decided(top)
    } else if union.is_empty() {
        Outcome::StrandsDiffer
    } else {
        flag()
    };
    (outcome, Some(spread))
}

/// An anchored repeat window: two runs of one base around a single anchor
/// base, 1-based reference positions, inclusive.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Window {
    /// The window's first reference position.
    pub first: usize,
    /// Its last.
    pub last: usize,
    /// The anchor's reference position.
    pub anchor: usize,
    /// The repeated base, on the reference strand.
    pub repeat: char,
    /// The anchor base, on the reference strand.
    pub anchor_base: char,
}

/// A reference-strand base on a read's strand: a reverse read's canonical
/// base is complemented, which is its own inverse.
const fn on_strand(base: char, forward: bool) -> char {
    let orientation = if forward {
        Orientation::Forward
    } else {
        Orientation::Reverse
    };
    orientation.reference_base(base)
}

/// The reference position (1-based) a read's call `call` is aligned to in
/// its alignment `columns`, or, on an inserted column, the position the
/// insertion follows on the reference, on either strand: an inserted base of
/// a run lies at the run's 3' end (SRS-ALN-012), after the run's last
/// position. The first position after it stands in only where no aligned
/// column precedes it. The two-strand rule and the review place a decode by
/// it alike.
#[must_use]
pub fn start_reference(columns: &[AlignmentColumn], call: usize) -> Option<usize> {
    let index = columns
        .iter()
        .position(|column| column.original_call_index_0based() == Some(call))?;
    let aligned = |column: &AlignmentColumn| column.reference_index_0based();
    columns[..=index]
        .iter()
        .rev()
        .find_map(aligned)
        .or_else(|| columns[index..].iter().find_map(aligned))
        .map(|position| position + 1)
}

/// The run of `window` a tract enters first: 1 before the anchor, 2 after
/// it, 1 when DNA calls one run (`single`); `None` when the tract is of
/// another base, or starts outside the window or at its anchor.
#[must_use]
pub fn entered(
    window: &Window,
    forward: bool,
    reference: Option<usize>,
    tract: &TractLengths,
    single: bool,
) -> Option<usize> {
    if tract.run != on_strand(window.repeat, forward) {
        return None;
    }
    entered_run(
        window.first,
        window.last,
        Some(window.anchor),
        reference?,
        single,
    )
}

/// The run of a window over the 1-based `first..=last` that a tract starting
/// at `position` enters first: 2 after the window's `anchor`, unless DNA
/// calls one run (`single`), and 1 otherwise, also in a window without an
/// anchor; `None` when it starts outside the window or at its anchor. The
/// two-strand rule and the review credit a decode by it alike.
#[must_use]
pub fn entered_run(
    first: usize,
    last: usize,
    anchor: Option<usize>,
    position: usize,
    single: bool,
) -> Option<usize> {
    if position < first || position > last || anchor == Some(position) {
        return None;
    }
    Some(match anchor {
        Some(anchor) if !single && position > anchor => 2,
        _ => 1,
    })
}

/// The anchor base on a read's strand.
#[must_use]
pub const fn anchor_on(window: &Window, forward: bool) -> char {
    on_strand(window.anchor_base, forward)
}

/// The rule's views of a two-run window's runs (run 1 first) from the
/// reads' tract decodes and, on a strand no decode views a run from, the
/// reads that show the run in phase (`shown_runs`, by read, as
/// [`read_runs_over`] gives them over the window).
#[must_use]
pub fn views(
    window: &Window,
    reads: &[&CalledRead],
    shown_runs: &[Option<Vec<ReadRun>>],
) -> [Vec<View>; 2] {
    let mut runs: [Vec<View>; 2] = [Vec::new(), Vec::new()];
    for read in reads {
        let forward = read.alignment.orientation() == Orientation::Forward;
        let anchor = anchor_on(window, forward);
        for tract in read.evidence.tracts() {
            let reference = start_reference(read.alignment.columns(), tract.call_start);
            let Some(run) = entered(window, forward, reference, tract, false) else {
                continue;
            };
            if tract.mark.is_some_and(|mark| mark != anchor) {
                continue;
            }
            let Some(view) = direct(forward, &tract.first) else {
                continue;
            };
            runs[run - 1].push(view);
            let next = if forward { 2 } else { 1 };
            if run != next
                && tract.mark == Some(anchor)
                && let Some(view) = derived(forward, &tract.first, &tract.total)
            {
                runs[next - 1].push(view);
            }
        }
    }
    for (read, shown) in reads.iter().zip(shown_runs) {
        let Some(shown) = shown else {
            continue;
        };
        let forward = read.alignment.orientation() == Orientation::Forward;
        let Some(index) = anchored_pair(shown, window.repeat, window.anchor_base) else {
            continue;
        };
        for (run, view) in [(0, &shown[index]), (1, &shown[index + 2])] {
            if view.in_phase
                && !runs[run]
                    .iter()
                    .any(|seen| seen.forward == forward && seen.source != Source::InPhase)
            {
                runs[run].push(in_phase(forward, view.length));
            }
        }
    }
    runs
}

/// Decides by the two-strand rule over the reads' tract decodes (SRS-CONS-013,
/// PROP-0012 increment 4) the runs a stretch's decision left as an estimate
/// (an anchored end, a phase loss, or the reference frame), in every pair of
/// runs of one base the stretch splits by a single-base anchor, where its
/// call keeps the reference's run composition. A run the rule decides takes
/// its length, with `two_strand` evidence. A run a read shows in phase, or
/// whose ends two reads confirm, keeps its length whatever the rule says;
/// the review flags their disagreement. The pattern is the stretch's own
/// structure, never a position: the rule names no region. `shown` gives each
/// read's runs over a 0-based reference stretch, as [`read_runs_over`] does.
pub(crate) fn two_strand(
    decided: &mut BTreeMap<usize, (usize, Decision)>,
    reads: &[&CalledRead],
    bases: &[u8],
    params: &RuleParams,
    shown: impl Fn(usize, usize) -> Vec<Option<Vec<ReadRun>>>,
) {
    for (&start, (end, decision)) in decided.iter_mut() {
        if !decision
            .runs
            .iter()
            .any(|run| run.evidence.yields_to_two_strand())
        {
            continue;
        }
        let reference_runs = runs_of(bases[start..=*end].iter().map(|&base| char::from(base)));
        if reference_runs.len() != decision.runs.len()
            || reference_runs
                .iter()
                .zip(&decision.runs)
                .any(|((base, _), run)| *base != run.base)
        {
            continue;
        }
        // Each reference run's first position, 1-based.
        let firsts = reference_runs
            .iter()
            .scan(start + 1, |position, (_, length)| {
                let first = *position;
                *position += length;
                Some(first)
            })
            .collect::<Vec<_>>();
        let mut changed = false;
        for index in 0..reference_runs.len().saturating_sub(2) {
            let (repeat, anchor) = (reference_runs[index], reference_runs[index + 1]);
            if reference_runs[index + 2].0 != repeat.0
                || anchor.1 != 1
                || decision.runs[index + 1].length != 1
                || !(decision.runs[index].evidence.yields_to_two_strand()
                    || decision.runs[index + 2].evidence.yields_to_two_strand())
            {
                continue;
            }
            let window = Window {
                first: firsts[index],
                last: firsts[index + 2] + reference_runs[index + 2].1 - 1,
                anchor: firsts[index + 1],
                repeat: repeat.0,
                anchor_base: anchor.0,
            };
            let shown = shown(window.first - 1, window.last - 1);
            let views = views(&window, reads, &shown);
            for (run, slot) in [(1, index), (2, index + 2)] {
                if !decision.runs[slot].evidence.yields_to_two_strand() {
                    continue;
                }
                let (Outcome::Decided(length), _) = decide(&views[run - 1], run, params) else {
                    continue;
                };
                decision.runs[slot].length = length;
                decision.runs[slot].evidence = RunEvidence::TwoStrand;
                changed = true;
            }
        }
        if changed {
            decision.value = Some(
                decision
                    .runs
                    .iter()
                    .map(|run| run.base.to_string().repeat(run.length))
                    .collect(),
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn column(base: char, call: Option<usize>, reference: Option<usize>) -> AlignmentColumn {
        AlignmentColumn::new(base, 'C', call, reference)
    }

    // SRS-API-020, SRS-CONS-013
    #[test]
    fn an_inserted_first_call_is_credited_to_the_position_it_follows() {
        // Reference 315-316 (0-based 314-315) read C G; a reverse read's run
        // of C starts, in its order, on the C inserted after 315.
        let columns = [
            column('C', Some(9), Some(313)),
            column('C', Some(8), Some(314)),
            column('C', Some(7), None),
            column('G', Some(6), Some(315)),
        ];
        assert_eq!(start_reference(&columns, 7), Some(315));
        assert_eq!(start_reference(&columns, 6), Some(316));
        // With no aligned column before it, the next one stands in.
        let leading = [column('C', Some(0), None), column('C', Some(1), Some(4))];
        assert_eq!(start_reference(&leading, 0), Some(5));
    }
}
