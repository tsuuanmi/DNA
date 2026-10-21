//! Run lengths in the profile's anchored repeat windows and length tracts,
//! judged from the reads' tract decodes (SRS-API-020, PROP-0012 increment 3).
//!
//! An anchored window holds two runs of one base around an anchor base
//! (HVS-II's `C7 T C5`), or one run where the anchor is lost; a length tract
//! outside every anchored window (HVS-III's `C6`) is one plain run. Each
//! admitted read's tract decode (SRS-CALL-025) is credited to the run its
//! first call is aligned to. A read gives the length of the run it enters
//! first (its first-run marginal) and, when it chose the window's anchor as
//! its mark and both its marginals are clear, the length of the next run it
//! reads (the total less the first and the mark). A read that shows a run in
//! phase, both its ends resolved by unmasked calls as the consensus reads it
//! (ADR-0074), gives that run's length too. A strand gives a run one clear
//! length when every length its reads give is clear and equal. The reads
//! give a run its length when both strands give it one clear length, the
//! same, or when one strand gives none and the other gives one clear length
//! that a read of it shows in phase; a run whose only views are decodes
//! needs both strands. A run is confirmed when that length is DNA's call, or
//! when the consensus decided it by the two-strand rule (PROP-0012
//! increment 4) and no read's structure differs from DNA's.
//! DNA reports a run's length only when every position it must show to
//! measure the run ([`run_spans`]) is callable; otherwise its
//! decided length is an estimate, never a reading, and the run instead
//! carries the reads' floor reading, the least length any read gives as its
//! top or next. A run of a length tract whose dominant
//! length the reads determine is noted, not checked ([`status`]). Nothing
//! here decides a length; it labels them.

use dna_consensus::ReadRun;
use dna_consensus::length::{Run, entered_run};
use dna_consensus::model::consensus::{ConsensusRun, RunEvidence};
use dna_consensus::runs::{anchored_pair, runs_of};
use dna_kernel::read_evidence::{LengthMarginal, Margin, TractLengths};
use dna_normalize::length_convention::{RULE as CONVENTION_RULE, ShortestDominant};

use crate::model::review_result::{
    LengthConventionResult, LengthReadResult, LengthRunResult, MarginalResult, StrandsResult,
    TractReadResult,
};
use crate::sample_analysis::NotationCall;
pub(super) use crate::sample_analysis::windows::{LengthWindow, windows};
use dna_align::model::alignment::Orientation;

/// DNA's run lengths in a window from the decided `haplotype` over it: the
/// two runs of the repeat base on each side of the first anchor base between
/// two such runs, or the longest repeat run when no anchor splits it. `None`
/// when the haplotype holds an `N` or no repeat run.
pub(super) fn called_runs(haplotype: &str, window: &LengthWindow) -> Option<Vec<usize>> {
    if haplotype.contains('N') {
        return None;
    }
    let runs = runs_of(haplotype.chars());
    let split = window
        .anchor
        .and_then(|(_, anchor_base)| anchored_pair(&runs, window.repeat, anchor_base));
    match split {
        Some(index) => Some(vec![runs[index].1, runs[index + 2].1]),
        None => runs
            .iter()
            .filter(|(base, _)| *base == window.repeat)
            .map(|(_, length)| *length)
            .max()
            .map(|length| vec![length]),
    }
}

/// One read's tract decode, with the reference position its first call is
/// aligned to (on an inserted column, the position the insertion follows on
/// the reference).
#[derive(Clone)]
pub(super) struct Decode<'a> {
    pub(super) read: &'a str,
    pub(super) strand: Orientation,
    pub(super) reference: Option<usize>,
    pub(super) tract: &'a TractLengths,
}

/// One read's runs over a window as the consensus reads them.
pub(super) struct InPhase<'a> {
    pub(super) read: &'a str,
    pub(super) strand: Orientation,
    pub(super) runs: &'a [ReadRun],
}

/// The window's runs a read shows, as `(run, read run)`: the two runs around
/// the first anchor between two repeat runs, else the longest repeat run.
fn window_runs<'r>(runs: &'r [ReadRun], window: &LengthWindow) -> Vec<(usize, &'r ReadRun)> {
    if let Some((_, anchor_base)) = window.anchor
        && let Some(index) = anchored_pair(runs, window.repeat, anchor_base)
    {
        return vec![(1, &runs[index]), (2, &runs[index + 2])];
    }
    runs.iter()
        .filter(|run| run.base == window.repeat)
        .max_by_key(|run| run.length)
        .map(|run| vec![(1, run)])
        .unwrap_or_default()
}

/// The runs of `window` (1 or 2) the consensus decided by the two-strand
/// rule, from the runs of the site that covers it.
pub(super) fn two_strand_runs(runs: &[ConsensusRun], window: &LengthWindow) -> Vec<usize> {
    let Some((_, anchor_base)) = window.anchor else {
        return Vec::new();
    };
    anchored_pair(runs, window.repeat, anchor_base)
        .map(|index| {
            [1, 2]
                .into_iter()
                .filter(|run| runs[index + 2 * (run - 1)].evidence == RunEvidence::TwoStrand)
                .collect()
        })
        .unwrap_or_default()
}

/// The length convention's reading of a run (`policy`, ISFG 2014 3.3.2):
/// each read's reading is its in-phase length, or the policy's reading of
/// its decode; when every read gives one and they hold one adjacent pair
/// `n`, `n + 1`, the run reads `n`. `None` without a policy, and when every
/// read already tops at `n`.
fn convention(reads: &[&Contribution<'_>], policy: Option<ShortestDominant>) -> Option<usize> {
    let policy = policy?;
    let readings = reads
        .iter()
        .map(|contribution| contribution.reading(policy))
        .collect::<Option<Vec<_>>>()?;
    let shortest = *readings.iter().min()?;
    let pair = [shortest, shortest + 1];
    let holds = readings.iter().all(|reading| pair.contains(reading));
    let moves = reads
        .iter()
        .any(|contribution| contribution.result.length != shortest);
    (holds && moves).then_some(shortest)
}

/// The reads' floor reading of a run: the least length any read gives as
/// its top or next.
fn floor(reads: &[&Contribution<'_>]) -> Option<usize> {
    reads
        .iter()
        .flat_map(|contribution| {
            let result = &contribution.result;
            std::iter::once(result.length).chain(result.next)
        })
        .min()
}

const fn margin_label(margin: Margin) -> &'static str {
    match margin {
        Margin::Clear => "clear",
        Margin::NearEqual => "near_equal",
        Margin::Undefined => "undefined",
    }
}

/// Where a read's length of a run comes from.
#[derive(Clone, Copy)]
enum Source<'a> {
    /// The read shows the run in phase.
    InPhase,
    /// The first-run marginal of the read's decode.
    Decode(&'a TractLengths),
    /// The read's total less its first run (the given length) and the mark.
    Derived(&'a TractLengths, usize),
}

/// One read's length for one run of a window.
struct Contribution<'a> {
    run: usize,
    result: LengthReadResult,
    clear: bool,
    /// The read's structure is not DNA's: a mark where DNA calls one run, or
    /// none where it calls two.
    structure_differs: bool,
    source: Source<'a>,
}

impl Contribution<'_> {
    /// The length convention's reading of this read.
    fn reading(&self, policy: ShortestDominant) -> Option<usize> {
        match self.source {
            Source::InPhase => Some(self.result.length),
            Source::Decode(tract) => policy.reading(&tract.first.lengths, &tract.first.shares),
            Source::Derived(tract, first) => policy
                .reading(&tract.total.lengths, &tract.total.shares)?
                .checked_sub(first + 1)
                .filter(|&derived| derived > 0),
        }
    }

    /// The read with its decode's components, for a note.
    fn detail(&self) -> TractReadResult {
        let marginal = |marginal: &LengthMarginal| MarginalResult {
            lengths: marginal.lengths.clone(),
            shares: marginal
                .shares
                .iter()
                .map(|&share| super::review::round(share))
                .collect(),
            z: marginal.z.map(super::review::round),
        };
        let (first, total) = match self.source {
            Source::InPhase => (None, None),
            Source::Decode(tract) | Source::Derived(tract, _) => {
                (Some(marginal(&tract.first)), Some(marginal(&tract.total)))
            }
        };
        TractReadResult {
            read: self.result.read.clone(),
            strand: self.result.strand,
            length: self.result.length,
            source: self.result.source,
            first,
            total,
        }
    }
}

/// The reference span DNA must show callable to report each run it calls
/// in `window` (1-based and inclusive), in run order. Where DNA calls two
/// runs `anchored` at the profile's anchor (its calls in the window change
/// only the runs' lengths, by insertions or deletions of the repeat base),
/// each run is measured between its bounds: run 1 spans the position before
/// the window to the anchor, run 2 the anchor to the position after the
/// window. So a run 1 whose 3' end is only estimated leaves run 2 reported
/// while the anchor stays callable, and a run whose own positions are
/// callable but whose anchor is not stays an estimate. Otherwise (one run, or
/// a structure another call moves) every run spans the whole window. Empty
/// when DNA calls no run there.
pub(super) fn run_spans(
    window: &LengthWindow,
    called: Option<&[usize]>,
    anchored: bool,
) -> Vec<[usize; 2]> {
    let count = called.map_or(0, <[usize]>::len);
    match window.anchor {
        Some((anchor, _)) if anchored && count == 2 => vec![
            [window.first.saturating_sub(1).max(1), anchor],
            [anchor, window.last + 1],
        ],
        _ => vec![[window.first, window.last]; count],
    }
}

/// Whether DNA's `calls` leave `window`'s runs around its anchor: every call
/// in the window inserts or deletes the repeat base, and none lies at the
/// anchor.
pub(super) fn anchored(window: &LengthWindow, calls: &[NotationCall]) -> bool {
    window.anchor.is_some_and(|(anchor, _)| {
        calls
            .iter()
            .filter(|call| (window.first..=window.last).contains(&call.position))
            .all(|call| {
                let repeat = Some(window.repeat);
                call.position != anchor
                    && match (call.reference, call.alternate) {
                        (None, alternate) => alternate == repeat,
                        (reference, None) => reference == repeat,
                        _ => false,
                    }
            })
    })
}

/// One judged run of a window: what `length_unconfirmed` lists, with each
/// strand's one clear length and each read's decode components for a note.
pub(super) struct Judged {
    pub(super) result: LengthRunResult,
    pub(super) strands: StrandsResult,
    pub(super) reads: Vec<TractReadResult>,
}

/// The runs of `window` the `decodes` give a length for, each with DNA's
/// length, the reads' lengths, and how they agree. `called` is DNA's
/// decided runs over the window, which set its structure; a run's length is
/// DNA's only when the run is in `reported` (its span callable,
/// [`run_spans`]), and otherwise the run carries the reads' floor reading
/// instead. A run in
/// `two_strand`, decided by the two-strand rule, is confirmed unless a
/// read's structure differs from DNA's. `policy` is the length convention,
/// when configured. Decodes of another base, or credited to no run of the
/// window, give none.
pub(super) fn judge(
    window: &LengthWindow,
    called: Option<&[usize]>,
    reported: &[usize],
    two_strand: &[usize],
    decodes: &[Decode<'_>],
    in_phase: &[InPhase<'_>],
    policy: Option<ShortestDominant>,
) -> Vec<Judged> {
    let mut contributions = Vec::new();
    for read in in_phase {
        let runs = window_runs(read.runs, window);
        let structure_differs = called.is_some_and(|called| called.len() != runs.len());
        for (run, view) in runs.into_iter().filter(|(_, view)| view.in_phase) {
            contributions.push(Contribution {
                run,
                result: LengthReadResult {
                    read: read.read.to_owned(),
                    strand: read.strand,
                    length: view.length,
                    share: 1.0,
                    next: None,
                    margin: margin_label(Margin::Clear),
                    source: "in_phase",
                },
                clear: true,
                structure_differs,
                source: Source::InPhase,
            });
        }
    }
    for decode in decodes {
        let tract = decode.tract;
        let forward = decode.strand == Orientation::Forward;
        let orient = |base: char| decode.strand.reference_base(base);
        let Some(position) = decode.reference else {
            continue;
        };
        if tract.run != orient(window.repeat) {
            continue;
        }
        let single = called.is_some_and(|runs| runs.len() == 1);
        let Some(entered) = entered_run(
            window.first,
            window.last,
            window.anchor.map(|(anchor, _)| anchor),
            position,
            single,
        )
        .map(Run::number) else {
            continue;
        };
        // A plain run has no anchor: a mark the read chooses lies past it.
        let structure_differs = called.is_some()
            && window
                .anchor
                .is_some_and(|(_, anchor_base)| match tract.mark {
                    Some(mark) => single || mark != orient(anchor_base),
                    None => !single,
                });
        let Some(length) = tract.first.dominant else {
            continue;
        };
        let share = |marginal: &LengthMarginal, length: usize| {
            marginal
                .lengths
                .iter()
                .position(|&kept| kept == length)
                .and_then(|index| marginal.shares.get(index))
                .copied()
                .unwrap_or(0.0)
        };
        contributions.push(Contribution {
            run: entered,
            result: LengthReadResult {
                read: decode.read.to_owned(),
                strand: decode.strand,
                length,
                share: super::review::round(share(&tract.first, length)),
                next: tract.first.next,
                margin: margin_label(tract.first.margin),
                source: "decode",
            },
            clear: tract.first.margin == Margin::Clear,
            structure_differs,
            source: Source::Decode(tract),
        });
        // The next run the read reads, from its total: only with the
        // window's anchor as its mark and both marginals clear.
        let next = if forward { 2 } else { 1 };
        if let Some((_, anchor_base)) = window.anchor
            && !single
            && entered != next
            && tract.mark == Some(orient(anchor_base))
            && tract.first.margin == Margin::Clear
            && tract.total.margin == Margin::Clear
            && let Some(total) = tract.total.dominant
            && let Some(derived) = total.checked_sub(length + 1).filter(|&derived| derived > 0)
        {
            contributions.push(Contribution {
                run: next,
                result: LengthReadResult {
                    read: decode.read.to_owned(),
                    strand: decode.strand,
                    length: derived,
                    share: super::review::round(share(&tract.total, total)),
                    next: None,
                    margin: margin_label(Margin::Clear),
                    source: "derived",
                },
                clear: true,
                structure_differs,
                source: Source::Derived(tract, length),
            });
        }
    }
    let runs = contributions
        .iter()
        .map(|contribution| contribution.run)
        .collect::<std::collections::BTreeSet<_>>();
    runs.into_iter()
        .map(|run| {
            let reads = contributions
                .iter()
                .filter(|contribution| contribution.run == run)
                .collect::<Vec<_>>();
            let called_length = called
                .filter(|_| reported.contains(&run))
                .and_then(|runs| runs.get(run - 1))
                .copied();
            // Each strand's one clear length: every length its reads give
            // clear and equal; `Err` when they are not.
            let strand = |forward: bool| -> Option<Result<usize, ()>> {
                let given = reads
                    .iter()
                    .filter(|contribution| {
                        (contribution.result.strand == Orientation::Forward) == forward
                    })
                    .collect::<Vec<_>>();
                if given.is_empty() {
                    return None;
                }
                let lengths = given
                    .iter()
                    .map(|contribution| contribution.result.length)
                    .collect::<std::collections::BTreeSet<_>>();
                Some(
                    match (
                        given.iter().all(|contribution| contribution.clear),
                        lengths.len(),
                    ) {
                        (true, 1) => Ok(given[0].result.length),
                        _ => Err(()),
                    },
                )
            };
            // A read in phase through the run shows its length directly
            // (ADR-0074): the other strand giving none is no contradiction.
            let in_phase = reads
                .iter()
                .any(|contribution| contribution.result.source == "in_phase");
            let read_length = match (strand(true), strand(false)) {
                _ if reads
                    .iter()
                    .any(|contribution| contribution.structure_differs) =>
                {
                    Err("structure_differs")
                }
                // The two-strand rule decided DNA's length.
                _ if called_length.is_some() && two_strand.contains(&run) => {
                    called_length.ok_or("undecided_call")
                }
                (Some(Ok(forward)), Some(Ok(reverse))) if forward != reverse => Err("disagree"),
                (Some(Err(())), _) | (_, Some(Err(()))) => Err("not_clear"),
                (Some(Ok(length)), None) | (None, Some(Ok(length))) if in_phase => Ok(length),
                (None, _) | (_, None) => Err("one_strand"),
                (Some(Ok(length)), Some(Ok(_))) => Ok(length),
            };
            let agreement = match (read_length, called_length) {
                (Err(agreement), _) => agreement,
                (Ok(_), None) => "undecided_call",
                (Ok(length), Some(called)) if called != length => "differs_from_call",
                (Ok(_), Some(_)) => "confirmed",
            };
            let one = |forward: bool| strand(forward).and_then(Result::ok);
            Judged {
                strands: StrandsResult {
                    forward: one(true),
                    reverse: one(false),
                },
                reads: reads
                    .iter()
                    .map(|contribution| contribution.detail())
                    .collect(),
                result: LengthRunResult {
                    run,
                    base: window.repeat,
                    called: called_length,
                    agreement,
                    convention: (!matches!(agreement, "confirmed" | "structure_differs"))
                        .then(|| convention(&reads, policy))
                        .flatten()
                        .map(|length| LengthConventionResult {
                            length,
                            rule: CONVENTION_RULE,
                        }),
                    floor: called_length.is_none().then(|| floor(&reads)).flatten(),
                    reads: reads
                        .into_iter()
                        .map(|contribution| contribution.result.clone())
                        .collect(),
                },
            }
        })
        .collect()
}

/// How a judged run of a length tract stands, with the length a reviewer
/// would record: its dominant length determined, so the run is noted, not
/// checked; `None` when it stays a check.
///
/// - `confirmed`: the reads confirm DNA's reported length;
/// - `dominant`: every read tops at DNA's reported length;
/// - `convention`: the reads hold one adjacent pair `n`, `n + 1` that the
///   length convention reads as `n`, and DNA's length is one of the pair;
/// - `unconfirmed`: DNA reports no length (the window is not callable), so
///   the run carries the reads' floor reading, unconfirmed.
///
/// `dominant` and `convention` need the run's structure shown: a read in
/// phase, or decodes on both strands. A read structure other than DNA's (its
/// anchor base read otherwise, as with `T310C` or `16189C`) always stays a
/// check, as do decodes on one strand alone and reads that top at lengths
/// no convention resolves, among them an in-phase read of another length.
pub(super) fn status(run: &LengthRunResult) -> Option<(&'static str, Option<usize>)> {
    if run.agreement == "structure_differs" {
        return None;
    }
    let Some(called) = run.called else {
        return Some(("unconfirmed", run.floor));
    };
    if run.agreement == "confirmed" {
        return Some(("confirmed", Some(called)));
    }
    // The structure the decodes rest on must be shown: by a read in phase,
    // or by decodes on both strands. A decode alone on one strand may rest
    // on a spurious or dropped mark, and its margin is over-confident under
    // noise (trace-simulator benchmark, 2026-10-21).
    let strand = |strand: Orientation| run.reads.iter().any(|read| read.strand == strand);
    let shown = run.reads.iter().any(|read| read.source == "in_phase")
        || (strand(Orientation::Forward) && strand(Orientation::Reverse));
    if !shown {
        return None;
    }
    if run.reads.iter().all(|read| read.length == called) {
        return Some(("dominant", Some(called)));
    }
    let length = run.convention.as_ref()?.length;
    let pair = [length, length + 1];
    (pair.contains(&called) && run.reads.iter().all(|read| pair.contains(&read.length)))
        .then_some(("convention", Some(length)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sample_analysis::windows::TractRuns;
    use dna_kernel::read_evidence::{LengthMarginal, StructureFit};

    const HVS2: LengthWindow = LengthWindow {
        first: 303,
        last: 315,
        repeat: 'C',
        anchor: Some((310, 'T')),
        recorded: true,
        tract: Some(TractRuns {
            index: 0,
            runs: [true, false],
            single: false,
        }),
    };

    /// The judged runs as `length_unconfirmed` lists them, under the
    /// shipped length convention (ratio 0.6).
    fn judge(
        window: &LengthWindow,
        called: Option<&[usize]>,
        reported: bool,
        two_strand: &[usize],
        decodes: &[Decode<'_>],
        in_phase: &[InPhase<'_>],
    ) -> Vec<LengthRunResult> {
        // Every run DNA calls is reported, or none.
        let reported = if reported {
            (1..=called.map_or(0, <[usize]>::len)).collect::<Vec<_>>()
        } else {
            Vec::new()
        };
        super::judge(
            window,
            called,
            &reported,
            two_strand,
            decodes,
            in_phase,
            ShortestDominant::new(0.6),
        )
        .into_iter()
        .map(|judged| judged.result)
        .collect()
    }

    fn marginal(lengths: &[usize], dominant: usize, margin: Margin) -> LengthMarginal {
        LengthMarginal {
            lengths: lengths.to_vec(),
            shares: lengths
                .iter()
                .map(|&length| if length == dominant { 1.0 } else { 0.0 })
                .collect(),
            errors: vec![0.01; lengths.len()],
            dominant: Some(dominant),
            next: None,
            z: None,
            z_adjusted: None,
            margin,
        }
    }

    fn tract(
        run: char,
        mark: Option<char>,
        first: usize,
        total: usize,
        margin: Margin,
    ) -> TractLengths {
        TractLengths {
            onset: 20,
            call_start: 10,
            extent: total,
            window: total + 5,
            run,
            mark,
            candidates: vec![StructureFit {
                mark,
                residual: 0.01,
            }],
            residual: 0.01,
            sigma: 0.01,
            autocorrelation: 0.0,
            first: marginal(&[first], first, margin),
            total: marginal(&[total], total, Margin::Clear),
            cleaned: None,
        }
    }

    #[test]
    fn reads_dna_s_run_lengths_around_the_anchor() {
        assert_eq!(called_runs("CCCCCCCTCCCCC", &HVS2), Some(vec![7, 5]));
        assert_eq!(called_runs("ACCCCCCCCTCCCCCG", &HVS2), Some(vec![8, 5]));
        assert_eq!(called_runs("CCCCCCCCCCCCC", &HVS2), Some(vec![13]));
        assert_eq!(called_runs("CCCCCCCNCCCCC", &HVS2), None);
    }

    // SRS-API-020
    #[test]
    fn confirms_a_run_both_strands_read_as_called() {
        let forward = tract('C', Some('T'), 8, 14, Margin::Clear);
        let reverse = tract('G', Some('A'), 5, 14, Margin::Clear);
        let decodes = [
            Decode {
                read: "F",
                strand: Orientation::Forward,
                reference: Some(303),
                tract: &forward,
            },
            Decode {
                read: "R",
                strand: Orientation::Reverse,
                reference: Some(315),
                tract: &reverse,
            },
        ];
        let runs = judge(&HVS2, Some(&[8, 5]), true, &[], &decodes, &[]);
        assert_eq!(runs.len(), 2);
        // Run 1: F's first run and R's derived one agree with the call.
        assert_eq!(runs[0].agreement, "confirmed");
        assert_eq!(runs[0].reads.len(), 2);
        assert_eq!(runs[0].reads[1].source, "derived");
        assert_eq!(runs[1].agreement, "confirmed");
        // The call one longer than both strands read it.
        let runs = judge(&HVS2, Some(&[9, 5]), true, &[], &decodes, &[]);
        assert_eq!(runs[0].agreement, "differs_from_call");
    }

    #[test]
    fn flags_disagreement_near_equal_and_one_strand() {
        let forward = tract('C', Some('T'), 8, 14, Margin::Clear);
        let reverse = tract('G', Some('A'), 5, 15, Margin::Clear);
        let near = tract('G', Some('A'), 5, 14, Margin::NearEqual);
        let decode = |read, strand, reference, tract| Decode {
            read,
            strand,
            reference: Some(reference),
            tract,
        };
        // R's derived run 1 is 9 where F reads 8.
        let runs = judge(
            &HVS2,
            Some(&[8, 5]),
            true,
            &[],
            &[
                decode("F", Orientation::Forward, 303, &forward),
                decode("R", Orientation::Reverse, 315, &reverse),
            ],
            &[],
        );
        assert_eq!(runs[0].agreement, "disagree");
        let runs = judge(
            &HVS2,
            Some(&[8, 5]),
            true,
            &[],
            &[
                decode("F", Orientation::Forward, 303, &forward),
                decode("R", Orientation::Reverse, 315, &near),
            ],
            &[],
        );
        assert_eq!(runs[0].agreement, "one_strand");
        assert_eq!(runs[1].agreement, "not_clear");
        // A run the consensus decided by the two-strand rule is confirmed
        // whatever the decodes' own agreement.
        let runs = judge(
            &HVS2,
            Some(&[8, 5]),
            true,
            &[1, 2],
            &[
                decode("F", Orientation::Forward, 303, &forward),
                decode("R", Orientation::Reverse, 315, &near),
            ],
            &[],
        );
        assert_eq!(runs[0].agreement, "confirmed");
        assert_eq!(runs[1].agreement, "confirmed");
        // A single run where DNA calls two is a structure that differs.
        let single = tract('C', None, 13, 13, Margin::Clear);
        let runs = judge(
            &HVS2,
            Some(&[8, 5]),
            true,
            &[],
            &[decode("F", Orientation::Forward, 303, &single)],
            &[],
        );
        assert_eq!(runs[0].agreement, "structure_differs");
        // Even then, a read whose structure differs keeps the run flagged.
        let runs = judge(
            &HVS2,
            Some(&[8, 5]),
            true,
            &[1],
            &[decode("F", Orientation::Forward, 303, &single)],
            &[],
        );
        assert_eq!(runs[0].agreement, "structure_differs");
        // Another base, or a start outside the window, gives nothing.
        let other = tract('A', None, 6, 6, Margin::Clear);
        assert!(
            judge(
                &HVS2,
                Some(&[8, 5]),
                true,
                &[],
                &[
                    decode("F", Orientation::Forward, 303, &other),
                    decode("F", Orientation::Forward, 290, &forward),
                ],
                &[],
            )
            .is_empty()
        );
    }

    // SRS-API-020
    #[test]
    fn suggests_the_shortest_dominant_of_a_near_equal_pair() {
        let mut forward = tract('C', Some('T'), 8, 14, Margin::Clear);
        forward.first = LengthMarginal {
            lengths: vec![8, 9],
            shares: vec![0.55, 0.45],
            errors: vec![0.05, 0.05],
            dominant: Some(8),
            next: Some(9),
            z: Some(0.7),
            z_adjusted: Some(0.5),
            margin: Margin::NearEqual,
        };
        let decodes = [Decode {
            read: "F",
            strand: Orientation::Forward,
            reference: Some(303),
            tract: &forward,
        }];
        let through = |length| {
            [
                ReadRun {
                    base: 'C',
                    length,
                    in_phase: true,
                },
                ReadRun {
                    base: 'T',
                    length: 1,
                    in_phase: true,
                },
                ReadRun {
                    base: 'C',
                    length: 5,
                    in_phase: true,
                },
            ]
        };
        let judged = |length| {
            let runs = through(length);
            judge(
                &HVS2,
                Some(&[9, 5]),
                true,
                &[],
                &decodes,
                &[InPhase {
                    read: "R",
                    strand: Orientation::Reverse,
                    runs: &runs,
                }],
            )
        };
        // Both strands read C8/C9: the convention reads the shortest
        // dominant, as a suggestion beside DNA's own call, which stays
        // unconfirmed.
        let runs = judged(9);
        assert_eq!(runs[0].agreement, "not_clear");
        assert_eq!(
            runs[0].convention,
            Some(LengthConventionResult {
                length: 8,
                rule: "shortest_dominant",
            })
        );
        assert_eq!(runs[0].called, Some(9));
        // A read outside the pair leaves no suggestion.
        assert_eq!(judged(10)[0].convention, None);
    }

    // SRS-API-020
    #[test]
    fn notes_a_tract_run_whose_dominant_length_the_reads_determine() {
        let read = |strand, length, source| LengthReadResult {
            read: "R".to_owned(),
            strand,
            length,
            share: 0.6,
            next: None,
            margin: "near_equal",
            source,
        };
        let run = |called, agreement, reads, convention: Option<usize>| LengthRunResult {
            run: 1,
            base: 'C',
            called,
            agreement,
            convention: convention.map(|length| LengthConventionResult {
                length,
                rule: CONVENTION_RULE,
            }),
            floor: called.is_none().then_some(7),
            reads,
        };
        let (forward, reverse) = (Orientation::Forward, Orientation::Reverse);
        // Every read tops at DNA's C8, on both strands: the dominant.
        let both = vec![read(forward, 8, "decode"), read(reverse, 8, "derived")];
        assert_eq!(
            status(&run(Some(8), "not_clear", both.clone(), None)),
            Some(("dominant", Some(8)))
        );
        // A decode alone on one strand rests on an unshown structure.
        let one = vec![read(forward, 8, "decode")];
        assert_eq!(status(&run(Some(8), "one_strand", one, None)), None);
        // A read in phase shows it.
        let shown = vec![read(forward, 8, "decode"), read(forward, 8, "in_phase")];
        assert_eq!(
            status(&run(Some(8), "one_strand", shown, None)),
            Some(("dominant", Some(8)))
        );
        // A near-equal pair C8/C9 that the convention reads as C8, DNA C9.
        let pair = vec![read(forward, 9, "decode"), read(reverse, 8, "in_phase")];
        assert_eq!(
            status(&run(Some(9), "disagree", pair.clone(), Some(8))),
            Some(("convention", Some(8)))
        );
        // Without a convention reading, the dominant is undecided.
        assert_eq!(status(&run(Some(9), "disagree", pair, None)), None);
        // Reads beyond the pair, or another structure, stay a check.
        let wide = vec![read(forward, 10, "decode"), read(reverse, 8, "decode")];
        assert_eq!(status(&run(Some(9), "disagree", wide, Some(8))), None);
        assert_eq!(status(&run(Some(8), "structure_differs", both, None)), None);
        // DNA reports no length: the floor, unconfirmed.
        assert_eq!(
            status(&run(None, "undecided_call", Vec::new(), None)),
            Some(("unconfirmed", Some(7)))
        );
    }

    // SRS-API-020
    #[test]
    fn judges_a_plain_tract_as_one_run_whatever_mark_a_read_chooses() {
        let plain = LengthWindow {
            first: 568,
            last: 573,
            repeat: 'C',
            anchor: None,
            recorded: true,
            tract: Some(TractRuns {
                index: 0,
                runs: [true, false],
                single: true,
            }),
        };
        // A forward decode that chose the A after the run as its mark, and
        // a reverse one entering from 573.
        let forward = tract('C', Some('A'), 7, 9, Margin::Clear);
        let reverse = tract('G', None, 7, 7, Margin::Clear);
        let decodes = [
            Decode {
                read: "F",
                strand: Orientation::Forward,
                reference: Some(568),
                tract: &forward,
            },
            Decode {
                read: "R",
                strand: Orientation::Reverse,
                reference: Some(573),
                tract: &reverse,
            },
        ];
        let runs = judge(&plain, Some(&[7]), true, &[], &decodes, &[]);
        assert_eq!(runs.len(), 1);
        assert_eq!((runs[0].run, runs[0].agreement), (1, "confirmed"));
        assert!(plain.in_tract(1, true) && !HVS2.in_tract(2, false) && HVS2.in_tract(1, false));
    }

    // SRS-API-020
    #[test]
    fn gives_the_reads_floor_where_dna_reports_no_length() {
        let mut forward = tract('C', Some('T'), 9, 15, Margin::Clear);
        forward.first = LengthMarginal {
            lengths: vec![9, 10],
            shares: vec![0.55, 0.45],
            errors: vec![0.05, 0.05],
            dominant: Some(9),
            next: Some(10),
            z: Some(0.7),
            z_adjusted: Some(0.5),
            margin: Margin::NearEqual,
        };
        let reverse = tract('G', Some('A'), 5, 14, Margin::Clear);
        let decodes = [
            Decode {
                read: "F",
                strand: Orientation::Forward,
                reference: Some(303),
                tract: &forward,
            },
            Decode {
                read: "R",
                strand: Orientation::Reverse,
                reference: Some(315),
                tract: &reverse,
            },
        ];
        // F reads run 2 in phase after a dephased run 1.
        let through = [
            ReadRun {
                base: 'C',
                length: 10,
                in_phase: false,
            },
            ReadRun {
                base: 'T',
                length: 1,
                in_phase: true,
            },
            ReadRun {
                base: 'C',
                length: 5,
                in_phase: true,
            },
        ];
        let in_phase = [InPhase {
            read: "F",
            strand: Orientation::Forward,
            runs: &through,
        }];
        // DNA's estimate, C10 then C5, is not reported: no run carries it,
        // and each carries the least length a read gives instead.
        let runs = judge(&HVS2, Some(&[10, 5]), false, &[], &decodes, &in_phase);
        assert_eq!(runs.len(), 2);
        assert_eq!(runs[0].called, None);
        assert_eq!(runs[0].floor, Some(8));
        assert_eq!(runs[0].agreement, "not_clear");
        // Both strands give run 2 one clear length; DNA reports none.
        assert_eq!(runs[1].called, None);
        assert_eq!(runs[1].floor, Some(5));
        assert_eq!(runs[1].agreement, "undecided_call");
        // Reported, the same runs carry DNA's lengths and no floor.
        let runs = judge(&HVS2, Some(&[10, 5]), true, &[], &decodes, &in_phase);
        assert_eq!(runs[0].called, Some(10));
        assert_eq!(runs[0].floor, None);
        assert_eq!(runs[1].agreement, "confirmed");
        // Run 2 alone reported (run 1's 3' end an estimate, the anchor
        // callable): run 2 is confirmed, run 1 keeps the floor reading.
        let runs = super::judge(
            &HVS2,
            Some(&[10, 5]),
            &[2],
            &[],
            &decodes,
            &in_phase,
            ShortestDominant::new(0.6),
        );
        assert_eq!(runs[0].result.called, None);
        assert_eq!(runs[0].result.floor, Some(8));
        assert_eq!(runs[1].result.called, Some(5));
        assert_eq!(runs[1].result.agreement, "confirmed");
    }

    // SRS-API-020
    #[test]
    fn bounds_each_run_by_its_flank_and_the_anchor() {
        // Anchored: run 1 from 302 to the anchor 310, run 2 from 310 to 316.
        assert_eq!(
            run_spans(&HVS2, Some(&[8, 5]), true),
            vec![[302, 310], [310, 316]]
        );
        // A structure another call moves: each run spans the window.
        assert_eq!(
            run_spans(&HVS2, Some(&[8, 5]), false),
            vec![[303, 315], [303, 315]]
        );
        // One run where the anchor is lost: the window.
        assert_eq!(run_spans(&HVS2, Some(&[14]), true), vec![[303, 315]]);
        // No run without a call.
        assert!(run_spans(&HVS2, None, true).is_empty());
    }

    // SRS-API-020
    #[test]
    fn runs_lie_around_the_anchor_when_calls_change_only_their_lengths() {
        let call =
            |position: usize, ordinal: usize, reference: Option<char>, alternate: Option<char>| {
                NotationCall {
                    position,
                    ordinal,
                    reference,
                    alternate,
                    text: String::new(),
                    convention: None,
                }
            };
        // 309.1C, 309.2C, 315.1C, and a call outside the window.
        let lengths = [
            call(309, 1, None, Some('C')),
            call(309, 2, None, Some('C')),
            call(315, 1, None, Some('C')),
            call(263, 0, Some('A'), Some('G')),
        ];
        assert!(anchored(&HVS2, &lengths));
        // A deleted C of run 2.
        assert!(anchored(&HVS2, &[call(315, 0, Some('C'), None)]));
        // T310C changes the anchor; 305T splits run 1; an inserted T moves
        // the anchor.
        assert!(!anchored(&HVS2, &[call(310, 0, Some('T'), Some('C'))]));
        assert!(!anchored(&HVS2, &[call(305, 0, Some('C'), Some('T'))]));
        assert!(!anchored(&HVS2, &[call(309, 1, None, Some('T'))]));
    }

    // SRS-API-020
    #[test]
    fn counts_a_read_in_phase_as_its_strand_s_length() {
        let run = |base, length, in_phase| ReadRun {
            base,
            length,
            in_phase,
        };
        let forward = tract('C', Some('T'), 8, 14, Margin::Clear);
        let decodes = [Decode {
            read: "F",
            strand: Orientation::Forward,
            reference: Some(303),
            tract: &forward,
        }];
        let through = [
            run('A', 1, true),
            run('C', 8, true),
            run('T', 1, true),
            run('C', 5, true),
        ];
        let reverse = |runs| InPhase {
            read: "R",
            strand: Orientation::Reverse,
            runs,
        };
        let runs = judge(
            &HVS2,
            Some(&[8, 5]),
            true,
            &[],
            &decodes,
            &[reverse(&through)],
        );
        assert_eq!(
            runs.iter().map(|run| run.agreement).collect::<Vec<_>>(),
            ["confirmed", "confirmed"]
        );
        assert!(runs[0].reads.iter().any(|read| read.source == "in_phase"));
        // A strand whose reads give two lengths gives none.
        let near = tract('G', Some('A'), 6, 14, Margin::Clear);
        let mixed = [
            Decode {
                read: "R2",
                strand: Orientation::Reverse,
                reference: Some(315),
                tract: &near,
            },
            decodes[0].clone(),
        ];
        let runs = judge(
            &HVS2,
            Some(&[8, 5]),
            true,
            &[],
            &mixed,
            &[reverse(&through)],
        );
        assert_eq!(runs[1].agreement, "not_clear");
        // A read in phase through one run where DNA calls two differs.
        let one = [run('C', 13, true)];
        let runs = judge(&HVS2, Some(&[8, 5]), true, &[], &decodes, &[reverse(&one)]);
        assert_eq!(runs[0].agreement, "structure_differs");
    }

    // SRS-API-020
    #[test]
    fn one_strand_confirms_only_with_a_read_in_phase() {
        let run = |base, length, in_phase| ReadRun {
            base,
            length,
            in_phase,
        };
        // F reads run 1 as 8 and, from its total, run 2 as 5; no reverse
        // read gives either run a length.
        let forward = tract('C', Some('T'), 8, 14, Margin::Clear);
        let decodes = [Decode {
            read: "F",
            strand: Orientation::Forward,
            reference: Some(303),
            tract: &forward,
        }];
        let agreements = |called: &[usize], in_phase: &[InPhase<'_>]| {
            judge(&HVS2, Some(called), true, &[], &decodes, in_phase)
                .iter()
                .map(|run| run.agreement)
                .collect::<Vec<_>>()
        };
        // A decode alone on one strand does not confirm.
        assert_eq!(agreements(&[8, 5], &[]), ["one_strand", "one_strand"]);
        // A forward read in phase through both runs confirms them.
        let through = [
            run('A', 1, true),
            run('C', 8, true),
            run('T', 1, true),
            run('C', 5, true),
        ];
        let read = |runs| InPhase {
            read: "F2",
            strand: Orientation::Forward,
            runs,
        };
        assert_eq!(
            agreements(&[8, 5], &[read(&through)]),
            ["confirmed", "confirmed"]
        );
        // The same length, not DNA's, differs from the call.
        assert_eq!(
            agreements(&[9, 5], &[read(&through)]),
            ["differs_from_call", "confirmed"]
        );
        // In phase through run 2 only: run 1 still rests on the decode.
        let second = [
            run('A', 1, true),
            run('C', 8, false),
            run('T', 1, true),
            run('C', 5, true),
        ];
        assert_eq!(
            agreements(&[8, 5], &[read(&second)]),
            ["one_strand", "confirmed"]
        );
        // A read in phase with another length than the decode's: the strand
        // gives no one length.
        let longer = [
            run('A', 1, true),
            run('C', 9, true),
            run('T', 1, true),
            run('C', 5, true),
        ];
        assert_eq!(
            agreements(&[8, 5], &[read(&longer)]),
            ["not_clear", "confirmed"]
        );
    }
}
