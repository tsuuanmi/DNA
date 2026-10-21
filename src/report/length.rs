//! Run lengths in the profile's anchored repeat windows, judged from the
//! reads' tract decodes (SRS-API-020, PROP-0012 increment 3).
//!
//! An anchored window holds two runs of one base around an anchor base
//! (HVS-II's `C7 T C5`), or one run where the anchor is lost. Each admitted
//! read's tract decode (SRS-CALL-025) is credited to the run its first call
//! is aligned to. A read gives the length of the run it enters first (its
//! first-run marginal) and, when it chose the window's anchor as its mark
//! and both its marginals are clear, the length of the next run it reads
//! (the total less the first and the mark). A read that shows a run in
//! phase, both its ends resolved by unmasked calls as the consensus reads it
//! (ADR-0074), gives that run's length too. A strand gives a run one clear
//! length when every length its reads give is clear and equal. The reads
//! give a run its length when both strands give it one clear length, the
//! same, or when one strand gives none and the other gives one clear length
//! that a read of it shows in phase; a run whose only views are decodes
//! needs both strands. A run is confirmed when that length is DNA's call, or
//! when the consensus decided it by the two-strand rule (PROP-0012
//! increment 4) and no read's structure differs from DNA's.
//! DNA reports a window's lengths only when every position of it is
//! callable; otherwise its decided lengths are estimates, never a reading,
//! and a run instead carries the reads' floor reading, the least length any
//! read gives as its top or next. Nothing here decides a length; it labels
//! them.

use dna_consensus::ReadRun;
use dna_consensus::model::consensus::{ConsensusRun, RunEvidence};
use dna_kernel::profile::{ConventionKind, ConventionLabel, Profile, WindowRule};
use dna_kernel::read_evidence::{Margin, TractLengths};

use crate::model::review_result::{LengthConventionResult, LengthReadResult, LengthRunResult};
use dna_align::model::alignment::Orientation;

/// One anchored repeat window of the profile, 1-based and inclusive.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct LengthWindow {
    pub(super) first: usize,
    pub(super) last: usize,
    /// The anchor's reference position.
    pub(super) anchor: usize,
    pub(super) repeat: char,
    pub(super) anchor_base: char,
    /// Whether the profile records the window's run lengths: a run-length
    /// convention of its repeat base over it leaves them unrecorded, so the
    /// decodes are context there, never a check.
    pub(super) recorded: bool,
}

/// The profile's anchored repeat windows: those whose rules split two runs
/// of one base by an anchor base, each with whether the profile records its
/// run lengths.
pub(super) fn windows(profile: &Profile) -> Vec<LengthWindow> {
    profile
        .windows
        .iter()
        .filter_map(|window| {
            let anchor = window.rules.iter().find_map(|rule| match rule {
                WindowRule::AnchoredRunLengths(anchor)
                | WindowRule::AnchorDuplication(anchor)
                | WindowRule::AnchorDeletion(anchor)
                | WindowRule::TerminalRepeatInsertion(anchor)
                | WindowRule::AnchoredTract { anchor, .. } => Some(*anchor),
                WindowRule::CanonicalHaplotype(_) | WindowRule::MotifShift { .. } => None,
            })?;
            let first = window.start + 1;
            let last = window.end();
            let recorded = !profile.conventions.iter().any(|convention| {
                let ConventionKind::RunLength { base, .. } = convention.kind;
                convention.label == ConventionLabel::NotRecorded
                    && base == anchor.repeat
                    && convention.window[0] <= last
                    && first <= convention.window[1]
            });
            Some(LengthWindow {
                first,
                last,
                anchor: first + anchor.index,
                repeat: char::from(anchor.repeat),
                anchor_base: char::from(anchor.base),
                recorded,
            })
        })
        .collect()
}

/// DNA's run lengths in a window from the decided `haplotype` over it: the
/// two runs of the repeat base on each side of the first anchor base between
/// two such runs, or the longest repeat run when no anchor splits it. `None`
/// when the haplotype holds an `N` or no repeat run.
pub(super) fn called_runs(haplotype: &str, window: &LengthWindow) -> Option<Vec<usize>> {
    if haplotype.contains('N') {
        return None;
    }
    let runs = super::review::runs_of(haplotype);
    let split = runs.windows(3).find(|triple| {
        triple[0].0 == window.repeat
            && triple[1] == (window.anchor_base, 1)
            && triple[2].0 == window.repeat
    });
    match split {
        Some(triple) => Some(vec![triple[0].1, triple[2].1]),
        None => runs
            .iter()
            .filter(|(base, _)| *base == window.repeat)
            .map(|(_, length)| *length)
            .max()
            .map(|length| vec![length]),
    }
}

/// One read's tract decode, with the reference position its first call is
/// aligned to (the nearest aligned call before it in the read's order when
/// it sits on an inserted column).
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
    if let Some(index) = (0..runs.len().saturating_sub(2)).find(|&index| {
        runs[index].base == window.repeat
            && runs[index + 1].base == window.anchor_base
            && runs[index + 1].length == 1
            && runs[index + 2].base == window.repeat
    }) {
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
    (0..runs.len().saturating_sub(2))
        .find(|&index| {
            runs[index].base == window.repeat
                && runs[index + 1].base == window.anchor_base
                && runs[index + 1].length == 1
                && runs[index + 2].base == window.repeat
        })
        .map(|index| {
            [1, 2]
                .into_iter()
                .filter(|run| runs[index + 2 * (run - 1)].evidence == RunEvidence::TwoStrand)
                .collect()
        })
        .unwrap_or_default()
}

/// The shorter length of the adjacent pair `n`, `n + 1` both strands read a
/// run as, when at least one read's top two are near-equal: every read with
/// a next length holds the pair as its top two, and every other read gives
/// one of the two. The lab reads such a run as its shorter member (the run
/// ends at the first mixed call of the reverse read); DNA only shows it.
fn near_equal_pair(reads: &[&Contribution]) -> Option<usize> {
    let shorter = reads.iter().find_map(|contribution| {
        let result = &contribution.result;
        (result.margin == margin_label(Margin::NearEqual))
            .then(|| result.next.map(|next| result.length.min(next)))
            .flatten()
    })?;
    let pair = [shorter, shorter + 1];
    let holds = reads.iter().all(|contribution| {
        let result = &contribution.result;
        pair.contains(&result.length)
            && result
                .next
                .is_none_or(|next| pair.contains(&next) && next != result.length)
    });
    let strands = |forward: bool| {
        reads
            .iter()
            .any(|contribution| (contribution.result.strand == Orientation::Forward) == forward)
    };
    (holds && strands(true) && strands(false)).then_some(shorter)
}

/// The reads' floor reading of a run: the least length any read gives as
/// its top or next.
fn floor(reads: &[&Contribution]) -> Option<usize> {
    reads
        .iter()
        .flat_map(|contribution| {
            let result = &contribution.result;
            std::iter::once(result.length).chain(result.next)
        })
        .min()
}

const fn complement(base: char) -> char {
    match base {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        other => other,
    }
}

const fn margin_label(margin: Margin) -> &'static str {
    match margin {
        Margin::Clear => "clear",
        Margin::NearEqual => "near_equal",
        Margin::Undefined => "undefined",
    }
}

/// One read's length for one run of a window.
struct Contribution {
    run: usize,
    result: LengthReadResult,
    clear: bool,
    /// The read's structure is not DNA's: a mark where DNA calls one run, or
    /// none where it calls two.
    structure_differs: bool,
}

/// The runs of `window` the `decodes` give a length for, each with DNA's
/// length, the reads' lengths, and how they agree. `called` is DNA's
/// decided runs over the window, which set its structure; their lengths are
/// DNA's only when it `reported` them (every window position callable), and
/// otherwise each run carries the reads' floor reading instead. A run in
/// `two_strand`, decided by the two-strand rule, is confirmed unless a
/// read's structure differs from DNA's. Decodes of another base, or
/// credited to no run of the window, give none.
pub(super) fn judge(
    window: &LengthWindow,
    called: Option<&[usize]>,
    reported: bool,
    two_strand: &[usize],
    decodes: &[Decode<'_>],
    in_phase: &[InPhase<'_>],
) -> Vec<LengthRunResult> {
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
            });
        }
    }
    for decode in decodes {
        let tract = decode.tract;
        let forward = decode.strand == Orientation::Forward;
        let (run_base, anchor_base) = if forward {
            (window.repeat, window.anchor_base)
        } else {
            (complement(window.repeat), complement(window.anchor_base))
        };
        let Some(position) = decode.reference else {
            continue;
        };
        if tract.run != run_base || position < window.first || position > window.last {
            continue;
        }
        if position == window.anchor {
            continue;
        }
        let single = called.is_some_and(|runs| runs.len() == 1);
        let entered = if single || position < window.anchor {
            1
        } else {
            2
        };
        let structure_differs = called.is_some()
            && match tract.mark {
                Some(mark) => single || mark != anchor_base,
                None => !single,
            };
        let Some(length) = tract.first.dominant else {
            continue;
        };
        let share = |marginal: &dna_kernel::read_evidence::LengthMarginal, length: usize| {
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
        });
        // The next run the read reads, from its total: only with the
        // window's anchor as its mark and both marginals clear.
        let next = if forward { 2 } else { 1 };
        if !single
            && entered != next
            && tract.mark == Some(anchor_base)
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
                .filter(|_| reported)
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
            LengthRunResult {
                run,
                base: window.repeat,
                called: called_length,
                agreement,
                convention: (!matches!(agreement, "confirmed" | "structure_differs"))
                    .then(|| near_equal_pair(&reads))
                    .flatten()
                    .map(|shorter| LengthConventionResult {
                        length: shorter,
                        rule: "shorter_of_near_equal_pair",
                    }),
                floor: called_length.is_none().then(|| floor(&reads)).flatten(),
                reads: reads
                    .into_iter()
                    .map(|contribution| contribution.result.clone())
                    .collect(),
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use dna_kernel::read_evidence::{LengthMarginal, StructureFit};

    const HVS2: LengthWindow = LengthWindow {
        first: 303,
        last: 315,
        anchor: 310,
        repeat: 'C',
        anchor_base: 'T',
        recorded: true,
    };

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
    fn suggests_the_shorter_member_of_a_near_equal_pair() {
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
        // Both strands read C8/C9: the convention reads the shorter, as a
        // suggestion beside DNA's own call, which stays unconfirmed.
        let runs = judged(9);
        assert_eq!(runs[0].agreement, "not_clear");
        assert_eq!(
            runs[0].convention,
            Some(LengthConventionResult {
                length: 8,
                rule: "shorter_of_near_equal_pair",
            })
        );
        assert_eq!(runs[0].called, Some(9));
        // A read outside the pair leaves no suggestion.
        assert_eq!(judged(10)[0].convention, None);
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
