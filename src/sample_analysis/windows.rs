//! The profile's judged repeat windows (anchored windows and the length
//! tracts outside them), as the review judges their run lengths, and the
//! anchored ones a run context gathers its length noise over.

use dna_call::model::called_read::CalledRead;
use dna_consensus::noise::CrossStrand;
use dna_kernel::model::reference::Reference;
use dna_kernel::profile::{ConventionKind, ConventionLabel, Profile, WindowRule};

use crate::config::Config;

/// One repeat window of the profile whose run lengths the review judges,
/// 1-based and inclusive: an anchored window, or a length tract outside
/// every anchored window.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct LengthWindow {
    pub(crate) first: usize,
    pub(crate) last: usize,
    pub(crate) repeat: char,
    /// The anchor's reference position and base; `None` for a plain run.
    pub(crate) anchor: Option<(usize, char)>,
    /// Whether the profile records the window's run lengths: a run-length
    /// convention of its repeat base over it leaves them unrecorded, so the
    /// decodes are context there, never a check.
    pub(crate) recorded: bool,
    /// The profile's length tract over the window, and which runs it holds.
    pub(crate) tract: Option<TractRuns>,
}

/// Which runs of a window lie in its length tract.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct TractRuns {
    /// Index of the tract in the profile's `length_tracts`.
    pub(crate) index: usize,
    /// Whether run 1 and run 2 lie in the tract.
    pub(crate) runs: [bool; 2],
    /// Whether the one run left where the anchor is lost lies in it.
    pub(crate) single: bool,
}

impl LengthWindow {
    /// Whether run `run` lies in the window's length tract, when DNA calls
    /// one run (`single`) or two.
    pub(crate) fn in_tract(&self, run: usize, single: bool) -> bool {
        self.tract.is_some_and(|tract| {
            if single {
                run == 1 && tract.single
            } else {
                run.checked_sub(1)
                    .and_then(|index| tract.runs.get(index))
                    .copied()
                    .unwrap_or(false)
            }
        })
    }
}

/// The profile's judged windows, in reference order: its anchored repeat
/// windows (whose rules split two runs of one base by an anchor base), each
/// with whether the profile records its run lengths and the length tract it
/// holds, and each length tract no anchored window holds, as one plain run.
pub(crate) fn windows(profile: &Profile) -> Vec<LengthWindow> {
    let inside =
        |span: [usize; 2], window: [usize; 2]| window[0] <= span[0] && span[1] <= window[1];
    let mut windows = profile
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
            // Each run's reference span: the repeat bases on either side of
            // the anchor in the window's reference sequence.
            let bases = window.sequence.as_bytes();
            let mut start = anchor.index;
            while start > 0 && bases.get(start - 1) == Some(&anchor.repeat) {
                start -= 1;
            }
            let mut end = anchor.index + 1;
            while bases.get(end) == Some(&anchor.repeat) {
                end += 1;
            }
            let run1 = [first + start, first + anchor.index - 1];
            let run2 = [first + anchor.index + 1, first + end - 1];
            let tract = profile
                .length_tracts
                .iter()
                .position(|tract| {
                    tract.base == anchor.repeat
                        && tract.window[0] <= last
                        && first <= tract.window[1]
                })
                .map(|index| {
                    let span = profile.length_tracts[index].window;
                    TractRuns {
                        index,
                        runs: [inside(run1, span), inside(run2, span)],
                        single: inside([run1[0], run2[1]], span),
                    }
                });
            Some(LengthWindow {
                first,
                last,
                repeat: char::from(anchor.repeat),
                anchor: Some((first + anchor.index, char::from(anchor.base))),
                recorded,
                tract,
            })
        })
        .collect::<Vec<_>>();
    for (index, tract) in profile.length_tracts.iter().enumerate() {
        let [first, last] = tract.window;
        if windows
            .iter()
            .any(|window| window.first <= last && first <= window.last)
        {
            continue;
        }
        windows.push(LengthWindow {
            first,
            last,
            repeat: char::from(tract.base),
            anchor: None,
            recorded: true,
            tract: Some(TractRuns {
                index,
                runs: [true, false],
                single: true,
            }),
        });
    }
    windows.sort_by_key(|window| window.first);
    windows
}

/// The profile's recorded anchored windows, each `[first, last]`: those
/// whose cross-strand pairs a run context's length section holds.
pub(crate) fn recorded_windows(profile: &Profile) -> Vec<[usize; 2]> {
    windows(profile)
        .iter()
        .filter(|window| window.recorded && window.anchor.is_some())
        .map(|window| [window.first, window.last])
        .collect()
}

/// The cross-strand pairs one sample's admitted `reads`, in its evidence's
/// order, give the two-strand rule's length noise over the profile's
/// recorded anchored windows (PROP-0014, SRS-CONS-014): per window run, the
/// forward and reverse decodes' difference in the rule's share and its
/// variance, where each strand decodes the run once with the same top.
pub(crate) fn sample_pairs(
    profile: &Profile,
    reference: &Reference,
    config: &Config,
    reads: &[&CalledRead],
) -> Vec<CrossStrand> {
    let mut pairs = Vec::new();
    for window in windows(profile).iter().filter(|window| window.recorded) {
        let Some((anchor, anchor_base)) = window.anchor else {
            continue;
        };
        let rule_window = dna_consensus::length::Window {
            first: window.first,
            last: window.last,
            anchor,
            repeat: window.repeat,
            anchor_base,
        };
        let shown = dna_consensus::read_runs_over(
            reads,
            reference,
            &config.core,
            &config.consensus,
            window.first - 1,
            window.last - 1,
        );
        let views = dna_consensus::length::views(&rule_window, reads, &shown);
        pairs.extend(
            views
                .iter()
                .enumerate()
                .filter_map(|(index, views)| dna_consensus::noise::cross_strand(views, index + 1)),
        );
    }
    pairs
}
