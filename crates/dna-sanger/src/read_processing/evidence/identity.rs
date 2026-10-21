//! The calls a read's own signal vouches for by their base: the leading
//! calls (SRS-CALL-014), the edge calls (SRS-CALL-017), and the start calls
//! (SRS-CALL-018).

use super::{DEPHASED_SIGNAL, POST_HOMOPOLYMER, placement};
use crate::model::basecalls::BaseCalls;
use crate::model::callability::ReadCallability;
use crate::model::locus_evidence::{LocusEvidence, channel_width};
use crate::read_processing::SangerEvidenceConfig;
use dna_kernel::model::nucleotide::{Nucleotide, is_canonical};
use dna_kernel::read_evidence::{CallEvidence, Placement};

/// The primary peak height of a call with a canonical primary base.
fn primary_height(call: &crate::model::basecalls::BaseCall) -> Option<i32> {
    let channel = Nucleotide::from_char(call.primary)?.channel_index();
    call.peaks.get(channel).map(|peak| peak.height)
}

/// The calls before the read's callable span that its signal vouches for by
/// their base: a canonical primary with no other qualifying channel and a
/// primary peak at least `leading_height_fraction` of the read's median
/// primary peak. A read starts with
/// noise its masked stretch only leaves where a whole window is clean, so
/// strong, unambiguous calls inside that stretch are masked with it; the
/// rule judges each call against the read itself, never a position. Such a
/// call shows its base only, not a run's length or end. A call masked after a
/// repeat run is the read's dephased tail, not start noise, and stays masked.
pub(super) fn leading_identity(
    calls: &BaseCalls,
    callability: &ReadCallability,
    config: &SangerEvidenceConfig,
) -> Vec<usize> {
    let mut heights = calls
        .calls
        .iter()
        .filter_map(primary_height)
        .collect::<Vec<_>>();
    heights.sort_unstable();
    let Some(&median) = heights.get(heights.len() / 2) else {
        return Vec::new();
    };
    let start = callability.callable_start_0based.min(calls.calls.len());
    (0..start)
        .filter(|&index| {
            let call = &calls.calls[index];
            // A masked call, never one where the read loses phase after a
            // repeat run: that is the dephased tail, not start noise.
            placement(callability, index)
                .is_some_and(|placement| placement.reason() != Some(POST_HOMOPOLYMER))
                && call.qualifying_channels.len() == 1
                && primary_height(call).is_some_and(|height| {
                    f64::from(height) >= config.leading_height_fraction * f64::from(median)
                })
        })
        .collect()
}

/// The upper median of `values`: the middle value, the upper of the two
/// middle ones for an even count; `None` for none.
fn median(mut values: Vec<f64>) -> Option<f64> {
    values.sort_by(f64::total_cmp);
    values.get(values.len() / 2).copied()
}

/// The medians a call's own signal is judged against: primary peak, width,
/// and spacing.
#[derive(Debug, Clone, Copy)]
struct Norms {
    height: f64,
    width: f64,
    spacing: f64,
}

impl Norms {
    /// The norms, when each median exists.
    fn new(height: Option<f64>, width: Option<f64>, spacing: Option<f64>) -> Option<Self> {
        Some(Self {
            height: height?,
            width: width?,
            spacing: spacing?,
        })
    }

    /// The medians of the primary peaks and widths of `indices` (in read
    /// order), and of the spacings between those of them that are adjacent.
    fn of(signal: Signal<'_>, indices: &[usize]) -> Option<Self> {
        Self::new(
            median(
                indices
                    .iter()
                    .filter_map(|&index| signal.calls.calls.get(index).and_then(primary_height))
                    .map(f64::from)
                    .collect(),
            ),
            median(
                indices
                    .iter()
                    .filter_map(|&index| signal.width(index))
                    .collect(),
            ),
            median(
                indices
                    .windows(2)
                    .filter(|pair| pair[1] == pair[0] + 1)
                    .filter_map(|pair| signal.gap(pair[0], pair[1]))
                    .collect(),
            ),
        )
    }
}

/// A read's calls and their locus evidence, as the edge and start rules
/// judge a call's own signal.
#[derive(Clone, Copy)]
struct Signal<'a> {
    calls: &'a BaseCalls,
    loci: &'a [LocusEvidence],
}

impl Signal<'_> {
    /// The width of the call's primary peak, when positive.
    fn width(self, index: usize) -> Option<f64> {
        let call = self.calls.calls.get(index)?;
        let channel = Nucleotide::from_char(call.primary)?.channel_index();
        channel_width(self.loci, index, channel).map(|width| width as f64)
    }

    /// Whether the call has a canonical primary and no other qualifying
    /// channel.
    fn single(self, index: usize) -> bool {
        self.calls
            .calls
            .get(index)
            .is_some_and(|call| is_canonical(call.primary) && call.qualifying_channels.len() == 1)
    }

    /// The distance between two calls' loci, in samples.
    fn gap(self, left: usize, right: usize) -> Option<f64> {
        Some(
            self.loci.get(right)?.locus_position_0based as f64
                - self.loci.get(left)?.locus_position_0based as f64,
        )
    }

    /// Whether the call lies `edge_spacing_minimum` to
    /// `edge_spacing_maximum` times `spacing` from `other`; true when either
    /// has no locus.
    fn spaced(
        self,
        index: usize,
        other: usize,
        spacing: f64,
        config: &SangerEvidenceConfig,
    ) -> bool {
        self.gap(index.min(other), index.max(other))
            .is_none_or(|gap| {
                gap >= config.edge_spacing_minimum * spacing
                    && gap <= config.edge_spacing_maximum * spacing
            })
    }

    /// Whether the call's own signal shows its base alone against `norms`:
    /// single-channel, a primary peak of at least `height_fraction` of the
    /// norm, a second-channel share below `second_share`, a primary peak no
    /// wider than `edge_width_factor` times the norm, and spaced from each
    /// neighbour.
    fn unambiguous(
        self,
        evidence: &[CallEvidence],
        index: usize,
        norms: Norms,
        (height_fraction, second_share): (f64, f64),
        config: &SangerEvidenceConfig,
    ) -> bool {
        self.single(index)
            && self
                .calls
                .calls
                .get(index)
                .and_then(primary_height)
                .is_some_and(|primary| f64::from(primary) >= height_fraction * norms.height)
            && evidence[index]
                .second
                .is_none_or(|second| second.share < second_share)
            && self
                .width(index)
                .is_some_and(|width| width <= config.edge_width_factor * norms.width)
            && (index == 0 || self.spaced(index, index - 1, norms.spacing, config))
            && (index + 1 >= evidence.len() || self.spaced(index, index + 1, norms.spacing, config))
    }
}

/// The calls at a signal edge that their own signal vouches for by their
/// base (SRS-CALL-017): scored, with a canonical primary and no other
/// qualifying channel, a second-channel share below `edge_second_share`, a
/// primary peak no wider than `edge_width_factor` times the median width,
/// and a spacing to each neighbour within `edge_spacing_minimum` to
/// `edge_spacing_maximum` times the median. Only calls within
/// `edge_identity_margin` of a stop of the read's in-phase signal qualify,
/// with every call between them and that stop and every call within that
/// margin on their inner side scored:
///
/// - at either informative end, against the read's medians, with a peak of
///   at least `edge_height_fraction` of its median, and not reading the base
///   of its neighbour towards that end (that neighbour's signal one position
///   off);
/// - before a run of one base, at least two calls, that closes a stop and
///   that the call precedes, against the medians of the last
///   `end_norm_calls` scored single-channel calls before that stop, with a
///   peak of at least `end_height_fraction` of that median: a read's signal
///   decays over its last calls, and its peaks vary with their dye context.
///   The stops are the informative end and every call masked
///   `post_homopolymer` right after such a run inside the interval: every
///   molecule shows the same sequence before a run the read loses phase
///   after. The call right before the run (its bound) must be vouched for
///   itself before any call further in qualifies: an ambiguous bound leaves
///   the read's phase at the run in doubt. Judged from the stop inwards, a
///   call reading the base of its
///   neighbour towards the stop qualifies only inside a run that a vouched
///   call of another base closes before the stop, behind a vouched
///   neighbour: a run that the stop or an unvouched call closes may hold
///   one call too many.
///
/// At the informative end no call within the margin beyond it may be masked
/// where the read loses phase (`dephased_signal`, `post_homopolymer`), unless
/// the interval ends with a run that starts after the call and the first
/// call beyond it is masked `post_homopolymer`. The rule judges each call
/// against its own read only, never a position.
pub(super) fn edge_identity(
    calls: &BaseCalls,
    loci: &[LocusEvidence],
    evidence: &[CallEvidence],
    informative: &std::ops::Range<usize>,
    config: &SangerEvidenceConfig,
) -> Vec<usize> {
    let margin = config.edge_identity_margin;
    let count = evidence.len();
    let signal = Signal { calls, loci };
    let single = |index: usize| signal.single(index);
    let read_norms = Norms::new(
        median(
            calls
                .calls
                .iter()
                .filter_map(primary_height)
                .map(f64::from)
                .collect(),
        ),
        median(
            informative
                .clone()
                .filter(|&index| single(index))
                .filter_map(|index| signal.width(index))
                .collect(),
        ),
        median(
            loci.windows(2)
                .map(|pair| {
                    pair[1].locus_position_0based as f64 - pair[0].locus_position_0based as f64
                })
                .collect(),
        ),
    );
    // The norms of the signal before a stop: its last scored single-channel
    // calls.
    let stop_norms = |stop: usize| {
        let mut last = (informative.start..stop)
            .rev()
            .filter(|&index| evidence[index].placement == Placement::Scored && single(index))
            .take(config.end_norm_calls)
            .collect::<Vec<_>>();
        last.reverse();
        Norms::of(signal, &last)
    };
    let scored = |range: std::ops::Range<usize>| {
        range
            .filter(|index| *index < count)
            .all(|index| evidence[index].placement == Placement::Scored)
    };
    let phase_loss = |range: std::ops::Range<usize>| {
        range.filter(|index| *index < count).any(|index| {
            matches!(
                evidence[index].placement.reason(),
                Some(reason) if reason == DEPHASED_SIGNAL || reason == POST_HOMOPOLYMER
            )
        })
    };
    // The first call of the run of one base that ends right before `stop`.
    let run_start = |stop: usize| {
        let base = stop
            .checked_sub(1)
            .and_then(|last| calls.calls.get(last))
            .map(|call| call.primary);
        let mut first = stop;
        while first > informative.start
            && base.is_some_and(|base| {
                calls
                    .calls
                    .get(first - 1)
                    .is_some_and(|call| call.primary == base)
            })
        {
            first -= 1;
        }
        first
    };
    // The read loses phase right after a repeat run that the call precedes:
    // a run of one base, at least two calls, ending at `stop`, starting
    // after the call, with the first call at `stop` masked where the read
    // loses phase after a repeat run. Every molecule shows the same
    // sequence before the run, so the call's signal is in phase.
    let before_phase_loss = |stop: usize, index: usize| {
        let first = run_start(stop);
        evidence
            .get(stop)
            .is_some_and(|beyond| beyond.placement.reason() == Some(POST_HOMOPOLYMER))
            && index < first
            && stop - first >= 2
    };
    let unambiguous = |index: usize, norms: Norms, height_fraction: f64| {
        signal.unambiguous(
            evidence,
            index,
            norms,
            (height_fraction, config.edge_second_share),
            config,
        )
    };
    // A call reading the base of its neighbour towards a stop may be that
    // neighbour's signal one position off.
    let same_base = |index: usize, outer: Option<usize>| {
        outer
            .and_then(|outer| calls.calls.get(outer))
            .zip(calls.calls.get(index))
            .is_some_and(|(outer, call)| outer.primary == call.primary)
    };
    let edge_call = |index: usize| {
        evidence[index].placement == Placement::Scored && evidence[index].identity.edge
    };
    // At the informative start, against the read's norms.
    let start_side = |index: usize| {
        index >= informative.start + margin
            || (scored(informative.start..index)
                && scored(index + 1..index + 1 + margin)
                && !phase_loss(informative.start.saturating_sub(margin)..informative.start)
                && !same_base(index, index.checked_sub(1))
                && read_norms
                    .is_some_and(|norms| unambiguous(index, norms, config.edge_height_fraction)))
    };
    let trailing_from = informative
        .start
        .max(informative.end.saturating_sub(margin));
    let mut accepted = (informative.start..trailing_from)
        .filter(|&index| {
            index < informative.start + margin && edge_call(index) && start_side(index)
        })
        .collect::<Vec<_>>();
    // The stops: the informative end, and every post-repeat phase loss
    // inside the interval right after a run of two calls or more.
    let mut stops = vec![informative.end];
    stops.extend(informative.clone().filter(|&index| {
        index > informative.start
            && evidence[index].placement.reason() == Some(POST_HOMOPOLYMER)
            && evidence[index - 1].placement == Placement::Scored
            && index - run_start(index) >= 2
    }));
    let mut vouched = evidence
        .iter()
        .map(|call| call.placement == Placement::Scored && !call.identity.edge)
        .collect::<Vec<_>>();
    // A call reading its outer neighbour's base, inside a run that a call of
    // another base closes before the stop: both that call and the outer
    // neighbour are vouched for, so the peaks read the run's calls. A run
    // the stop or an unvouched call closes may hold one call too many.
    let bounded_run = |index: usize, stop: usize, vouched: &[bool]| {
        let base = calls.calls.get(index).map(|call| call.primary);
        vouched.get(index + 1).copied().unwrap_or(false)
            && (index + 1..stop)
                .find(|&outer| calls.calls.get(outer).map(|call| call.primary) != base)
                .is_some_and(|outer| vouched[outer])
    };
    for stop in stops {
        let norms = stop_norms(stop);
        let interior = stop < informative.end;
        // The stop closes a run of one base, two calls or more, that the
        // call precedes.
        let first = run_start(stop);
        let before_run = |index: usize| index < first && stop - first >= 2;
        let from = informative.start.max(stop.saturating_sub(margin));
        for index in (from..stop).rev() {
            let placed = edge_call(index)
                && start_side(index)
                && scored(index + 1..stop)
                && scored(index.saturating_sub(margin)..index)
                && if interior {
                    before_phase_loss(stop, index)
                } else {
                    !phase_loss(stop..stop + margin) || before_phase_loss(stop, index)
                };
            // At the informative end against the read's norms, never
            // reading the outer neighbour's base; before a run that closes
            // the stop, behind the run's vouched bound (the call right
            // before it: an ambiguous bound leaves the read's phase at the
            // run in doubt), against the stop's own norms.
            let qualifies = placed
                && ((!interior
                    && !same_base(index, Some(index + 1))
                    && read_norms.is_some_and(|norms| {
                        unambiguous(index, norms, config.edge_height_fraction)
                    }))
                    || (before_run(index)
                        && (index + 1 == first || vouched[first - 1])
                        && (!same_base(index, Some(index + 1))
                            || bounded_run(index, stop, &vouched))
                        && norms.is_some_and(|norms| {
                            unambiguous(index, norms, config.end_height_fraction)
                        })));
            if qualifies {
                vouched[index] = true;
                accepted.push(index);
            }
        }
    }
    accepted
}

/// The calls at a signal edge of the read's start that their own signal
/// vouches for by their base (SRS-CALL-018). The start reaches
/// `start_identity_margin` calls past the callable span's start; its masked calls before the callable span
/// are start noise, so a call may sit behind them, but not behind a masked
/// call of the callable span. A call qualifies under the rule of
/// [`edge_identity`] judged against the start's own norms (the median primary
/// peak, width, and spacing of its first `start_norm_calls` in-phase,
/// single-channel calls; a primary peak of at least `start_height_fraction` of
/// that median, since a read's signal rises over its first calls), with every
/// call within the margin on its inner side scored and not reading the base
/// of the call before it.
pub(super) fn start_identity(
    calls: &BaseCalls,
    loci: &[LocusEvidence],
    evidence: &[CallEvidence],
    informative: &std::ops::Range<usize>,
    callable_start: usize,
    config: &SangerEvidenceConfig,
) -> Vec<usize> {
    let margin = config.start_identity_margin;
    let count = evidence.len();
    let signal = Signal { calls, loci };
    let norm = (callable_start..count)
        .filter(|&index| evidence[index].placement == Placement::Scored && signal.single(index))
        .take(config.start_norm_calls)
        .collect::<Vec<_>>();
    let Some(norms) = Norms::of(signal, &norm) else {
        return Vec::new();
    };
    let end = (callable_start + margin).min(count);
    // A call vouched for its base: scored and off the edge, or accepted by
    // this rule, in read order.
    let mut vouched = evidence
        .iter()
        .map(|call| call.placement == Placement::Scored && !call.identity.edge)
        .collect::<Vec<_>>();
    let mut accepted = Vec::new();
    for index in informative.start..end {
        let call = &evidence[index];
        // A call reading the base of the call before it may be that call's
        // signal one position off, unless the call before is itself vouched
        // for: then the two peaks read two calls of a run.
        let same_base_vouched = index.checked_sub(1).is_none_or(|before| {
            calls
                .calls
                .get(before)
                .zip(calls.calls.get(index))
                .is_none_or(|(before, call)| before.primary != call.primary)
                || vouched[before]
        });
        let qualifies = call.placement == Placement::Scored
            && call.identity.edge
            && (callable_start..index)
                .all(|before| evidence[before].placement == Placement::Scored)
            && (index + 1..(index + 1 + margin).min(count))
                .all(|inner| evidence[inner].placement == Placement::Scored)
            && same_base_vouched
            && signal.unambiguous(
                evidence,
                index,
                norms,
                (config.start_height_fraction, config.start_second_share),
                config,
            );
        if qualifies {
            vouched[index] = true;
            accepted.push(index);
        }
    }
    accepted
}
