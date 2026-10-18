//! The Sanger adapter to the core: one [`ReadEvidence`] per processed read.
//!
//! The core never sees chromatogram types. This adapter turns the Sanger
//! products of one read into the modality-neutral record: the primary calls,
//! each call's basecall-independent locus profile, its placement (scored, or
//! masked with its phase-state reason), the trim interval, and its verdicts:
//! whether it shows its base (peak floor, relative quality, mixed signal),
//! an indel (peak floor, relative quality), and a run's end, which calls lie
//! at an edge of the read's usable signal, the leading calls the signal
//! vouches for by their base alone, the bounded runs before a post-repeat
//! phase loss (ADR-0076), and each call's second base.

use super::SangerEvidenceConfig;
use crate::error::{ProductsError, Result};
use crate::model::basecalls::BaseCalls;
use crate::model::callability::{PhaseState, ReadCallability};
use crate::model::quality::QualityControlResult;
use crate::model::signal::SignalAnalysis;
use dna_kernel::model::nucleotide::Nucleotide;
use dna_kernel::model::nucleotide::is_canonical;
use dna_kernel::read_evidence::{
    Bound, CallEvidence, EvidenceReason, Placement, ReadEvidence, Reasons, SecondBase, Verdict,
};

/// A supporting call's highest A/C/G/T peak is below the configured floor.
pub(crate) const PEAK_BELOW_MINIMUM: EvidenceReason = EvidenceReason::new("peak_below_minimum");
/// A supporting call's relative quality does not strictly exceed the threshold.
pub(crate) const RELATIVE_QUALITY_NOT_ABOVE_THRESHOLD: EvidenceReason =
    EvidenceReason::new("relative_quality_not_above_threshold");
/// An SNV supporting call retains more than one co-localized qualifying channel.
pub(crate) const MIXED_SUPPORTING_DNA: EvidenceReason = EvidenceReason::new("mixed_supporting_dna");
/// A masked call in a segment that starts right after a long repeat run.
pub(crate) const POST_HOMOPOLYMER: EvidenceReason = EvidenceReason::new("post_homopolymer");
/// A masked call in a dephased segment.
pub(crate) const DEPHASED_SIGNAL: EvidenceReason = EvidenceReason::new("dephased_signal");
/// A masked call in a mixed-signal segment.
pub(crate) const MIXED_SIGNAL: EvidenceReason = EvidenceReason::new("mixed_signal");
/// A masked call in a weak-signal segment.
pub(crate) const WEAK_SIGNAL: EvidenceReason = EvidenceReason::new("weak_signal");
/// A masked call in an irregular-spacing segment.
pub(crate) const IRREGULAR_SPACING: EvidenceReason = EvidenceReason::new("irregular_spacing");

/// The Sanger verdict reasons, in reporting order.
const REASON_ORDER: [EvidenceReason; 3] = [
    PEAK_BELOW_MINIMUM,
    RELATIVE_QUALITY_NOT_ABOVE_THRESHOLD,
    MIXED_SUPPORTING_DNA,
];

/// Informative calls between a run's bound and a post-repeat edge for the
/// run to count as bounded before that phase loss (ADR-0076).
const POST_REPEAT_MINIMUM_CALLS: usize = 2;

/// Largest second-channel share of an edge call its own signal vouches for
/// (SRS-CALL-017).
const EDGE_SECOND_SHARE: f64 = 0.10;

/// Largest primary peak width of such a call, as a multiple of the read's
/// median width: a dye blob is broader than the read's peaks.
const EDGE_WIDTH_FACTOR: f64 = 1.5;

/// The spacing to each neighbouring call of such a call, as a band of
/// multiples of the read's median spacing.
const EDGE_SPACING: (f64, f64) = (0.5, 1.5);

/// Builds one read's core evidence from its Sanger products. A call within
/// `edge_margin` calls of the trim interval's ends, of a masked call that
/// does not anchor the alignment, or of the read's physical ends is at an
/// edge of the read's usable signal. The leading calls the signal vouches for
/// by their base (SRS-CALL-014) are then unmasked, trusted for their base
/// only, and freed of the relative-quality reason, and the informative
/// interval starts at the first of them. The edge calls at the interval's
/// ends that their own signal vouches for (SRS-CALL-017) then have their
/// verdict on their base taken off the edge. Last, the bounded runs before a
/// post-repeat phase loss are marked (ADR-0076).
///
/// # Errors
///
/// Returns `ProductsError` when the Sanger products disagree
/// on call count, trim, or call identity, and `EvidenceError` when the evidence
/// is invalid.
pub fn read_evidence(
    calls: &BaseCalls,
    signal: &SignalAnalysis,
    callability: &ReadCallability,
    quality: &QualityControlResult,
    config: &SangerEvidenceConfig,
    edge_margin: usize,
) -> Result<ReadEvidence> {
    let count = quality.per_call.len();
    if callability.mask.len() != count {
        return Err(ProductsError::MaskCountMismatch.into());
    }
    if signal.loci.len() != count {
        return Err(ProductsError::CallCountMismatch {
            loci: signal.loci.len(),
            qualities: count,
        }
        .into());
    }
    let (start, end) = (quality.trim_start_0based, quality.trim_end_0based_exclusive);
    if start > end || end > count {
        return Err(ProductsError::InvalidTrim {
            start,
            end,
            profiles: count,
        }
        .into());
    }
    if quality.retained_sequence.len() != end - start {
        return Err(ProductsError::RetainedLengthMismatch {
            bases: quality.retained_sequence.len(),
            profiles: end - start,
        }
        .into());
    }
    let evidence = (0..count)
        .map(|index| {
            let call = calls
                .calls
                .get(index)
                .ok_or(ProductsError::MissingCall { index })?;
            let score = &quality.per_call[index];
            if call.index_0based != index || score.index_0based != index {
                return Err(ProductsError::CallMismatch { index }.into());
            }
            let highest_peak = call
                .peaks
                .iter()
                .map(|peak| peak.height)
                .max()
                .ok_or(ProductsError::NoChannelPeaks { index })?;
            let mut evidence = CallEvidence {
                second: second_base(call),
                ..CallEvidence::scored(call.primary, signal.loci[index].profile)
            };
            // Peak floor and relative quality withhold the base and an
            // indel; mixed signal withholds the base only.
            if highest_peak < config.minimum_peak_height {
                evidence.identity.withheld.insert(PEAK_BELOW_MINIMUM)?;
                evidence.indel.withheld.insert(PEAK_BELOW_MINIMUM)?;
            }
            if score.relative_quality_score <= config.relative_quality_threshold {
                evidence
                    .identity
                    .withheld
                    .insert(RELATIVE_QUALITY_NOT_ABOVE_THRESHOLD)?;
                evidence
                    .indel
                    .withheld
                    .insert(RELATIVE_QUALITY_NOT_ABOVE_THRESHOLD)?;
            }
            if call.qualifying_channels.len() > 1 {
                evidence.identity.withheld.insert(MIXED_SUPPORTING_DNA)?;
            }
            Ok(match placement(callability, index) {
                Some(placement) => evidence.masked(placement),
                None => evidence,
            })
        })
        .collect::<Result<Vec<_>>>()?;
    let edged = ReadEvidence::new(evidence, start..end, REASON_ORDER.to_vec())?
        .with_signal_edges(edge_margin);
    let leading = leading_identity(calls, callability);
    let mut evidence = edged.calls().to_vec();
    for &index in &leading {
        let call = &mut evidence[index];
        let mut withheld = Reasons::default();
        if call.identity.withheld.contains(PEAK_BELOW_MINIMUM) {
            withheld.insert(PEAK_BELOW_MINIMUM)?;
        }
        *call = CallEvidence {
            placement: Placement::Scored,
            identity: Verdict {
                edge: false,
                withheld,
            },
            indel: Verdict {
                edge: false,
                withheld,
            },
            bound: Bound::Withheld,
            ..*call
        };
    }
    let informative = leading.first().map_or(start, |&first| first.min(start))..end;
    let mut lifted = edge_identity(calls, &signal.loci, &evidence, &informative, edge_margin);
    lifted.extend(start_identity(
        calls,
        &signal.loci,
        &evidence,
        &informative,
        callability.callable_start_0based,
        edge_margin,
    ));
    lifted.sort_unstable();
    lifted.dedup();
    for index in lifted {
        let call = &mut evidence[index];
        let mut withheld = Reasons::default();
        if call.identity.withheld.contains(PEAK_BELOW_MINIMUM) {
            withheld.insert(PEAK_BELOW_MINIMUM)?;
        }
        call.identity = Verdict {
            edge: false,
            withheld,
        };
    }
    let after_repeat = (0..count)
        .map(|index| {
            evidence[index].placement != Placement::Scored
                && callability
                    .segment_at(index)
                    .is_some_and(|segment| segment.after_repeat)
        })
        .collect::<Vec<_>>();
    mark_post_repeat_runs(&mut evidence, &informative, &after_repeat, edge_margin);
    Ok(ReadEvidence::new(
        evidence,
        informative,
        REASON_ORDER.to_vec(),
    )?)
}

/// Marks every run of one base the read bounds right before it loses phase
/// after a repeat run (ADR-0076): a run of at least two canonical calls,
/// scored and counting towards runs, bounded on both sides by calls that
/// show a run's end and read another base, where on each side either more
/// than `margin` informative calls separate the bound from the nearest
/// uninformative call, or at least [`POST_REPEAT_MINIMUM_CALLS`] do and that
/// call is masked after a repeat run (`after_repeat`), while some call of the
/// run or its bounds is at a signal edge or within `margin` of the read's
/// ends. The core trusts such a run's length; it takes the run from the read's
/// own calls, and the alignment-dependent half of the rule stays there.
fn mark_post_repeat_runs(
    calls: &mut [CallEvidence],
    informative: &std::ops::Range<usize>,
    after_repeat: &[bool],
    margin: usize,
) {
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
            || (run >= POST_REPEAT_MINIMUM_CALLS
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

/// The primary peak height of a call with a canonical primary base.
fn primary_height(call: &crate::model::basecalls::BaseCall) -> Option<i32> {
    let channel = Nucleotide::from_char(call.primary)?.channel_index();
    call.peaks.get(channel).map(|peak| peak.height)
}

/// The calls before the read's callable span that its signal vouches for by
/// their base: a canonical primary with no other qualifying channel and a
/// primary peak at least the read's median primary peak. A read starts with
/// noise its masked stretch only leaves where a whole window is clean, so
/// strong, unambiguous calls inside that stretch are masked with it; the
/// rule judges each call against the read itself, never a position. Such a
/// call shows its base only, not a run's length or end.
fn leading_identity(calls: &BaseCalls, callability: &ReadCallability) -> Vec<usize> {
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
            callability.mask.get(index).copied().flatten().is_some()
                && call.qualifying_channels.len() == 1
                && primary_height(call).is_some_and(|height| height >= median)
        })
        .collect()
}

/// The calls at a signal edge of the read's informative ends that their own
/// signal vouches for by their base (SRS-CALL-017): scored, with a canonical
/// primary and no other qualifying channel, a primary peak at least the read's
/// median, a second-channel share below [`EDGE_SECOND_SHARE`], a primary peak
/// no wider than [`EDGE_WIDTH_FACTOR`] times the read's median width, and a
/// spacing to each neighbour within [`EDGE_SPACING`] of the read's median.
/// Only calls within `margin` of an informative end qualify, with every call
/// between them and that end and every call within `margin` on their inner
/// side scored, and no call within `margin` beyond that end masked where the
/// read loses phase (`dephased_signal`, `post_homopolymer`). The rule judges
/// each call against its own read only, never a position.
fn edge_identity(
    calls: &BaseCalls,
    loci: &[crate::model::locus_evidence::LocusEvidence],
    evidence: &[CallEvidence],
    informative: &std::ops::Range<usize>,
    margin: usize,
) -> Vec<usize> {
    let count = evidence.len();
    let median = |mut values: Vec<f64>| {
        values.sort_by(f64::total_cmp);
        values.get(values.len() / 2).copied()
    };
    let width = |index: usize| {
        let call = calls.calls.get(index)?;
        let channel = Nucleotide::from_char(call.primary)?.channel_index();
        let width = *loci.get(index)?.channel_widths.get(channel)?;
        (width > 0).then_some(width as f64)
    };
    let single = |index: usize| {
        calls
            .calls
            .get(index)
            .is_some_and(|call| is_canonical(call.primary) && call.qualifying_channels.len() == 1)
    };
    let (Some(height), Some(median_width), Some(spacing)) = (
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
                .filter_map(width)
                .collect(),
        ),
        median(
            loci.windows(2)
                .map(|pair| {
                    pair[1].locus_position_0based as f64 - pair[0].locus_position_0based as f64
                })
                .collect(),
        ),
    ) else {
        return Vec::new();
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
    let spaced = |index: usize, other: usize| {
        let (Some(left), Some(right)) = (loci.get(index.min(other)), loci.get(index.max(other)))
        else {
            return true;
        };
        let gap = right.locus_position_0based as f64 - left.locus_position_0based as f64;
        gap >= EDGE_SPACING.0 * spacing && gap <= EDGE_SPACING.1 * spacing
    };
    informative
        .clone()
        .filter(|&index| {
            let call = &evidence[index];
            let leading = index < informative.start + margin;
            let trailing = index + margin >= informative.end;
            let ends = (!leading
                || (scored(informative.start..index)
                    && scored(index + 1..index + 1 + margin)
                    && !phase_loss(informative.start.saturating_sub(margin)..informative.start)))
                && (!trailing
                    || (scored(index + 1..informative.end)
                        && scored(index.saturating_sub(margin)..index)
                        && !phase_loss(informative.end..informative.end + margin)));
            // A call reading the base of its neighbour towards the read's
            // end may be that neighbour's signal one position off.
            let outer_same = |outer: Option<usize>| {
                outer
                    .and_then(|outer| calls.calls.get(outer))
                    .zip(calls.calls.get(index))
                    .is_some_and(|(outer, call)| outer.primary == call.primary)
            };
            call.placement == Placement::Scored
                && call.identity.edge
                && (leading || trailing)
                && !(leading && outer_same(index.checked_sub(1)))
                && !(trailing && outer_same(Some(index + 1)))
                && ends
                && single(index)
                && calls
                    .calls
                    .get(index)
                    .and_then(primary_height)
                    .is_some_and(|primary| f64::from(primary) >= height)
                && call
                    .second
                    .is_none_or(|second| second.share < EDGE_SECOND_SHARE)
                && width(index).is_some_and(|width| width <= EDGE_WIDTH_FACTOR * median_width)
                && (index == 0 || spaced(index, index - 1))
                && (index + 1 >= count || spaced(index, index + 1))
        })
        .collect()
}

/// Calls of the read's start whose norms are the start's own: the first
/// in-phase calls of the callable span (SRS-CALL-018).
const START_NORM_CALLS: usize = 30;

/// Smallest primary peak of a start call, as a share of the start's median:
/// a read's signal rises over its first calls.
const START_HEIGHT_SHARE: f64 = 0.5;

/// The calls at a signal edge of the read's start that their own signal
/// vouches for by their base (SRS-CALL-018). The start reaches `margin` calls
/// past the callable span's start; its masked calls before the callable span
/// are start noise, so a call may sit behind them, but not behind a masked
/// call of the callable span. A call qualifies under the rule of
/// [`edge_identity`] judged against the start's own norms (the median primary
/// peak, width, and spacing of its first [`START_NORM_CALLS`] in-phase,
/// single-channel calls), with every call within `margin` on its inner side
/// scored and not reading the base of the call before it.
fn start_identity(
    calls: &BaseCalls,
    loci: &[crate::model::locus_evidence::LocusEvidence],
    evidence: &[CallEvidence],
    informative: &std::ops::Range<usize>,
    callable_start: usize,
    margin: usize,
) -> Vec<usize> {
    let count = evidence.len();
    let median = |mut values: Vec<f64>| {
        values.sort_by(f64::total_cmp);
        values.get(values.len() / 2).copied()
    };
    let single = |index: usize| {
        calls
            .calls
            .get(index)
            .is_some_and(|call| is_canonical(call.primary) && call.qualifying_channels.len() == 1)
    };
    let width = |index: usize| {
        let call = calls.calls.get(index)?;
        let channel = Nucleotide::from_char(call.primary)?.channel_index();
        let width = *loci.get(index)?.channel_widths.get(channel)?;
        (width > 0).then_some(width as f64)
    };
    let norm = (callable_start..count)
        .filter(|&index| evidence[index].placement == Placement::Scored && single(index))
        .take(START_NORM_CALLS)
        .collect::<Vec<_>>();
    let (Some(height), Some(median_width), Some(spacing)) = (
        median(
            norm.iter()
                .filter_map(|&index| calls.calls.get(index).and_then(primary_height))
                .map(f64::from)
                .collect(),
        ),
        median(norm.iter().filter_map(|&index| width(index)).collect()),
        median(
            norm.windows(2)
                .filter(|pair| pair[1] == pair[0] + 1)
                .filter_map(|pair| {
                    Some(
                        loci.get(pair[1])?.locus_position_0based as f64
                            - loci.get(pair[0])?.locus_position_0based as f64,
                    )
                })
                .collect(),
        ),
    ) else {
        return Vec::new();
    };
    let spaced = |index: usize, other: usize| {
        let (Some(left), Some(right)) = (loci.get(index.min(other)), loci.get(index.max(other)))
        else {
            return true;
        };
        let gap = right.locus_position_0based as f64 - left.locus_position_0based as f64;
        gap >= EDGE_SPACING.0 * spacing && gap <= EDGE_SPACING.1 * spacing
    };
    let end = (callable_start + margin).min(count);
    (informative.start..end)
        .filter(|&index| {
            let call = &evidence[index];
            call.placement == Placement::Scored
                && call.identity.edge
                && (callable_start..index)
                    .all(|before| evidence[before].placement == Placement::Scored)
                && (index + 1..(index + 1 + margin).min(count))
                    .all(|inner| evidence[inner].placement == Placement::Scored)
                && index
                    .checked_sub(1)
                    .and_then(|before| calls.calls.get(before))
                    .zip(calls.calls.get(index))
                    .is_none_or(|(before, call)| before.primary != call.primary)
                && single(index)
                && calls
                    .calls
                    .get(index)
                    .and_then(primary_height)
                    .is_some_and(|primary| f64::from(primary) >= START_HEIGHT_SHARE * height)
                && call
                    .second
                    .is_none_or(|second| second.share < EDGE_SECOND_SHARE)
                && width(index).is_some_and(|width| width <= EDGE_WIDTH_FACTOR * median_width)
                && (index == 0 || spaced(index, index - 1))
                && (index + 1 >= count || spaced(index, index + 1))
        })
        .collect()
}

/// The second base of a call's signal: the highest peak of the other three
/// channels, with its share of the primary and that peak together.
fn second_base(call: &crate::model::basecalls::BaseCall) -> Option<SecondBase> {
    let primary = Nucleotide::from_char(call.primary)?.channel_index();
    let first = call.peaks.get(primary)?.height;
    let (channel, second) = call
        .peaks
        .iter()
        .enumerate()
        .filter(|(channel, _)| *channel != primary)
        .map(|(channel, peak)| (channel, peak.height))
        .max_by_key(|(channel, height)| (*height, std::cmp::Reverse(*channel)))?;
    let total = f64::from(first) + f64::from(second);
    (second > 0 && total > 0.0).then(|| SecondBase {
        base: Nucleotide::ALL[channel].as_char(),
        share: f64::from(second) / total,
    })
}

/// The placement of a masked call: `post_homopolymer` when its segment starts
/// right after a repeat run, otherwise its phase state; dephased calls still
/// anchor the alignment, every other masked call is unresolved. `None` for an
/// unmasked call.
fn placement(callability: &ReadCallability, index: usize) -> Option<Placement> {
    let state = callability.mask.get(index).copied().flatten()?;
    if state == PhaseState::InPhase {
        return None;
    }
    let after_repeat = callability
        .segment_at(index)
        .is_some_and(|segment| segment.after_repeat);
    let reason = if after_repeat {
        POST_HOMOPOLYMER
    } else {
        match state {
            PhaseState::Dephased => DEPHASED_SIGNAL,
            PhaseState::Mixed => MIXED_SIGNAL,
            PhaseState::Weak => WEAK_SIGNAL,
            PhaseState::Irregular => IRREGULAR_SPACING,
            PhaseState::InPhase => return None,
        }
    };
    Some(if state == PhaseState::Dephased {
        Placement::Anchoring(reason)
    } else {
        Placement::Unresolved(reason)
    })
}

#[cfg(test)]
mod tests {
    use crate::model::basecalls::{BaseCall, ChannelPeak, PeakSource};
    use crate::model::callability::PhaseSegment;
    use crate::model::locus_evidence::LocusEvidence;
    use crate::model::quality::CallQuality;
    use crate::model::signal::SangerIntegrity;
    use dna_kernel::model::nucleotide::Nucleotide;

    use super::*;

    fn config() -> SangerEvidenceConfig {
        SangerEvidenceConfig {
            minimum_peak_height: 150,
            relative_quality_threshold: 30,
        }
    }

    /// Sanger products for calls with the given highest peaks and relative
    /// scores; every call is in phase and retained.
    fn products(
        peaks: &[i32],
        scores: &[u8],
    ) -> (
        BaseCalls,
        SignalAnalysis,
        ReadCallability,
        QualityControlResult,
    ) {
        let count = peaks.len();
        let calls = peaks
            .iter()
            .enumerate()
            .map(|(index, &height)| BaseCall {
                index_0based: index,
                locus_position_0based: index * 4,
                window_start_0based: index * 4,
                window_end_0based_exclusive: index * 4 + 1,
                peaks: std::array::from_fn(|channel| ChannelPeak {
                    base: Nucleotide::ALL[channel],
                    height,
                    position_0based: index * 4,
                    source: PeakSource::LocalMaximum,
                }),
                primary_peak_evidence: None,
                primary: 'A',
                ambiguity: 'A',
                qualifying_channels: vec![Nucleotide::A],
                vendor_agrees: None,
            })
            .collect::<Vec<_>>();
        let loci = (0..count)
            .map(|index| LocusEvidence {
                call_index_0based: index,
                locus_position_0based: index * 4,
                window_start_0based: index * 4,
                window_end_0based_exclusive: index * 4 + 1,
                context_call_start_0based: index,
                context_call_end_0based_exclusive: index + 1,
                context_sample_start_0based: index * 4,
                context_sample_end_0based_exclusive: index * 4 + 1,
                event_position_0based: index * 4,
                channel_heights: [100, 0, 0, 0],
                channel_baselines: [0.0; 4],
                channel_noise_sigmas: [1.0; 4],
                corrected_amplitudes: [100.0, 0.0, 0.0, 0.0],
                snrs: [100.0, 0.0, 0.0, 0.0],
                profile: None,
                channel_widths: [4, 0, 0, 0],
            })
            .collect();
        let quality = QualityControlResult {
            per_call: scores
                .iter()
                .enumerate()
                .map(|(index, &score)| CallQuality {
                    index_0based: index,
                    penalty: 0,
                    relative_quality_score: score,
                    vendor_quality_applies: false,
                })
                .collect(),
            trim_start_0based: 0,
            trim_end_0based_exclusive: count,
            retained_sequence: "A".repeat(count),
        };
        (
            BaseCalls {
                primary_sequence: "A".repeat(count),
                calls,
            },
            SignalAnalysis {
                integrity: SangerIntegrity {
                    locus_count: count,
                    duplicate_locus_count: 0,
                    vendor_primary_count: None,
                    vendor_quality_count: None,
                    minimum_locus_spacing: None,
                    median_locus_spacing: None,
                    maximum_locus_spacing: None,
                    clipped_channel_samples: 0,
                    maximum_to_median_event_signal_ratio: None,
                },
                loci,
                windows: Vec::new(),
                noisy_regions: Vec::new(),
            },
            ReadCallability::in_phase(count),
            quality,
        )
    }

    #[test]
    fn raises_peak_and_quality_vetoes_at_their_boundaries() -> Result<()> {
        let (calls, signal, callability, quality) = products(&[149, 150, 150], &[31, 30, 31]);
        let evidence = read_evidence(&calls, &signal, &callability, &quality, &config(), 0)?;
        let vetoes = evidence
            .calls()
            .iter()
            .map(|call| {
                (
                    call.identity.withheld.contains(PEAK_BELOW_MINIMUM)
                        && call.indel.withheld.contains(PEAK_BELOW_MINIMUM),
                    call.identity
                        .withheld
                        .contains(RELATIVE_QUALITY_NOT_ABOVE_THRESHOLD)
                        && call
                            .indel
                            .withheld
                            .contains(RELATIVE_QUALITY_NOT_ABOVE_THRESHOLD),
                )
            })
            .collect::<Vec<_>>();
        assert_eq!(vetoes, [(true, false), (false, true), (false, false)]);
        Ok(())
    }

    #[test]
    fn raises_the_mixed_veto_for_substitutions_only() -> Result<()> {
        let (mut calls, signal, callability, quality) = products(&[200], &[31]);
        calls.calls[0].qualifying_channels = vec![Nucleotide::A, Nucleotide::C];
        let evidence = read_evidence(&calls, &signal, &callability, &quality, &config(), 0)?;
        let call = evidence.calls()[0];
        assert!(call.identity.withheld.contains(MIXED_SUPPORTING_DNA));
        assert!(!call.indel.withheld.contains(MIXED_SUPPORTING_DNA));
        Ok(())
    }

    #[test]
    fn maps_phase_segments_to_placements() -> Result<()> {
        let (calls, signal, mut callability, quality) = products(&[200; 25], &[31; 25]);
        let segments = [
            (0, 5, PhaseState::Weak, false),
            (5, 10, PhaseState::InPhase, false),
            (10, 15, PhaseState::Dephased, false),
            (15, 20, PhaseState::Mixed, true),
            (20, 25, PhaseState::Irregular, false),
        ];
        callability.segments = segments
            .iter()
            .map(|&(start, end, state, after_repeat)| PhaseSegment {
                call_start_0based: start,
                call_end_0based_exclusive: end,
                state,
                after_repeat,
                shadow: None,
            })
            .collect();
        for &(start, end, state, _) in &segments {
            callability.mask[start..end].fill((state != PhaseState::InPhase).then_some(state));
        }
        let evidence = read_evidence(&calls, &signal, &callability, &quality, &config(), 0)?;
        let placements = [0, 7, 12, 17, 22].map(|index| {
            (
                evidence.calls()[index].placement,
                evidence.calls()[index].bound,
            )
        });
        assert_eq!(
            placements,
            [
                (Placement::Unresolved(WEAK_SIGNAL), Bound::Withheld),
                (Placement::Scored, Bound::InPhase),
                (Placement::Anchoring(DEPHASED_SIGNAL), Bound::Anchored),
                (Placement::Unresolved(POST_HOMOPOLYMER), Bound::Withheld),
                (Placement::Unresolved(IRREGULAR_SPACING), Bound::Withheld),
            ]
        );
        Ok(())
    }

    #[test]
    fn vouches_for_strong_unambiguous_leading_calls_by_their_base() -> Result<()> {
        // Calls 0-3 lead a mixed stretch; the read's median primary peak is 200.
        let (mut calls, signal, mut callability, mut quality) = products(
            &[300, 100, 300, 300, 200, 200, 200, 200],
            &[25, 25, 25, 25, 40, 40, 40, 40],
        );
        calls.calls[2].qualifying_channels = vec![Nucleotide::A, Nucleotide::C];
        callability.segments = vec![
            PhaseSegment {
                call_start_0based: 0,
                call_end_0based_exclusive: 4,
                state: PhaseState::Mixed,
                after_repeat: false,
                shadow: None,
            },
            PhaseSegment {
                call_start_0based: 4,
                call_end_0based_exclusive: 8,
                state: PhaseState::InPhase,
                after_repeat: false,
                shadow: None,
            },
        ];
        callability.mask[..4].fill(Some(PhaseState::Mixed));
        callability.callable_start_0based = 4;
        quality.trim_start_0based = 4;
        quality.retained_sequence = "A".repeat(4);
        let evidence = read_evidence(&calls, &signal, &callability, &quality, &config(), 1)?;
        let vouched = evidence
            .calls()
            .iter()
            .map(CallEvidence::base_only)
            .collect::<Vec<_>>();
        assert_eq!(
            vouched,
            [true, false, false, true, false, false, false, false]
        );
        let first = evidence.calls()[0];
        assert!(first.placement == Placement::Scored && !first.identity.edge && !first.indel.edge);
        assert!(
            !first
                .identity
                .withheld
                .contains(RELATIVE_QUALITY_NOT_ABOVE_THRESHOLD)
        );
        assert_ne!(evidence.calls()[1].placement, Placement::Scored);
        assert_eq!(evidence.informative(), 0..8);
        // Every other call keeps the edge it had before: call 4 sits next to
        // the trim start.
        assert!(evidence.calls()[4].identity.edge && !evidence.calls()[5].identity.edge);
        Ok(())
    }

    #[test]
    fn vouches_for_unambiguous_edge_calls_by_their_base_only() -> Result<()> {
        // Ten distinct, single-channel calls of equal height; calls 8-9 lie
        // outside the trim, so with a margin of 2 calls 0-1 and 6-7 sit at an
        // edge of the informative interval 0..8.
        let (mut calls, mut signal, callability, mut quality) = products(&[200; 10], &[40; 10]);
        let bases = ['A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'T', 'A'];
        for (call, base) in calls.calls.iter_mut().zip(bases) {
            let channel = Nucleotide::from_char(base).map_or(0, Nucleotide::channel_index);
            call.primary = base;
            call.ambiguity = base;
            call.qualifying_channels = vec![Nucleotide::ALL[channel]];
            for (index, peak) in call.peaks.iter_mut().enumerate() {
                peak.height = if index == channel { 200 } else { 0 };
            }
        }
        for (locus, base) in signal.loci.iter_mut().zip(bases) {
            let channel = Nucleotide::from_char(base).map_or(0, Nucleotide::channel_index);
            locus.channel_widths =
                std::array::from_fn(|index| if index == channel { 4 } else { 0 });
        }
        // Call 6 is as broad as a dye blob; call 7 reads the base of call 8,
        // its neighbour towards the read's end.
        signal.loci[6].channel_widths = [0, 0, 10, 0];
        quality.trim_end_0based_exclusive = 8;
        quality.retained_sequence = "A".repeat(8);
        let evidence = read_evidence(&calls, &signal, &callability, &quality, &config(), 2)?;
        let lifted = evidence
            .calls()
            .iter()
            .map(|call| !call.identity.edge)
            .collect::<Vec<_>>();
        assert_eq!(
            lifted,
            [
                true, true, true, true, true, true, false, false, false, false
            ]
        );
        let first = evidence.calls()[0];
        // Only the verdict on the base is lifted: the call still counts as it
        // did, and supports no indel at the edge.
        assert!(first.indel.edge && first.bound != Bound::Withheld && !first.base_only());
        Ok(())
    }

    #[test]
    fn vouches_for_unambiguous_start_calls_against_the_start_norms() -> Result<()> {
        // Calls 0-3 are start noise before the callable span, which starts at
        // 4; with a margin of 2 the start reaches call 6. The read's signal
        // rises: the start calls stand at half the later calls' height.
        let bases = [
            'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
        ];
        let heights = [
            200, 200, 200, 200, 200, 200, 200, 400, 400, 400, 400, 400, 400, 400, 400, 400,
        ];
        let (mut calls, mut signal, mut callability, quality) = products(&heights, &[40; 16]);
        for ((call, base), height) in calls.calls.iter_mut().zip(bases).zip(heights) {
            let channel = Nucleotide::from_char(base).map_or(0, Nucleotide::channel_index);
            call.primary = base;
            call.ambiguity = base;
            call.qualifying_channels = vec![Nucleotide::ALL[channel]];
            for (index, peak) in call.peaks.iter_mut().enumerate() {
                peak.height = if index == channel { height } else { 0 };
            }
        }
        for (locus, base) in signal.loci.iter_mut().zip(bases) {
            let channel = Nucleotide::from_char(base).map_or(0, Nucleotide::channel_index);
            locus.channel_widths =
                std::array::from_fn(|index| if index == channel { 4 } else { 0 });
        }
        callability.segments = vec![
            PhaseSegment {
                call_start_0based: 0,
                call_end_0based_exclusive: 4,
                state: PhaseState::Weak,
                after_repeat: false,
                shadow: None,
            },
            PhaseSegment {
                call_start_0based: 4,
                call_end_0based_exclusive: 16,
                state: PhaseState::InPhase,
                after_repeat: false,
                shadow: None,
            },
        ];
        callability.mask[..4].fill(Some(PhaseState::Weak));
        callability.callable_start_0based = 4;
        // Call 5 is broader than the start's peaks.
        signal.loci[5].channel_widths = [0, 10, 0, 0];
        let evidence = read_evidence(&calls, &signal, &callability, &quality, &config(), 2)?;
        let lifted = (4..8)
            .map(|index| !evidence.calls()[index].identity.edge)
            .collect::<Vec<_>>();
        // Call 4 sits behind start noise only and reads at half the start's
        // median; call 5 is too broad; calls 6 and 7 lie beyond the edge
        // margin of the noise and were never at an edge.
        assert_eq!(lifted, [true, false, true, true]);
        assert!(evidence.calls()[4].indel.edge);
        Ok(())
    }

    #[test]
    fn gives_each_call_its_second_base() -> Result<()> {
        let (mut calls, signal, callability, quality) = products(&[200], &[31]);
        calls.calls[0].peaks[2].height = 50;
        calls.calls[0].peaks[1].height = 10;
        calls.calls[0].peaks[3].height = 0;
        let evidence = read_evidence(&calls, &signal, &callability, &quality, &config(), 0)?;
        assert_eq!(
            evidence.calls()[0].second,
            Some(SecondBase {
                base: 'G',
                share: 0.2
            })
        );
        Ok(())
    }

    /// Calls of `ACGTACGTAC G CCCCC T AG`, then, unless `tail` is false, six
    /// masked calls after a repeat run (`after_repeat`) or not, and clean
    /// calls; `close` masks the two calls right after the run's bound too.
    fn post_repeat_marks(tail: bool, after_repeat: bool, close: bool) -> Vec<(bool, bool)> {
        let sequence = if tail {
            "ACGTACGTACGCCCCCTAGAAAAAACGTACGTACGTACG"
        } else {
            "ACGTACGTACGCCCCCTAG"
        };
        let masked = |index: usize| {
            tail && ((19..26).contains(&index) || (close && (17..19).contains(&index)))
        };
        let calls = sequence
            .chars()
            .enumerate()
            .map(|(index, base)| {
                let call = CallEvidence::scored(base, None);
                if masked(index) {
                    call.masked(Placement::Unresolved(POST_HOMOPOLYMER))
                } else {
                    call
                }
            })
            .collect::<Vec<_>>();
        let count = calls.len();
        let Ok(evidence) = ReadEvidence::new(calls, 0..count, REASON_ORDER.to_vec()) else {
            return Vec::new();
        };
        let mut calls = evidence.with_signal_edges(8).calls().to_vec();
        let flags = (0..count)
            .map(|index| masked(index) && after_repeat)
            .collect::<Vec<_>>();
        mark_post_repeat_runs(&mut calls, &(0..count), &flags, 8);
        calls
            .iter()
            .map(|call| (call.post_repeat.call, call.post_repeat.bound))
            .collect()
    }

    #[test]
    fn marks_a_bounded_run_before_a_post_repeat_phase_loss() {
        let marks = post_repeat_marks(true, true, false);
        assert!((11..=15).all(|index| marks[index] == (true, false)));
        assert_eq!((marks[10], marks[16]), ((false, true), (false, true)));
        // Too few informative calls between the bound and the edge.
        assert!(!post_repeat_marks(true, true, true)[13].0);
    }

    #[test]
    fn other_edges_mark_no_run() {
        // A masked segment not after a repeat run, or the read's own end.
        assert!(!post_repeat_marks(true, false, false)[13].0);
        assert!(!post_repeat_marks(false, false, false)[13].0);
    }

    /// Pins the published labels of every Sanger reason.
    #[test]
    fn publishes_the_sanger_reason_labels() {
        let labels = [
            PEAK_BELOW_MINIMUM,
            RELATIVE_QUALITY_NOT_ABOVE_THRESHOLD,
            MIXED_SUPPORTING_DNA,
            POST_HOMOPOLYMER,
            DEPHASED_SIGNAL,
            MIXED_SIGNAL,
            WEAK_SIGNAL,
            IRREGULAR_SPACING,
        ]
        .map(EvidenceReason::label);
        assert_eq!(
            labels,
            [
                "peak_below_minimum",
                "relative_quality_not_above_threshold",
                "mixed_supporting_dna",
                "post_homopolymer",
                "dephased_signal",
                "mixed_signal",
                "weak_signal",
                "irregular_spacing",
            ]
        );
    }
}
