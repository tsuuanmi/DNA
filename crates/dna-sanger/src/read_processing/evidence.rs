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
use crate::model::second_cause::{SecondCause, SecondPeak};
use crate::model::signal::SignalAnalysis;
use crate::second_cause::second_base;
use dna_kernel::model::nucleotide::Nucleotide;
use dna_kernel::model::nucleotide::is_canonical;
use dna_kernel::read_evidence::{
    Bound, CallEvidence, EvidenceReason, PhaseMixture, Placement, ReadEvidence, Reasons,
    SecondBase, Verdict,
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

/// Every reason label the Sanger modality raises: the verdict reasons, then
/// the mask labels.
const REASONS: [EvidenceReason; 8] = [
    PEAK_BELOW_MINIMUM,
    RELATIVE_QUALITY_NOT_ABOVE_THRESHOLD,
    MIXED_SUPPORTING_DNA,
    POST_HOMOPOLYMER,
    DEPHASED_SIGNAL,
    MIXED_SIGNAL,
    WEAK_SIGNAL,
    IRREGULAR_SPACING,
];

/// The Sanger reason a published `label` names, or `None` when the modality
/// raises no such reason: a recorded read evidence document is read back
/// onto the modality's own labels (PROP-0011).
#[must_use]
pub fn reason(label: &str) -> Option<EvidenceReason> {
    REASONS.into_iter().find(|reason| reason.label() == label)
}

/// Builds one read's core evidence from its Sanger products. A call within
/// `signal_edge_margin` calls of the trim interval's ends, of a masked call that
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
    second: &[Option<SecondPeak>],
    mixtures: &[PhaseMixture],
    config: &SangerEvidenceConfig,
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
            // The call's second peak and its cause (SRS-CALL-021): a cause
            // other than a neighbour's tail or a dye blob leaves it
            // unexplained.
            let peak = second.get(index).copied().flatten();
            let stands_out = peak.is_some_and(|peak| peak.stands_out);
            let artefact = peak.is_some_and(|peak| {
                matches!(
                    peak.cause,
                    SecondCause::NeighbourTail | SecondCause::DyeBlob
                )
            });
            let mut evidence = CallEvidence {
                second: second_base(call).map(|second| SecondBase {
                    stands_out,
                    ..second
                }),
                second_unexplained: peak.is_some() && !artefact,
                ..CallEvidence::scored(call.primary, signal.loci[index].profile)
            };
            // Peak floor and relative quality withhold the base and an
            // indel; mixed signal withholds the base only. A call below the
            // peak floor whose relative quality reaches the rescue level still
            // shows its base: a low absolute peak (a weak read end, a weak
            // dye context) is not a doubtful call when it stands clear of its
            // neighbourhood and of a second base: a low peak with a strong
            // second base is the split signal of two templates.
            if highest_peak < config.minimum_peak_height
                && (score.relative_quality_score < config.peak_rescue_quality
                    || evidence
                        .second
                        .is_some_and(|second| second.share >= config.peak_rescue_second_share))
            {
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
                // A mixed call withheld for nothing else, in phase between
                // in-phase single-channel calls, whose second peak is a
                // neighbour's tail or a dye blob, may still confirm the
                // reference (SRS-CALL-022): its base is withheld all the same.
                let single = |at: Option<usize>| {
                    at.is_some_and(|at| {
                        callability.mask.get(at).is_some_and(Option::is_none)
                            && calls
                                .calls
                                .get(at)
                                .is_some_and(|call| call.qualifying_channels.len() == 1)
                    })
                };
                evidence.reference_only = artefact
                    && evidence.identity.withheld.is_empty()
                    && callability.mask.get(index).is_some_and(Option::is_none)
                    && single(index.checked_sub(1))
                    && single(Some(index + 1));
                evidence.identity.withheld.insert(MIXED_SUPPORTING_DNA)?;
            }
            Ok(match placement(callability, index) {
                Some(placement) => evidence.masked(placement),
                None => evidence,
            })
        })
        .collect::<Result<Vec<_>>>()?;
    let edged = ReadEvidence::new(evidence, start..end, REASON_ORDER.to_vec())?
        .with_signal_edges(config.signal_edge_margin);
    let leading = leading_identity(calls, callability, config);
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
    let mut lifted = edge_identity(calls, &signal.loci, &evidence, &informative, config);
    lifted.extend(start_identity(
        calls,
        &signal.loci,
        &evidence,
        &informative,
        callability.callable_start_0based,
        config,
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
    mark_post_repeat_runs(&mut evidence, &informative, &after_repeat, config);
    reference_only_calls(
        &mut evidence,
        calls,
        callability,
        (start, end),
        config.signal_edge_margin,
    );
    Ok(
        ReadEvidence::new(evidence, informative, REASON_ORDER.to_vec())?
            .with_mixtures(mixtures.to_vec())?,
    )
}

/// Marks the read's other reference-only calls (SRS-CALL-022, ADR-0080), once
/// its signal edges are known, anywhere in the read. A scored single-channel
/// call without a second base standing out or unexplained is reference-only
/// when it is withheld only by the peak floor or relative quality, or lies at
/// a signal edge every cause of which is benign: the untrimmed start of the
/// read, or a weak-signal mask clear of every repeat run by more than the
/// margin (beside a repeat run a weak segment is the dephased tail). The
/// untrimmed 3' end is not benign: there the read's signal decays. A mixed
/// reference-only call at a signal edge loses the mark: the edge stands
/// against it.
fn reference_only_calls(
    evidence: &mut [CallEvidence],
    calls: &BaseCalls,
    callability: &ReadCallability,
    (start, end): (usize, usize),
    margin: usize,
) {
    let count = evidence.len();
    let informative = |evidence: &[CallEvidence], at: usize| {
        (start..end).contains(&at) && !matches!(evidence[at].placement, Placement::Unresolved(_))
    };
    // A weak segment is benign only clear of every repeat run: beside one it
    // is the read's dephased tail, not a weak start.
    let clear_of_repeats = |at: usize| {
        callability.segment_at(at).is_some_and(|segment| {
            callability.repeats.iter().all(|repeat| {
                repeat.call_end_0based_exclusive + margin <= segment.call_start_0based
                    || segment.call_end_0based_exclusive + margin <= repeat.call_start_0based
            })
        })
    };
    let benign = |evidence: &[CallEvidence], at: usize| {
        at < start
            || (evidence[at].placement == Placement::Unresolved(WEAK_SIGNAL)
                && clear_of_repeats(at))
    };
    for index in 0..count {
        let call = evidence[index];
        if call.identity.withheld.contains(MIXED_SUPPORTING_DNA) {
            if call.identity.edge {
                evidence[index].reference_only = false;
            }
            continue;
        }
        let quality = call.identity.withheld.iter().all(|reason| {
            reason == PEAK_BELOW_MINIMUM || reason == RELATIVE_QUALITY_NOT_ABOVE_THRESHOLD
        });
        let single = calls
            .calls
            .get(index)
            .is_some_and(|base| base.qualifying_channels.len() == 1);
        if call.placement != Placement::Scored
            || !single
            || !quality
            || call.second_unexplained
            || call.second.is_some_and(|second| second.stands_out)
            || (call.identity.withheld.is_empty() && !call.identity.edge)
        {
            continue;
        }
        if call.identity.edge {
            // The nearest uninformative call within the margin on each side.
            let left = (index.saturating_sub(margin)..index)
                .rev()
                .find(|&at| !informative(evidence, at));
            let right =
                (index + 1..(index + margin + 1).min(count)).find(|&at| !informative(evidence, at));
            let left_ok = match left {
                Some(at) => benign(evidence, at),
                None => index >= margin,
            };
            let right_ok = match right {
                Some(at) => benign(evidence, at),
                None => index + margin < count,
            };
            if !(left_ok && right_ok) {
                continue;
            }
        }
        evidence[index].reference_only = true;
    }
}

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
fn mark_post_repeat_runs(
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
fn leading_identity(
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

/// The calls at a signal edge of the read's informative ends that their own
/// signal vouches for by their base (SRS-CALL-017): scored, with a canonical
/// primary and no other qualifying channel, a primary peak at least the read's
/// median times `edge_height_fraction`, a second-channel share below
/// `edge_second_share`, a primary peak no wider than `edge_width_factor`
/// times the read's median width, and a spacing to each neighbour within
/// `edge_spacing_minimum` to `edge_spacing_maximum` times the read's median.
/// Only calls within `edge_identity_margin` of an informative end qualify,
/// with every call between them and that end and every call within that
/// margin on their inner side scored, and no call within it beyond that end
/// masked where the read loses phase (`dephased_signal`, `post_homopolymer`),
/// unless the informative interval ends with a run of one base, at least two
/// calls, that starts after the call, and the first call beyond that end is
/// masked `post_homopolymer`: the read loses phase after that run, and
/// before it every molecule shows the same sequence. The rule judges each
/// call against its own read only, never a position.
fn edge_identity(
    calls: &BaseCalls,
    loci: &[crate::model::locus_evidence::LocusEvidence],
    evidence: &[CallEvidence],
    informative: &std::ops::Range<usize>,
    config: &SangerEvidenceConfig,
) -> Vec<usize> {
    let margin = config.edge_identity_margin;
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
    // The read loses phase right after a repeat run that the call precedes:
    // the informative interval ends with a run of one base, at least two
    // calls, starting after the call, and the first call beyond its end is
    // masked where the read loses phase after a repeat run. Every molecule
    // shows the same sequence before the run, so the call's signal is in
    // phase.
    let run_start = {
        let end = informative.end;
        let base = end
            .checked_sub(1)
            .and_then(|last| calls.calls.get(last))
            .map(|call| call.primary);
        let mut first = end;
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
    let before_phase_loss = |index: usize| {
        evidence
            .get(informative.end)
            .is_some_and(|beyond| beyond.placement.reason() == Some(POST_HOMOPOLYMER))
            && index < run_start
            && informative.end - run_start >= 2
    };
    let spaced = |index: usize, other: usize| {
        let (Some(left), Some(right)) = (loci.get(index.min(other)), loci.get(index.max(other)))
        else {
            return true;
        };
        let gap = right.locus_position_0based as f64 - left.locus_position_0based as f64;
        gap >= config.edge_spacing_minimum * spacing && gap <= config.edge_spacing_maximum * spacing
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
                        && (!phase_loss(informative.end..informative.end + margin)
                            || before_phase_loss(index))));
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
                    .is_some_and(|primary| {
                        f64::from(primary) >= config.edge_height_fraction * height
                    })
                && call
                    .second
                    .is_none_or(|second| second.share < config.edge_second_share)
                && width(index)
                    .is_some_and(|width| width <= config.edge_width_factor * median_width)
                && (index == 0 || spaced(index, index - 1))
                && (index + 1 >= count || spaced(index, index + 1))
        })
        .collect()
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
fn start_identity(
    calls: &BaseCalls,
    loci: &[crate::model::locus_evidence::LocusEvidence],
    evidence: &[CallEvidence],
    informative: &std::ops::Range<usize>,
    callable_start: usize,
    config: &SangerEvidenceConfig,
) -> Vec<usize> {
    let margin = config.start_identity_margin;
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
        .take(config.start_norm_calls)
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
        gap >= config.edge_spacing_minimum * spacing && gap <= config.edge_spacing_maximum * spacing
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
            && single(index)
            && calls
                .calls
                .get(index)
                .and_then(primary_height)
                .is_some_and(|primary| f64::from(primary) >= config.start_height_fraction * height)
            && call
                .second
                .is_none_or(|second| second.share < config.start_second_share)
            && width(index).is_some_and(|width| width <= config.edge_width_factor * median_width)
            && (index == 0 || spaced(index, index - 1))
            && (index + 1 >= count || spaced(index, index + 1));
        if qualifies {
            vouched[index] = true;
            accepted.push(index);
        }
    }
    accepted
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
    use crate::model::callability::{PhaseSegment, RepeatRun, RepeatUnit};
    use crate::model::locus_evidence::LocusEvidence;
    use crate::model::quality::CallQuality;
    use crate::model::signal::SangerIntegrity;
    use dna_kernel::model::nucleotide::Nucleotide;

    use super::*;

    #[test]
    fn reads_back_every_label_the_modality_raises_and_no_other() {
        for known in REASONS {
            assert_eq!(reason(known.label()), Some(known));
        }
        assert_eq!(reason("read_end"), None);
        assert_eq!(reason("Weak_Signal"), None);
        assert_eq!(reason(""), None);
    }

    /// The evidence settings with every read-edge margin `margin`.
    fn config(margin: usize) -> SangerEvidenceConfig {
        SangerEvidenceConfig {
            minimum_peak_height: 150,
            relative_quality_threshold: 30,
            peak_rescue_quality: 50,
            peak_rescue_second_share: 0.20,
            signal_edge_margin: margin,
            leading_height_fraction: 1.0,
            post_repeat_margin: margin,
            post_repeat_minimum_calls: 2,
            edge_identity_margin: margin,
            edge_height_fraction: 1.0,
            edge_second_share: 0.10,
            edge_width_factor: 1.5,
            edge_spacing_minimum: 0.5,
            edge_spacing_maximum: 1.5,
            start_identity_margin: margin,
            start_norm_calls: 30,
            start_height_fraction: 0.5,
            start_second_share: 0.10,
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

    // SRS-CALL-022
    #[test]
    fn a_quality_withheld_or_benign_edge_call_is_reference_only() -> Result<()> {
        // Twenty single-channel calls; the trim starts at 4.
        let (calls, _, mut callability, _) = products(&[300; 20], &[40; 20]);
        let scored = |base| CallEvidence::scored(base, None);
        let fresh = || vec![scored('A'); 20];
        let mark = |evidence: &mut Vec<CallEvidence>, callability: &ReadCallability| {
            reference_only_calls(evidence, &calls, callability, (4, 20), 3);
        };
        // Withheld by relative quality only: reference-only.
        let mut evidence = fresh();
        evidence[10]
            .identity
            .withheld
            .insert(RELATIVE_QUALITY_NOT_ABOVE_THRESHOLD)?;
        mark(&mut evidence, &callability);
        assert!(evidence[10].reference_only);
        // A second base standing out keeps it off.
        let mut evidence = fresh();
        evidence[10]
            .identity
            .withheld
            .insert(RELATIVE_QUALITY_NOT_ABOVE_THRESHOLD)?;
        evidence[10].second = Some(SecondBase {
            base: 'C',
            share: 0.3,
            stands_out: true,
        });
        mark(&mut evidence, &callability);
        assert!(!evidence[10].reference_only);
        // At the edge the untrimmed start makes: reference-only.
        let mut evidence = fresh();
        evidence[5].identity.edge = true;
        mark(&mut evidence, &callability);
        assert!(evidence[5].reference_only);
        // Untrimmed from its first call, the read's edge at call 1 is its
        // physical end: not reference-only.
        let mut evidence = fresh();
        evidence[1].identity.edge = true;
        reference_only_calls(&mut evidence, &calls, &callability, (0, 20), 3);
        assert!(!evidence[1].reference_only);
        // Beside a weak-signal mask clear of repeat runs: reference-only;
        // beside one next to a repeat run, the dephased tail: not.
        let weak = Placement::Unresolved(WEAK_SIGNAL);
        callability.segments = vec![PhaseSegment {
            call_start_0based: 14,
            call_end_0based_exclusive: 15,
            state: PhaseState::Weak,
            after_repeat: false,
            shadow: None,
        }];
        let mut evidence = fresh();
        evidence[14] = evidence[14].masked(weak);
        evidence[12].identity.edge = true;
        mark(&mut evidence, &callability);
        assert!(evidence[12].reference_only);
        callability.repeats = vec![RepeatRun {
            call_start_0based: 15,
            call_end_0based_exclusive: 19,
            unit: RepeatUnit::Homopolymer(0),
        }];
        let mut evidence = fresh();
        evidence[14] = evidence[14].masked(weak);
        evidence[12].identity.edge = true;
        mark(&mut evidence, &callability);
        assert!(!evidence[12].reference_only);
        Ok(())
    }

    // SRS-CALL-022
    #[test]
    fn an_isolated_tail_or_blob_mixed_call_is_reference_only() -> Result<()> {
        use crate::model::second_cause::{SecondCause, SecondPeak};
        let (mut calls, signal, callability, quality) = products(&[300; 5], &[40; 5]);
        calls.calls[2].qualifying_channels = vec![Nucleotide::A, Nucleotide::C];
        let peak = |cause| {
            let mut second = vec![None; 5];
            second[2] = Some(SecondPeak {
                base: 'C',
                share: 0.3,
                stands_out: cause != SecondCause::MixedStretch,
                cause,
            });
            second
        };
        let evidence = |calls: &BaseCalls, second: &[Option<SecondPeak>]| {
            read_evidence(
                calls,
                &signal,
                &callability,
                &quality,
                second,
                &[],
                &config(0),
            )
        };
        // A neighbour's tail: withheld as mixed all the same, its run role
        // unchanged, and reference-only; its second base stands out and is
        // explained.
        let tail = evidence(&calls, &peak(SecondCause::NeighbourTail))?;
        let call = tail.calls()[2];
        assert!(call.reference_only && !call.second_unexplained);
        assert!(call.identity.withheld.contains(MIXED_SUPPORTING_DNA));
        assert_eq!(call.bound, Bound::InPhase);
        assert!(call.second.is_some_and(|second| second.stands_out));
        // An unexplained second peak: not reference-only, and flagged.
        let template = evidence(&calls, &peak(SecondCause::SecondTemplate))?;
        assert!(!template.calls()[2].reference_only && template.calls()[2].second_unexplained);
        // A mixed stretch: flagged, but its second base does not stand out.
        let stretch = evidence(&calls, &peak(SecondCause::MixedStretch))?;
        let call = stretch.calls()[2];
        assert!(call.second_unexplained && !call.second.is_some_and(|second| second.stands_out));
        // Beside another mixed call it is not isolated.
        calls.calls[3].qualifying_channels = vec![Nucleotide::A, Nucleotide::G];
        let crowded = evidence(&calls, &peak(SecondCause::DyeBlob))?;
        assert!(!crowded.calls()[2].reference_only);
        Ok(())
    }

    #[test]
    fn raises_peak_and_quality_vetoes_at_their_boundaries() -> Result<()> {
        let (calls, signal, callability, quality) = products(&[149, 150, 150], &[31, 30, 31]);
        let evidence = read_evidence(
            &calls,
            &signal,
            &callability,
            &quality,
            &[],
            &[],
            &config(0),
        )?;
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
    fn a_clear_call_below_the_peak_floor_keeps_its_base_and_indel() -> Result<()> {
        // Below the floor at quality 49, 50, and 50 beside a strong second
        // base; the rescue level is 50 and the second-share bound 0.20.
        let (mut calls, signal, callability, quality) = products(&[149, 149, 149], &[49, 50, 50]);
        for call in &mut calls.calls[..2] {
            for peak in &mut call.peaks[1..] {
                peak.height = 10;
            }
        }
        for peak in &mut calls.calls[2].peaks[1..] {
            peak.height = 60;
        }
        let evidence = read_evidence(
            &calls,
            &signal,
            &callability,
            &quality,
            &[],
            &[],
            &config(0),
        )?;
        let withheld = evidence
            .calls()
            .iter()
            .map(|call| {
                (
                    call.identity.withheld.contains(PEAK_BELOW_MINIMUM),
                    call.indel.withheld.contains(PEAK_BELOW_MINIMUM),
                )
            })
            .collect::<Vec<_>>();
        assert_eq!(withheld, [(true, true), (false, false), (true, true)]);
        Ok(())
    }

    #[test]
    fn raises_the_mixed_veto_for_substitutions_only() -> Result<()> {
        let (mut calls, signal, callability, quality) = products(&[200], &[31]);
        calls.calls[0].qualifying_channels = vec![Nucleotide::A, Nucleotide::C];
        let evidence = read_evidence(
            &calls,
            &signal,
            &callability,
            &quality,
            &[],
            &[],
            &config(0),
        )?;
        let call = evidence.calls()[0];
        assert!(call.identity.withheld.contains(MIXED_SUPPORTING_DNA));
        assert!(!call.indel.withheld.contains(MIXED_SUPPORTING_DNA));
        Ok(())
    }

    // SRS-CALL-011
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
        let evidence = read_evidence(
            &calls,
            &signal,
            &callability,
            &quality,
            &[],
            &[],
            &config(0),
        )?;
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

    // SRS-CALL-013, SRS-CALL-014
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
        let evidence = read_evidence(
            &calls,
            &signal,
            &callability,
            &quality,
            &[],
            &[],
            &config(1),
        )?;
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
        // Masked after a repeat run, the stretch is a dephased tail, not start
        // noise: no call of it is vouched for.
        callability.segments[0].after_repeat = true;
        let evidence = read_evidence(
            &calls,
            &signal,
            &callability,
            &quality,
            &[],
            &[],
            &config(1),
        )?;
        assert!(evidence.calls()[..4].iter().all(|call| {
            call.placement == Placement::Unresolved(POST_HOMOPOLYMER) && !call.base_only()
        }));
        Ok(())
    }

    // SRS-CALL-013, SRS-CALL-017
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
        let evidence = read_evidence(
            &calls,
            &signal,
            &callability,
            &quality,
            &[],
            &[],
            &config(2),
        )?;
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

    // SRS-CALL-017
    #[test]
    fn vouches_for_a_call_before_a_run_where_the_read_loses_phase() -> Result<()> {
        // Twelve single-channel calls; the informative interval 0..10 ends with
        // the run of three G at 7-9, and the read loses phase after it (10-11).
        // With a margin of 4, calls 6-9 sit at the trailing edge.
        let (mut calls, mut signal, mut callability, mut quality) = products(&[200; 12], &[40; 12]);
        let bases = ['A', 'C', 'G', 'T', 'A', 'C', 'A', 'G', 'G', 'G', 'T', 'T'];
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
        quality.trim_end_0based_exclusive = 10;
        quality.retained_sequence = "A".repeat(10);
        let lifted = |callability: &ReadCallability| -> Result<Vec<bool>> {
            let evidence =
                read_evidence(&calls, &signal, callability, &quality, &[], &[], &config(4))?;
            Ok(evidence.calls()[6..10]
                .iter()
                .map(|call| !call.identity.edge)
                .collect())
        };
        callability.segments = vec![
            PhaseSegment {
                call_start_0based: 0,
                call_end_0based_exclusive: 10,
                state: PhaseState::InPhase,
                after_repeat: false,
                shadow: None,
            },
            PhaseSegment {
                call_start_0based: 10,
                call_end_0based_exclusive: 12,
                state: PhaseState::Dephased,
                after_repeat: true,
                shadow: None,
            },
        ];
        callability.mask[10..].fill(Some(PhaseState::Dephased));
        // After a repeat run: the call before the run is in phase; the run's
        // own calls stay at the edge.
        assert_eq!(lifted(&callability)?, [true, false, false, false]);
        // A phase loss that is not after a repeat run lifts nothing.
        callability.segments[1].after_repeat = false;
        assert_eq!(lifted(&callability)?, [false, false, false, false]);
        Ok(())
    }

    /// Sanger products of single-channel calls of `bases` at `heights`, all
    /// four channels narrow, whose calls 0-3 are weak start noise before the
    /// callable span, which starts at 4.
    fn start_products(
        bases: [char; 16],
        heights: [i32; 16],
    ) -> (
        BaseCalls,
        SignalAnalysis,
        ReadCallability,
        QualityControlResult,
    ) {
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
        (calls, signal, callability, quality)
    }

    // SRS-CALL-018
    #[test]
    fn vouches_for_unambiguous_start_calls_against_the_start_norms() -> Result<()> {
        // With a margin of 2 the start reaches call 6. The read's signal
        // rises: the start calls stand at half the later calls' height.
        let (calls, mut signal, callability, quality) = start_products(
            [
                'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
            ],
            [
                200, 200, 200, 200, 200, 200, 200, 400, 400, 400, 400, 400, 400, 400, 400, 400,
            ],
        );
        // Call 5 is broader than the start's peaks.
        signal.loci[5].channel_widths = [0, 10, 0, 0];
        let evidence = read_evidence(
            &calls,
            &signal,
            &callability,
            &quality,
            &[],
            &[],
            &config(2),
        )?;
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

    // SRS-CALL-018
    #[test]
    fn a_start_call_reads_its_predecessors_base_only_behind_a_vouched_call() -> Result<()> {
        // Calls 4-6 read `TTT` at the start of the callable span; with a
        // margin of 3 the start reaches call 6.
        let bases = [
            'A', 'C', 'G', 'T', 'T', 'T', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A',
        ];
        // The start noise stands too low for the leading rule to vouch for it.
        let heights = [
            100, 100, 100, 100, 400, 400, 400, 400, 400, 400, 400, 400, 400, 400, 400, 400,
        ];
        let (calls, signal, callability, quality) = start_products(bases, heights);
        let lifted = |evidence: &ReadEvidence| {
            (4..7)
                .map(|index| !evidence.calls()[index].identity.edge)
                .collect::<Vec<_>>()
        };
        // Call 4 follows start noise reading `T` but not vouched for, so it
        // may be that call's signal one position off; calls 5 and 6 each
        // follow a call reading `T` that is.
        let evidence = read_evidence(
            &calls,
            &signal,
            &callability,
            &quality,
            &[],
            &[],
            &config(3),
        )?;
        assert_eq!(lifted(&evidence), [false, false, false]);
        // With call 3 reading `G`, call 4 qualifies, and the run chains on.
        let mut bases = bases;
        bases[3] = 'G';
        let (calls, signal, callability, quality) = start_products(bases, heights);
        let evidence = read_evidence(
            &calls,
            &signal,
            &callability,
            &quality,
            &[],
            &[],
            &config(3),
        )?;
        assert_eq!(lifted(&evidence), [true, true, true]);
        // A broad call 5 breaks the chain: call 6 no longer follows a vouched
        // call of its base.
        let (calls, mut signal, callability, quality) = start_products(bases, heights);
        signal.loci[5].channel_widths = [0, 0, 0, 10];
        let evidence = read_evidence(
            &calls,
            &signal,
            &callability,
            &quality,
            &[],
            &[],
            &config(3),
        )?;
        assert_eq!(lifted(&evidence), [true, false, false]);
        Ok(())
    }

    // SRS-CALL-015
    #[test]
    fn gives_each_call_its_second_base() -> Result<()> {
        let (mut calls, signal, callability, quality) = products(&[200], &[31]);
        calls.calls[0].peaks[2].height = 50;
        calls.calls[0].peaks[1].height = 10;
        calls.calls[0].peaks[3].height = 0;
        let evidence = read_evidence(
            &calls,
            &signal,
            &callability,
            &quality,
            &[],
            &[],
            &config(0),
        )?;
        assert_eq!(
            evidence.calls()[0].second,
            Some(SecondBase {
                base: 'G',
                share: 0.2,
                stands_out: false,
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
        mark_post_repeat_runs(&mut calls, &(0..count), &flags, &config(8));
        calls
            .iter()
            .map(|call| (call.post_repeat.call, call.post_repeat.bound))
            .collect()
    }

    // SRS-CALL-016
    #[test]
    fn marks_a_bounded_run_before_a_post_repeat_phase_loss() {
        let marks = post_repeat_marks(true, true, false);
        assert!((11..=15).all(|index| marks[index] == (true, false)));
        assert_eq!((marks[10], marks[16]), ((false, true), (false, true)));
        // Too few informative calls between the bound and the edge.
        assert!(!post_repeat_marks(true, true, true)[13].0);
    }

    // SRS-CALL-016
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
