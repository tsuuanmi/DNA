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

mod identity;
mod post_repeat;
mod reference_only;
#[cfg(test)]
mod tests;

use super::{ProcessedRead, SangerEvidenceConfig};
use crate::error::{ProductsError, Result};
use crate::model::callability::{PhaseState, ReadCallability};
use crate::model::second_cause::{SecondCause, SecondPeak};
use crate::second_cause::second_base;
use dna_kernel::read_evidence::{
    Bound, CallEvidence, EvidenceReason, Placement, ReadEvidence, Reasons, SecondBase,
    UnreliableSpan, Verdict,
};
use identity::{edge_identity, leading_identity, start_identity};
use post_repeat::mark_post_repeat_runs;
use reference_only::reference_only_calls;

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
/// The read's start, before its signal settles (SRS-CALL-024).
pub(crate) const POOR_START_SIGNAL: EvidenceReason = EvidenceReason::new("poor_start_signal");
/// A call inside an electrophoretic spike (SRS-CALL-027).
pub(crate) const SIGNAL_SPIKE: EvidenceReason = EvidenceReason::new("signal_spike");

/// The Sanger verdict reasons, in reporting order.
const REASON_ORDER: [EvidenceReason; 4] = [
    PEAK_BELOW_MINIMUM,
    RELATIVE_QUALITY_NOT_ABOVE_THRESHOLD,
    MIXED_SUPPORTING_DNA,
    SIGNAL_SPIKE,
];

/// Every reason label the Sanger modality raises: the verdict reasons, then
/// the mask labels.
const REASONS: [EvidenceReason; 10] = [
    PEAK_BELOW_MINIMUM,
    RELATIVE_QUALITY_NOT_ABOVE_THRESHOLD,
    MIXED_SUPPORTING_DNA,
    SIGNAL_SPIKE,
    POST_HOMOPOLYMER,
    DEPHASED_SIGNAL,
    MIXED_SIGNAL,
    WEAK_SIGNAL,
    IRREGULAR_SPACING,
    POOR_START_SIGNAL,
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
pub fn read_evidence(read: &ProcessedRead, config: &SangerEvidenceConfig) -> Result<ReadEvidence> {
    let ProcessedRead {
        calls,
        signal,
        callability,
        quality,
        second,
        mixtures,
        tracts,
        warnings: _,
    } = read;
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
            let stands_out = peak.is_some_and(|peak| peak.second.stands_out);
            let artefact = peak.is_some_and(explained);
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
                // A mixed call withheld for nothing else than the peak floor
                // or its relative quality, in phase between neighbours that
                // show one base, whose second peak is a neighbour's tail or a
                // dye blob, may still confirm the reference (SRS-CALL-022):
                // its base is withheld all the same. A neighbour shows one
                // base when it has one qualifying channel, or when its own
                // second peak is a neighbour's tail or a dye blob: a blob is
                // wider than one call, so its neighbours are mixed too. A
                // neighbour before the trim interval may be masked as weak:
                // the read's untrimmed start.
                let single = |at: Option<usize>| {
                    at.is_some_and(|at| {
                        callability.mask.get(at).is_some_and(|state| {
                            state.is_none() || (at < start && *state == Some(PhaseState::Weak))
                        }) && calls.calls.get(at).is_some_and(|call| {
                            call.qualifying_channels.len() == 1
                                || second.get(at).copied().flatten().is_some_and(explained)
                        })
                    })
                };
                evidence.reference_only = artefact
                    && quality_vetoes_only(&evidence.identity.withheld)
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
        let verdict = vouched(&call.identity)?;
        *call = CallEvidence {
            placement: Placement::Scored,
            identity: verdict,
            indel: verdict,
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
        call.identity = vouched(&call.identity)?;
    }
    // A call inside an electrophoretic spike reads the spike's channel, not
    // the template's base (SRS-CALL-027): it shows no base, places no indel,
    // shows no run's end, and confirms nothing. No earlier rule lifts this.
    for spike in &signal.spikes {
        for call in evidence
            .iter_mut()
            .take(spike.calls.end)
            .skip(spike.calls.start)
        {
            call.identity.withheld.insert(SIGNAL_SPIKE)?;
            call.indel.withheld.insert(SIGNAL_SPIKE)?;
            call.bound = Bound::Withheld;
            call.reference_only = false;
        }
    }
    let after_repeat = (0..count)
        .map(|index| {
            evidence[index].placement != Placement::Scored
                && callability
                    .segment_at(index)
                    .is_some_and(|segment| segment.after_repeat)
        })
        .collect::<Vec<_>>();
    // A junction whose spacing marks a missing or an extra peak
    // (SRS-CALL-026): a read that lost a peak there, or called one twice,
    // shows a deletion or an insertion that is not in the template, so
    // neither call beside it places an indel or shows a run's end.
    for &junction in &callability.gaps {
        for call in evidence.iter_mut().skip(junction).take(2) {
            call.bound = Bound::Withheld;
        }
    }
    mark_post_repeat_runs(&mut evidence, &informative, &after_repeat, config);
    reference_only_calls(
        &mut evidence,
        calls,
        callability,
        (start, end),
        config.signal_edge_margin,
    );
    let poor = poor_start(&evidence, config.start_identity_margin);
    Ok(
        ReadEvidence::new(evidence, informative, REASON_ORDER.to_vec())?
            .with_mixtures(mixtures.clone())?
            .with_unreliable(poor)?
            .with_tracts(tracts.clone())?,
    )
}

/// Whether a second peak's cause is the call's own signal: a neighbour's tail
/// or a dye blob (SRS-CALL-021).
fn explained(peak: SecondPeak) -> bool {
    matches!(peak.cause, SecondCause::NeighborTail | SecondCause::DyeBlob)
}

/// The verdict on a call its own signal vouches for: off the edge, and
/// withheld for nothing but the peak floor, if `verdict` was.
fn vouched(verdict: &Verdict) -> Result<Verdict> {
    let mut withheld = Reasons::default();
    if verdict.withheld.contains(PEAK_BELOW_MINIMUM) {
        withheld.insert(PEAK_BELOW_MINIMUM)?;
    }
    Ok(Verdict {
        edge: false,
        withheld,
    })
}

/// Whether a call is withheld for nothing but the peak floor or its relative
/// quality, if at all.
pub(super) fn quality_vetoes_only(withheld: &Reasons) -> bool {
    withheld.iter().all(|reason| {
        reason == PEAK_BELOW_MINIMUM || reason == RELATIVE_QUALITY_NOT_ABOVE_THRESHOLD
    })
}

/// The read's poor start (SRS-CALL-024): from its first call to the first run
/// of `margin` consecutive sound calls, a sound call being scored, clear of a
/// signal edge for its base, and withheld for nothing. The vouching rules
/// have already judged each call against the read's own heights, widths, and
/// spacing, so the span needs no level of its own. None when the read's first
/// `margin` calls are sound; the whole read when no such run exists.
fn poor_start(evidence: &[CallEvidence], margin: usize) -> Vec<UnreliableSpan> {
    let sound =
        |call: &CallEvidence| call.placement == Placement::Scored && call.identity.trusted();
    let mut run = 0;
    let mut settled = None;
    for (index, call) in evidence.iter().enumerate() {
        run = if sound(call) { run + 1 } else { 0 };
        if run >= margin.max(1) {
            settled = Some(index + 1 - run);
            break;
        }
    }
    let end = settled.unwrap_or(evidence.len());
    if end == 0 {
        return Vec::new();
    }
    vec![UnreliableSpan {
        calls: 0..end,
        reason: POOR_START_SIGNAL,
    }]
}

/// The placement of a masked call: `post_homopolymer` when its segment starts
/// right after a repeat run, otherwise its phase state; dephased calls still
/// anchor the alignment, every other masked call is unresolved. `None` for an
/// unmasked call.
pub(super) fn placement(callability: &ReadCallability, index: usize) -> Option<Placement> {
    let state = callability.mask.get(index).copied().flatten()?;
    let reason = match state {
        PhaseState::InPhase => return None,
        PhaseState::Dephased => DEPHASED_SIGNAL,
        PhaseState::Mixed => MIXED_SIGNAL,
        PhaseState::Weak => WEAK_SIGNAL,
        PhaseState::Irregular => IRREGULAR_SPACING,
    };
    let after_repeat = callability
        .segment_at(index)
        .is_some_and(|segment| segment.after_repeat);
    let reason = if after_repeat {
        POST_HOMOPOLYMER
    } else {
        reason
    };
    Some(if state == PhaseState::Dephased {
        Placement::Anchoring(reason)
    } else {
        Placement::Unresolved(reason)
    })
}
