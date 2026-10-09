//! The Sanger adapter to the core: one [`ReadEvidence`] per processed read.
//!
//! The core never sees chromatogram types. This adapter turns the Sanger
//! products of one read into the modality-neutral record: the primary calls,
//! each call's basecall-independent locus profile, the callability mask with
//! its phase-state reason, the trim interval, and the support vetoes the
//! Sanger evidence raises against each call (peak floor, relative quality, and
//! mixed signal for substitutions).

use crate::config::SangerEvidenceConfig;
use crate::error::{AlignmentError, Result, VariantError};
use crate::model::basecalls::BaseCalls;
use crate::model::callability::{PhaseState, ReadCallability};
use crate::model::quality::QualityControlResult;
use crate::model::signal::SignalAnalysis;
use crate::read_evidence::{
    CallEvidence, CallMask, EvidenceReason, MaskedAlignment, ReadEvidence, SupportVeto, VetoScope,
    VetoSet,
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

/// The Sanger support-veto vocabulary, in reporting order.
const SUPPORT_VETOES: [SupportVeto; 3] = [
    SupportVeto {
        reason: PEAK_BELOW_MINIMUM,
        scope: VetoScope::SubstitutionsAndInsertions,
    },
    SupportVeto {
        reason: RELATIVE_QUALITY_NOT_ABOVE_THRESHOLD,
        scope: VetoScope::SubstitutionsAndInsertions,
    },
    SupportVeto {
        reason: MIXED_SUPPORTING_DNA,
        scope: VetoScope::Substitutions,
    },
];
const PEAK: usize = 0;
const QUALITY: usize = 1;
const MIXED: usize = 2;

/// Builds one read's core evidence from its Sanger products.
pub(crate) fn read_evidence(
    calls: &BaseCalls,
    signal: &SignalAnalysis,
    callability: &ReadCallability,
    quality: &QualityControlResult,
    config: &SangerEvidenceConfig,
) -> Result<ReadEvidence> {
    let count = quality.per_call.len();
    if callability.mask.len() != count {
        return Err(AlignmentError::Inconsistent(
            "callability mask and quality evidence disagree in call count",
        )
        .into());
    }
    if signal.loci.len() != count {
        return Err(AlignmentError::CallCountMismatch {
            loci: signal.loci.len(),
            qualities: count,
        }
        .into());
    }
    let (start, end) = (quality.trim_start_0based, quality.trim_end_0based_exclusive);
    if start > end || end > count {
        return Err(AlignmentError::InvalidTrim {
            start,
            end,
            profiles: count,
        }
        .into());
    }
    if quality.retained_sequence.len() != end - start {
        return Err(AlignmentError::RetainedLengthMismatch {
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
                .ok_or(VariantError::MissingCall { index })?;
            let score = &quality.per_call[index];
            if call.index_0based != index || score.index_0based != index {
                return Err(VariantError::CallMismatch { index }.into());
            }
            let highest_peak = call
                .peaks
                .iter()
                .map(|peak| peak.height)
                .max()
                .ok_or(VariantError::NoChannelPeaks { index })?;
            let mut vetoes = VetoSet::default();
            if highest_peak < config.minimum_peak_height {
                vetoes.insert(PEAK);
            }
            if score.relative_quality_score <= config.relative_quality_threshold {
                vetoes.insert(QUALITY);
            }
            if call.qualifying_channels.len() > 1 {
                vetoes.insert(MIXED);
            }
            Ok(CallEvidence {
                base: call.primary,
                profile: signal.loci[index].profile,
                mask: call_mask(callability, index),
                vetoes,
            })
        })
        .collect::<Result<Vec<_>>>()?;
    ReadEvidence::new(evidence, start..end, SUPPORT_VETOES.to_vec())
}

/// The mask of a call: `post_homopolymer` when its segment starts right after
/// a repeat run, otherwise its phase state; dephased calls still anchor the
/// alignment, every other masked call is unresolved.
fn call_mask(callability: &ReadCallability, index: usize) -> Option<CallMask> {
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
    let alignment = if state == PhaseState::Dephased {
        MaskedAlignment::Anchoring
    } else {
        MaskedAlignment::Unresolved
    };
    Some(CallMask { alignment, reason })
}

#[cfg(test)]
mod tests {
    use crate::model::basecalls::{BaseCall, ChannelPeak, PeakSource};
    use crate::model::callability::PhaseSegment;
    use crate::model::locus_evidence::LocusEvidence;
    use crate::model::nucleotide::Nucleotide;
    use crate::model::quality::CallQuality;
    use crate::model::signal::SangerIntegrity;
    use crate::model::variant::VariantExclusionReason;

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
        let evidence = read_evidence(&calls, &signal, &callability, &quality, &config())?;
        let vetoes = evidence
            .calls()
            .iter()
            .map(|call| (call.vetoes.contains(PEAK), call.vetoes.contains(QUALITY)))
            .collect::<Vec<_>>();
        assert_eq!(vetoes, [(true, false), (false, true), (false, false)]);
        Ok(())
    }

    #[test]
    fn raises_the_mixed_veto_for_substitutions_only() -> Result<()> {
        let (mut calls, signal, callability, quality) = products(&[200], &[31]);
        calls.calls[0].qualifying_channels = vec![Nucleotide::A, Nucleotide::C];
        let evidence = read_evidence(&calls, &signal, &callability, &quality, &config())?;
        assert!(evidence.calls()[0].vetoes.contains(MIXED));
        assert_eq!(
            evidence.support_vetoes()[MIXED].scope,
            VetoScope::Substitutions
        );
        Ok(())
    }

    #[test]
    fn maps_phase_segments_to_masks() -> Result<()> {
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
        let evidence = read_evidence(&calls, &signal, &callability, &quality, &config())?;
        let masks = [0, 7, 12, 17, 22].map(|index| evidence.calls()[index].mask);
        let expected = [
            Some((MaskedAlignment::Unresolved, WEAK_SIGNAL)),
            None,
            Some((MaskedAlignment::Anchoring, DEPHASED_SIGNAL)),
            Some((MaskedAlignment::Unresolved, POST_HOMOPOLYMER)),
            Some((MaskedAlignment::Unresolved, IRREGULAR_SPACING)),
        ];
        assert_eq!(
            masks.map(|mask| mask.map(|mask| (mask.alignment, mask.reason))),
            expected
        );
        Ok(())
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
        .map(|reason| VariantExclusionReason::Evidence(reason).label());
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
