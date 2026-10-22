//! Tests of the Sanger evidence adapter.

use crate::model::basecalls::{BaseCall, BaseCalls, ChannelPeak, PeakSource};
use crate::model::callability::{PhaseSegment, RepeatRun, RepeatUnit};
use crate::model::locus_evidence::LocusEvidence;
use crate::model::quality::{CallQuality, QualityControlResult};
use crate::model::signal::{SangerIntegrity, SignalAnalysis};
use crate::read_processing::ReadWarnings;
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
/// The evidence of a read with these products and no mixtures or tracts.
fn evidence_of(
    calls: &BaseCalls,
    signal: &SignalAnalysis,
    callability: &ReadCallability,
    quality: &QualityControlResult,
    second: &[Option<SecondPeak>],
    config: &SangerEvidenceConfig,
) -> Result<ReadEvidence> {
    let read = ProcessedRead {
        calls: calls.clone(),
        signal: signal.clone(),
        callability: callability.clone(),
        quality: quality.clone(),
        second: second.to_vec(),
        mixtures: Vec::new(),
        tracts: Vec::new(),
        warnings: ReadWarnings {
            unresolved_primary_calls: 0,
            multi_channel_unresolved_calls: 0,
            vendor_disagreements: 0,
            locus_vendor_length_mismatches: 0,
            clipped_channel_samples: 0,
        },
    };
    read_evidence(&read, config)
}

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
        end_norm_calls: 30,
        end_height_fraction: 0.5,
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
            spikes: Vec::new(),
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
    // A second base standing out that the call's own signal explains
    // (a neighbour's tail or a dye blob) leaves it reference-only; an
    // unexplained one keeps it off.
    for unexplained in [false, true] {
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
        evidence[10].second_unexplained = unexplained;
        mark(&mut evidence, &callability);
        assert_eq!(evidence[10].reference_only, !unexplained);
    }
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
            second: SecondBase {
                base: 'C',
                share: 0.3,
                stands_out: cause != SecondCause::MixedStretch,
            },
            cause,
        });
        second
    };
    let evidence = |calls: &BaseCalls, second: &[Option<SecondPeak>]| {
        evidence_of(calls, &signal, &callability, &quality, second, &config(0))
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
    // A mixed neighbour whose own second peak is a neighbour's tail or a
    // dye blob shows one base all the same: a blob is wider than one call.
    let neighbour = |cause| {
        let mut second = peak(SecondCause::DyeBlob);
        second[3] = Some(SecondPeak {
            second: SecondBase {
                base: 'G',
                share: 0.3,
                stands_out: true,
            },
            cause,
        });
        second
    };
    for cause in [SecondCause::NeighbourTail, SecondCause::DyeBlob] {
        let blob = evidence(&calls, &neighbour(cause))?;
        assert!(blob.calls()[2].reference_only, "{cause:?}");
    }
    for cause in [
        SecondCause::SecondTemplate,
        SecondCause::PhaseShift,
        SecondCause::OffCentre,
    ] {
        let other = evidence(&calls, &neighbour(cause))?;
        assert!(!other.calls()[2].reference_only, "{cause:?}");
    }
    Ok(())
}

/// The causes whose stand-out second peak the call's own signal explains,
/// and those it does not.
const EXPLAINED: [SecondCause; 2] = [SecondCause::NeighbourTail, SecondCause::DyeBlob];
const UNEXPLAINED: [SecondCause; 3] = [
    SecondCause::OffCentre,
    SecondCause::SecondTemplate,
    SecondCause::PhaseShift,
];

/// A stand-out second peak of `cause` at call `index` of `count`.
fn stand_out(count: usize, index: usize, cause: SecondCause) -> Vec<Option<SecondPeak>> {
    let mut second = vec![None; count];
    second[index] = Some(SecondPeak {
        second: SecondBase {
            base: 'C',
            share: 0.3,
            stands_out: true,
        },
        cause,
    });
    second
}

// SRS-CALL-022
#[test]
fn a_quality_withheld_call_with_a_tail_or_blob_second_is_reference_only() -> Result<()> {
    // Call 3 of 7, single-channel, withheld by its relative quality only.
    let (calls, signal, callability, quality) = products(&[300; 7], &[40, 40, 40, 30, 40, 40, 40]);
    let evidence = |cause| {
        evidence_of(
            &calls,
            &signal,
            &callability,
            &quality,
            &stand_out(7, 3, cause),
            &config(0),
        )
    };
    for cause in EXPLAINED {
        let call = evidence(cause)?.calls()[3];
        assert!(call.second.is_some_and(|second| second.stands_out));
        assert!(call.reference_only && !call.second_unexplained, "{cause:?}");
        assert!(
            call.identity
                .withheld
                .contains(RELATIVE_QUALITY_NOT_ABOVE_THRESHOLD)
        );
    }
    // An off-centre, second-template, or phase-shift second never
    // qualifies: it may be a second sequence.
    for cause in UNEXPLAINED {
        let call = evidence(cause)?.calls()[3];
        assert!(!call.reference_only && call.second_unexplained, "{cause:?}");
    }
    Ok(())
}

// SRS-CALL-022
#[test]
fn a_mixed_reference_only_call_keeps_the_mark_only_at_a_benign_edge() -> Result<()> {
    // Twenty calls; call 5 or 1 is a mixed tail call already marked
    // reference-only by its own signal, at a signal edge.
    let (calls, _, mut callability, _) = products(&[300; 20], &[40; 20]);
    let mixed = |at: usize| -> Result<Vec<CallEvidence>> {
        let mut evidence = vec![CallEvidence::scored('A', None); 20];
        evidence[at]
            .identity
            .withheld
            .insert(MIXED_SUPPORTING_DNA)?;
        evidence[at].identity.edge = true;
        evidence[at].reference_only = true;
        Ok(evidence)
    };
    // Beside the untrimmed start: the mark stays.
    let mut evidence = mixed(5)?;
    reference_only_calls(&mut evidence, &calls, &callability, (4, 20), 3);
    assert!(evidence[5].reference_only);
    // At the read's physical end, never benign: the mark goes.
    let mut evidence = mixed(1)?;
    reference_only_calls(&mut evidence, &calls, &callability, (0, 20), 3);
    assert!(!evidence[1].reference_only);
    // Beside the untrimmed 3' end: the mark goes.
    let mut evidence = mixed(14)?;
    reference_only_calls(&mut evidence, &calls, &callability, (0, 16), 3);
    assert!(!evidence[14].reference_only);
    // Beside a weak mask clear of repeat runs the mark stays; beside one
    // next to a repeat run, the dephased tail, it goes.
    let weak = Placement::Unresolved(WEAK_SIGNAL);
    callability.segments = vec![PhaseSegment {
        call_start_0based: 14,
        call_end_0based_exclusive: 15,
        state: PhaseState::Weak,
        after_repeat: false,
        shadow: None,
    }];
    let mut evidence = mixed(12)?;
    evidence[14] = evidence[14].masked(weak);
    reference_only_calls(&mut evidence, &calls, &callability, (0, 20), 3);
    assert!(evidence[12].reference_only);
    callability.repeats = vec![RepeatRun {
        call_start_0based: 15,
        call_end_0based_exclusive: 19,
        unit: RepeatUnit::Homopolymer(0),
    }];
    let mut evidence = mixed(12)?;
    evidence[14] = evidence[14].masked(weak);
    reference_only_calls(&mut evidence, &calls, &callability, (0, 20), 3);
    assert!(!evidence[12].reference_only);
    // An edge never marks a mixed call its own signal did not mark.
    let mut evidence = mixed(5)?;
    evidence[5].reference_only = false;
    reference_only_calls(&mut evidence, &calls, &callability, (4, 20), 3);
    assert!(!evidence[5].reference_only);
    Ok(())
}

// SRS-CALL-022
#[test]
fn a_quality_withheld_mixed_tail_call_is_reference_only() -> Result<()> {
    let (mut calls, signal, mut callability, mut quality) =
        products(&[300, 300, 300, 300, 300], &[40, 40, 30, 40, 40]);
    calls.calls[2].qualifying_channels = vec![Nucleotide::A, Nucleotide::C];
    let evidence = |calls: &BaseCalls,
                    callability: &ReadCallability,
                    quality: &QualityControlResult,
                    cause| {
        evidence_of(
            calls,
            &signal,
            callability,
            quality,
            &stand_out(5, 2, cause),
            &config(0),
        )
    };
    // Withheld by its relative quality as well as mixed: reference-only
    // for a tail or blob, never for an unexplained second.
    for cause in EXPLAINED {
        let call = evidence(&calls, &callability, &quality, cause)?.calls()[2];
        assert!(call.reference_only, "{cause:?}");
        assert!(call.identity.withheld.contains(MIXED_SUPPORTING_DNA));
        assert!(
            call.identity
                .withheld
                .contains(RELATIVE_QUALITY_NOT_ABOVE_THRESHOLD)
        );
    }
    for cause in UNEXPLAINED {
        let call = evidence(&calls, &callability, &quality, cause)?.calls()[2];
        assert!(!call.reference_only && call.second_unexplained, "{cause:?}");
    }
    // Below the peak floor as well: reference-only all the same.
    let (mut low, _, _, clear) = products(&[300, 300, 100, 300, 300], &[40; 5]);
    low.calls[2].qualifying_channels = vec![Nucleotide::A, Nucleotide::C];
    let call = evidence(&low, &callability, &clear, SecondCause::NeighbourTail)?.calls()[2];
    assert!(call.reference_only && call.identity.withheld.contains(PEAK_BELOW_MINIMUM));
    // A neighbour masked as weak before the trim interval, the read's
    // untrimmed start, still lets it confirm the reference.
    callability.mask[1] = Some(PhaseState::Weak);
    quality.trim_start_0based = 2;
    quality.retained_sequence = "A".repeat(3);
    let call = evidence(&calls, &callability, &quality, SecondCause::NeighbourTail)?.calls()[2];
    assert!(call.reference_only);
    // Masked as weak inside the trim interval, or masked for another
    // reason before it, the neighbour does not.
    quality.trim_start_0based = 0;
    quality.retained_sequence = "A".repeat(5);
    let call = evidence(&calls, &callability, &quality, SecondCause::NeighbourTail)?.calls()[2];
    assert!(!call.reference_only);
    callability.mask[1] = Some(PhaseState::Mixed);
    quality.trim_start_0based = 2;
    quality.retained_sequence = "A".repeat(3);
    let call = evidence(&calls, &callability, &quality, SecondCause::NeighbourTail)?.calls()[2];
    assert!(!call.reference_only);
    Ok(())
}

// SRS-CALL-027
/// The calls inside a spike show no base, place no indel, show no run's
/// end, and confirm nothing; the calls beside it keep their verdicts.
#[test]
fn the_calls_inside_a_spike_show_nothing() -> Result<()> {
    let (calls, mut signal, callability, quality) = products(&[300; 8], &[40; 8]);
    signal.spikes = vec![crate::model::signal::SignalSpike {
        samples: 0..1,
        calls: 3..5,
    }];
    let evidence = evidence_of(&calls, &signal, &callability, &quality, &[], &config(0))?;
    for (index, call) in evidence.calls().iter().enumerate() {
        let inside = (3..5).contains(&index);
        assert_eq!(call.identity.withheld.contains(SIGNAL_SPIKE), inside);
        assert_eq!(call.indel.withheld.contains(SIGNAL_SPIKE), inside);
        assert_eq!(call.bound == Bound::Withheld, inside);
        assert!(!call.reference_only);
    }
    Ok(())
}

// SRS-CALL-026
#[test]
fn the_calls_beside_a_gap_show_their_base_only() -> Result<()> {
    let (calls, signal, mut callability, quality) = products(&[300; 8], &[40; 8]);
    callability.gaps = vec![3];
    let evidence = evidence_of(&calls, &signal, &callability, &quality, &[], &config(0))?;
    let base_only = evidence
        .calls()
        .iter()
        .map(CallEvidence::base_only)
        .collect::<Vec<_>>();
    assert_eq!(
        base_only,
        [false, false, false, true, true, false, false, false]
    );
    // Their bases stay theirs.
    assert!(evidence.calls()[3].identity.trusted());
    assert!(evidence.calls()[4].identity.trusted());
    Ok(())
}

#[test]
fn raises_peak_and_quality_vetoes_at_their_boundaries() -> Result<()> {
    let (calls, signal, callability, quality) = products(&[149, 150, 150], &[31, 30, 31]);
    let evidence = evidence_of(&calls, &signal, &callability, &quality, &[], &config(0))?;
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
    let evidence = evidence_of(&calls, &signal, &callability, &quality, &[], &config(0))?;
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
    let evidence = evidence_of(&calls, &signal, &callability, &quality, &[], &config(0))?;
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
    let evidence = evidence_of(&calls, &signal, &callability, &quality, &[], &config(0))?;
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
    let evidence = evidence_of(&calls, &signal, &callability, &quality, &[], &config(1))?;
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
    let evidence = evidence_of(&calls, &signal, &callability, &quality, &[], &config(1))?;
    assert!(evidence.calls()[..4].iter().all(|call| {
        call.placement == Placement::Unresolved(POST_HOMOPOLYMER) && !call.base_only()
    }));
    Ok(())
}

// SRS-CALL-024
#[test]
fn a_poor_start_runs_to_the_first_run_of_sound_calls() {
    let sound = CallEvidence::scored('A', None);
    let mut edge = sound;
    edge.identity.edge = true;
    let masked = sound.masked(Placement::Unresolved(WEAK_SIGNAL));
    // Two masked calls, three sound, one at an edge, then four sound.
    let calls = [
        masked, masked, sound, sound, sound, edge, sound, sound, sound, sound,
    ];
    let poor = poor_start(&calls, 4);
    assert_eq!(poor.len(), 1);
    assert_eq!(poor[0].calls, 0..6);
    assert_eq!(poor[0].reason, POOR_START_SIGNAL);
    // A read whose first calls settle at once has no poor start.
    assert!(poor_start(&calls[6..], 4).is_empty());
    // No run long enough: the whole read.
    assert_eq!(poor_start(&calls[..6], 4)[0].calls, 0..6);
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
        locus.channel_widths = std::array::from_fn(|index| if index == channel { 4 } else { 0 });
    }
    // Call 6 is as broad as a dye blob; call 7 reads the base of call 8,
    // its neighbour towards the read's end.
    signal.loci[6].channel_widths = [0, 0, 10, 0];
    quality.trim_end_0based_exclusive = 8;
    quality.retained_sequence = "A".repeat(8);
    let evidence = evidence_of(&calls, &signal, &callability, &quality, &[], &config(2))?;
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
        locus.channel_widths = std::array::from_fn(|index| if index == channel { 4 } else { 0 });
    }
    quality.trim_end_0based_exclusive = 10;
    quality.retained_sequence = "A".repeat(10);
    let lifted = |callability: &ReadCallability| -> Result<Vec<bool>> {
        let evidence = evidence_of(&calls, &signal, callability, &quality, &[], &config(4))?;
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

/// Sanger products of single-channel calls of `bases` at `heights`, each
/// call's own channel 4 wide.
fn based_products(
    bases: &[char],
    heights: &[i32],
) -> (
    BaseCalls,
    SignalAnalysis,
    ReadCallability,
    QualityControlResult,
) {
    let (mut calls, mut signal, callability, quality) = products(heights, &vec![40; heights.len()]);
    for ((call, &base), &height) in calls.calls.iter_mut().zip(bases).zip(heights) {
        let channel = Nucleotide::from_char(base).map_or(0, Nucleotide::channel_index);
        call.primary = base;
        call.ambiguity = base;
        call.qualifying_channels = vec![Nucleotide::ALL[channel]];
        for (index, peak) in call.peaks.iter_mut().enumerate() {
            peak.height = if index == channel { height } else { 0 };
        }
    }
    for (locus, &base) in signal.loci.iter_mut().zip(bases) {
        let channel = Nucleotide::from_char(base).map_or(0, Nucleotide::channel_index);
        locus.channel_widths = std::array::from_fn(|index| if index == channel { 4 } else { 0 });
    }
    (calls, signal, callability, quality)
}

// SRS-CALL-017
#[test]
fn judges_calls_before_a_closing_run_against_the_ends_own_norms() -> Result<()> {
    // Twenty calls; the trim ends at 16, so with a margin of 6 calls
    // 10-15 sit at the trailing edge. The read's signal decays: calls
    // 10-19 stand at half the read's median.
    let mut heights = [400; 20];
    heights[10..].fill(200);
    let lifted = |bases: &[char], fraction: f64, broad: Option<usize>| -> Result<Vec<bool>> {
        let (calls, mut signal, callability, mut quality) = based_products(bases, &heights);
        if let Some(broad) = broad {
            signal.loci[broad].channel_widths = [10; 4];
        }
        quality.trim_end_0based_exclusive = 16;
        quality.retained_sequence = "A".repeat(16);
        let config = SangerEvidenceConfig {
            end_norm_calls: 6,
            end_height_fraction: fraction,
            ..config(6)
        };
        let evidence = evidence_of(&calls, &signal, &callability, &quality, &[], &config)?;
        Ok(evidence.calls()[10..16]
            .iter()
            .map(|call| !call.identity.edge)
            .collect())
    };
    // The interval ends with the run `GG` (14-15): the calls before it
    // reach half the end's median, and `TT` (10-11) is read behind the
    // vouched `A` (12) that closes it; the run's own calls stay at the
    // edge.
    let closed = [
        'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'T', 'T', 'A', 'C', 'G', 'G', 'G', 'G',
        'T', 'A',
    ];
    assert_eq!(
        lifted(&closed, 0.5, None)?,
        [true, true, true, true, false, false]
    );
    // Below the end's height fraction, nothing is lifted.
    assert_eq!(lifted(&closed, 1.5, None)?, [false; 6]);
    // A run's bound (13) as broad as a dye blob holds back every call
    // before it; a broad call further in (11) holds back only itself
    // and the run behind it.
    assert_eq!(lifted(&closed, 0.5, Some(13))?, [false; 6]);
    assert_eq!(
        lifted(&closed, 0.5, Some(11))?,
        [false, false, true, true, false, false]
    );
    // An interval that no run closes keeps the read's norms, which these
    // decayed calls do not reach.
    let open = [
        'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T', 'A', 'C', 'T', 'G', 'A', 'C', 'G', 'T', 'G', 'G',
        'T', 'A',
    ];
    assert_eq!(lifted(&open, 0.5, None)?, [false; 6]);
    Ok(())
}

// SRS-CALL-017
#[test]
fn vouches_for_a_call_before_a_run_the_read_loses_phase_after_inside_the_read() -> Result<()> {
    // Twenty-four calls in the trim; the run `GGG` (10-12) is followed by
    // calls masked after the repeat run (13-16), and the read regains
    // phase at 17. With a margin of 4, calls 9-12 and 17-20 sit at an edge.
    let bases = [
        'A', 'C', 'G', 'T', 'A', 'C', 'A', 'T', 'C', 'A', 'G', 'G', 'G', 'T', 'T', 'A', 'C', 'A',
        'C', 'T', 'A', 'C', 'G', 'T',
    ];
    let (calls, signal, mut callability, quality) = based_products(&bases, &[200; 24]);
    callability.segments = vec![
        PhaseSegment {
            call_start_0based: 0,
            call_end_0based_exclusive: 13,
            state: PhaseState::InPhase,
            after_repeat: false,
            shadow: None,
        },
        PhaseSegment {
            call_start_0based: 13,
            call_end_0based_exclusive: 17,
            state: PhaseState::Mixed,
            after_repeat: true,
            shadow: None,
        },
        PhaseSegment {
            call_start_0based: 17,
            call_end_0based_exclusive: 24,
            state: PhaseState::InPhase,
            after_repeat: false,
            shadow: None,
        },
    ];
    callability.mask[13..17].fill(Some(PhaseState::Mixed));
    let evidence = evidence_of(&calls, &signal, &callability, &quality, &[], &config(4))?;
    let lifted = evidence.calls()[8..22]
        .iter()
        .map(|call| call.placement == Placement::Scored && !call.identity.edge)
        .collect::<Vec<_>>();
    // The call before the run is in phase; the run, the masked calls,
    // and the calls after the loss keep their verdicts.
    assert_eq!(
        lifted,
        [
            true, true, false, false, false, false, false, false, false, false, false, false,
            false, true
        ]
    );
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
        locus.channel_widths = std::array::from_fn(|index| if index == channel { 4 } else { 0 });
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
    let evidence = evidence_of(&calls, &signal, &callability, &quality, &[], &config(2))?;
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
    let evidence = evidence_of(&calls, &signal, &callability, &quality, &[], &config(3))?;
    assert_eq!(lifted(&evidence), [false, false, false]);
    // With call 3 reading `G`, call 4 qualifies, and the run chains on.
    let mut bases = bases;
    bases[3] = 'G';
    let (calls, signal, callability, quality) = start_products(bases, heights);
    let evidence = evidence_of(&calls, &signal, &callability, &quality, &[], &config(3))?;
    assert_eq!(lifted(&evidence), [true, true, true]);
    // A broad call 5 breaks the chain: call 6 no longer follows a vouched
    // call of its base.
    let (calls, mut signal, callability, quality) = start_products(bases, heights);
    signal.loci[5].channel_widths = [0, 0, 0, 10];
    let evidence = evidence_of(&calls, &signal, &callability, &quality, &[], &config(3))?;
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
    let evidence = evidence_of(&calls, &signal, &callability, &quality, &[], &config(0))?;
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
    let masked =
        |index: usize| tail && ((19..26).contains(&index) || (close && (17..19).contains(&index)));
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
        SIGNAL_SPIKE,
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
            "signal_spike",
            "post_homopolymer",
            "dephased_signal",
            "mixed_signal",
            "weak_signal",
            "irregular_spacing",
        ]
    );
}
