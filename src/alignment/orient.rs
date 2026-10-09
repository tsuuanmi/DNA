//! Forward/reverse profile alignment selection and circular coordinate projection.

use std::cmp::Ordering;

use crate::alignment::exact::{self, UpperBoundPlacement};
use crate::alignment::gotoh;
use crate::alignment::traceback::RawAlignment;
use crate::config::AlignmentConfig;
use crate::error::{AlignmentError, Result};
use crate::model::alignment::{Alignment, AlignmentColumn, Orientation, ReferenceSegment};
use crate::model::callability::{PhaseState, ReadCallability};
use crate::model::locus_evidence::EvidenceProfile;
use crate::model::nucleotide::reverse_complement;
use crate::model::quality::QualityControlResult;
use crate::model::reference::{Reference, ReferenceTopology};
use crate::model::signal::SignalAnalysis;

struct Candidate {
    orientation: Orientation,
    mapping: Vec<usize>,
    placements: Vec<RawAlignment>,
}

/// Aligns both evidence-profile orientations and returns one unique selected result.
///
/// Masked calls inside the trim interval enter the query as unresolved `N`
/// without an evidence profile, so they neither score nor count as callable;
/// dephased calls keep their call and profile, because they still read the
/// main ladder and anchor the alignment.
pub(crate) fn align_best(
    qc: &QualityControlResult,
    callability: &ReadCallability,
    signal: &SignalAnalysis,
    reference: &Reference,
    config: &AlignmentConfig,
) -> Result<Alignment> {
    if callability.mask.len() != qc.per_call.len() {
        return Err(AlignmentError::Inconsistent(
            "callability mask and quality evidence disagree in call count",
        )
        .into());
    }
    let masked = &callability.mask;
    let forward_query = qc
        .retained_sequence
        .chars()
        .zip(qc.trim_start_0based..)
        .map(|(base, index)| if unresolved(masked[index]) { 'N' } else { base })
        .collect::<String>();
    let reverse_query = reverse_complement(&forward_query);
    let forward_profiles = retained_profiles(qc, signal, masked)?;
    let reverse_profiles = reverse_profiles(&forward_profiles);
    let forward_mapping = (qc.trim_start_0based..qc.trim_end_0based_exclusive).collect::<Vec<_>>();
    let reverse_mapping = (qc.trim_start_0based..qc.trim_end_0based_exclusive)
        .rev()
        .collect::<Vec<_>>();
    let (working_reference, modulo_length) = match reference.topology {
        ReferenceTopology::Linear => (reference.sequence.clone(), None),
        ReferenceTopology::Circular => (
            format!("{}{}", reference.sequence, reference.sequence),
            Some(reference.len()),
        ),
    };
    let forward_proof = exact::classify(
        &forward_query,
        &forward_profiles,
        &working_reference,
        config,
        modulo_length,
    );
    let reverse_proof = exact::classify(
        &reverse_query,
        &reverse_profiles,
        &working_reference,
        config,
        modulo_length,
    );
    if let Some(selected) = select_proven_candidate(
        forward_proof,
        reverse_proof,
        &forward_mapping,
        &reverse_mapping,
    )? {
        return finish_alignment(&selected, reference, config, masked);
    }

    let mut forward_fast = exact::align_pruned(
        &forward_query,
        &forward_profiles,
        &working_reference,
        config,
        modulo_length,
    )?;
    let mut reverse_fast = exact::align_pruned(
        &reverse_query,
        &reverse_profiles,
        &working_reference,
        config,
        modulo_length,
    )?;

    if reverse_fast.is_none()
        && let Some(forward_placements) = forward_fast.as_ref()
    {
        let threshold = forward_placements[0].score;
        if let Some(reverse_placements) = exact::align_at_or_above(
            &reverse_query,
            &reverse_profiles,
            &working_reference,
            config,
            modulo_length,
            threshold,
        )? {
            if reverse_placements.is_empty() {
                let selected = Candidate {
                    orientation: Orientation::Forward,
                    mapping: forward_mapping,
                    placements: forward_fast.take().ok_or(AlignmentError::Inconsistent(
                        "missing proven forward placement",
                    ))?,
                };
                return finish_alignment(&selected, reference, config, masked);
            }
            reverse_fast = Some(reverse_placements);
        }
    }

    if forward_fast.is_none()
        && let Some(reverse_placements) = reverse_fast.as_ref()
    {
        let threshold = reverse_placements[0].score;
        if let Some(forward_placements) = exact::align_at_or_above(
            &forward_query,
            &forward_profiles,
            &working_reference,
            config,
            modulo_length,
            threshold,
        )? {
            if forward_placements.is_empty() {
                let selected = Candidate {
                    orientation: Orientation::Reverse,
                    mapping: reverse_mapping,
                    placements: reverse_fast.take().ok_or(AlignmentError::Inconsistent(
                        "missing proven reverse placement",
                    ))?,
                };
                return finish_alignment(&selected, reference, config, masked);
            }
            forward_fast = Some(forward_placements);
        }
    }

    let forward = Candidate {
        orientation: Orientation::Forward,
        mapping: forward_mapping,
        placements: match forward_fast {
            Some(placements) => placements,
            None => gotoh::align(
                &forward_query,
                &forward_profiles,
                &working_reference,
                config,
                modulo_length,
            )?,
        },
    };
    let reverse = Candidate {
        orientation: Orientation::Reverse,
        mapping: reverse_mapping,
        placements: match reverse_fast {
            Some(placements) => placements,
            None => gotoh::align(
                &reverse_query,
                &reverse_profiles,
                &working_reference,
                config,
                modulo_length,
            )?,
        },
    };
    let ordering = compare(&forward.placements[0], &reverse.placements[0]);
    let selected = match ordering {
        Ordering::Greater => &forward,
        Ordering::Less => &reverse,
        Ordering::Equal => {
            return Err(AlignmentError::OrientationTie.into());
        }
    };
    finish_alignment(selected, reference, config, masked)
}

fn select_proven_candidate(
    forward: UpperBoundPlacement,
    reverse: UpperBoundPlacement,
    forward_mapping: &[usize],
    reverse_mapping: &[usize],
) -> Result<Option<Candidate>> {
    use UpperBoundPlacement::{Ambiguous, Unattained, Unique, Unproven};

    match (forward, reverse) {
        (Unproven, _) | (_, Unproven) | (Unattained, Unattained) => Ok(None),
        (Unique(raw), Unattained) => Ok(Some(Candidate {
            orientation: Orientation::Forward,
            mapping: forward_mapping.to_vec(),
            placements: vec![raw],
        })),
        (Unattained, Unique(raw)) => Ok(Some(Candidate {
            orientation: Orientation::Reverse,
            mapping: reverse_mapping.to_vec(),
            placements: vec![raw],
        })),
        (Ambiguous, Unattained) | (Unattained, Ambiguous) => {
            Err(AlignmentError::AmbiguousPlacement.into())
        }
        (Unique(_) | Ambiguous, Unique(_) | Ambiguous) => {
            Err(AlignmentError::OrientationTie.into())
        }
    }
}

fn finish_alignment(
    selected: &Candidate,
    reference: &Reference,
    config: &AlignmentConfig,
    masked: &[Option<PhaseState>],
) -> Result<Alignment> {
    if selected.placements.len() != 1 {
        return Err(AlignmentError::AmbiguousPlacement.into());
    }
    let raw = &selected.placements[0];
    if raw.metrics.callable_columns < config.minimum_callable_bases {
        return Err(AlignmentError::TooFewCallableColumns {
            found: raw.metrics.callable_columns,
            minimum: config.minimum_callable_bases,
        }
        .into());
    }
    if raw.metrics.callable_identity < config.minimum_identity {
        return Err(AlignmentError::LowIdentity {
            identity: raw.metrics.callable_identity,
            minimum: config.minimum_identity,
        }
        .into());
    }
    let (segments, wraps_origin) = segments(raw, reference);
    let columns = raw
        .columns
        .iter()
        .map(|column| AlignmentColumn {
            query_base: column.query_base,
            reference_base: column.reference_base,
            original_call_index_0based: column
                .query_index
                .and_then(|index| selected.mapping.get(index).copied()),
            reference_index_0based: column
                .reference_index
                .map(|index| match reference.topology {
                    ReferenceTopology::Linear => index,
                    ReferenceTopology::Circular => index % reference.len(),
                }),
        })
        .collect::<Vec<_>>();
    let mask_of = |column: &AlignmentColumn| {
        column
            .original_call_index_0based
            .and_then(|index| masked.get(index).copied())
            .flatten()
    };
    let masked_query_bases = columns
        .iter()
        .filter(|column| mask_of(column).is_some())
        .count();
    let hidden = columns
        .iter()
        .filter(|column| unresolved(mask_of(column)))
        .count();
    let mut metrics = raw.metrics.clone();
    metrics.unresolved_query_bases =
        metrics
            .unresolved_query_bases
            .checked_sub(hidden)
            .ok_or(AlignmentError::Inconsistent(
                "masked query bases exceed unresolved query bases",
            ))?;
    metrics.masked_query_bases = masked_query_bases;
    let callable_segments = callable_segments(&columns, |column| mask_of(column).is_some());
    Ok(Alignment {
        orientation: selected.orientation,
        score: raw.score,
        reference_segments: segments,
        callable_segments,
        wraps_origin,
        metrics,
        columns,
    })
}

fn retained_profiles(
    qc: &QualityControlResult,
    signal: &SignalAnalysis,
    masked: &[Option<PhaseState>],
) -> Result<Vec<Option<EvidenceProfile>>> {
    if signal.loci.len() != qc.per_call.len() {
        return Err(AlignmentError::CallCountMismatch {
            loci: signal.loci.len(),
            qualities: qc.per_call.len(),
        }
        .into());
    }
    if qc.trim_start_0based > qc.trim_end_0based_exclusive
        || qc.trim_end_0based_exclusive > signal.loci.len()
    {
        return Err(AlignmentError::InvalidTrim {
            start: qc.trim_start_0based,
            end: qc.trim_end_0based_exclusive,
            profiles: signal.loci.len(),
        }
        .into());
    }
    let profiles = (qc.trim_start_0based..qc.trim_end_0based_exclusive)
        .map(|index| {
            if unresolved(masked[index]) {
                None
            } else {
                signal.loci[index].profile
            }
        })
        .collect::<Vec<_>>();
    if profiles.len() != qc.retained_sequence.len() {
        return Err(AlignmentError::RetainedLengthMismatch {
            bases: qc.retained_sequence.len(),
            profiles: profiles.len(),
        }
        .into());
    }
    Ok(profiles)
}

/// Reference segments observed by unmasked calls: a call column counts when
/// its call is not masked, a deletion column when the call columns on both
/// sides do; runs of consecutive reference indexes form one segment.
fn callable_segments(
    columns: &[AlignmentColumn],
    masked: impl Fn(&AlignmentColumn) -> bool,
) -> Vec<ReferenceSegment> {
    let calls = columns
        .iter()
        .enumerate()
        .filter(|(_, column)| column.original_call_index_0based.is_some())
        .map(|(index, column)| (index, !masked(column)))
        .collect::<Vec<_>>();
    let callable_call_around = |index: usize| {
        let before = calls.iter().rev().find(|(position, _)| *position < index);
        let after = calls.iter().find(|(position, _)| *position > index);
        before.is_some_and(|(_, callable)| *callable)
            && after.is_some_and(|(_, callable)| *callable)
    };
    let mut segments: Vec<ReferenceSegment> = Vec::new();
    for (index, column) in columns.iter().enumerate() {
        let Some(reference) = column.reference_index_0based else {
            continue;
        };
        let callable = if column.original_call_index_0based.is_some() {
            !masked(column)
        } else {
            callable_call_around(index)
        };
        if !callable {
            continue;
        }
        match segments.last_mut() {
            Some(last) if last.end_0based_exclusive == reference => {
                last.end_0based_exclusive += 1;
            }
            _ => segments.push(ReferenceSegment {
                start_0based: reference,
                end_0based_exclusive: reference + 1,
            }),
        }
    }
    segments
}

/// A masked call aligns as unresolved unless it is dephased.
fn unresolved(mask: Option<PhaseState>) -> bool {
    mask.is_some_and(|state| state != PhaseState::Dephased)
}

fn reverse_profiles(profiles: &[Option<EvidenceProfile>]) -> Vec<Option<EvidenceProfile>> {
    profiles
        .iter()
        .rev()
        .map(|profile| profile.map(EvidenceProfile::complemented))
        .collect()
}

fn compare(left: &RawAlignment, right: &RawAlignment) -> Ordering {
    left.score.cmp(&right.score)
}

fn segments(alignment: &RawAlignment, reference: &Reference) -> (Vec<ReferenceSegment>, bool) {
    match reference.topology {
        ReferenceTopology::Linear => (
            vec![ReferenceSegment {
                start_0based: alignment.start_reference,
                end_0based_exclusive: alignment.end_reference,
            }],
            false,
        ),
        ReferenceTopology::Circular => {
            let length = reference.len();
            let start = alignment.start_reference % length;
            let span = alignment.end_reference - alignment.start_reference;
            let unwrapped_end = start + span;
            if unwrapped_end <= length {
                (
                    vec![ReferenceSegment {
                        start_0based: start,
                        end_0based_exclusive: unwrapped_end,
                    }],
                    false,
                )
            } else {
                (
                    vec![
                        ReferenceSegment {
                            start_0based: start,
                            end_0based_exclusive: length,
                        },
                        ReferenceSegment {
                            start_0based: 0,
                            end_0based_exclusive: unwrapped_end - length,
                        },
                    ],
                    true,
                )
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::model::alignment::AlignmentMetrics;
    use crate::model::locus_evidence::LocusEvidence;
    use crate::model::quality::CallQuality;
    use crate::model::reference::ReferenceTopology;
    use crate::model::signal::SangerIntegrity;

    use super::*;

    fn raw(score: i64, exact_matches: usize, mismatches: usize, gap_opens: usize) -> RawAlignment {
        RawAlignment {
            score,
            start_reference: 0,
            end_reference: 1,
            columns: Vec::new(),
            metrics: AlignmentMetrics {
                exact_matches,
                mismatches,
                gap_opens,
                callable_columns: exact_matches + mismatches,
                callable_identity: 0.0,
                unresolved_query_bases: 0,
                masked_query_bases: 0,
            },
        }
    }

    #[test]
    fn reverse_profiles_reverse_order_and_complement_channels() {
        let profiles = vec![
            Some(EvidenceProfile {
                weights: [1.0, 0.0, 0.0, 0.0],
            }),
            None,
            Some(EvidenceProfile {
                weights: [0.0, 1.0, 0.0, 0.0],
            }),
        ];
        let reversed = reverse_profiles(&profiles);
        assert_eq!(
            reversed[0].map(|profile| profile.weights),
            Some([0.0, 0.0, 1.0, 0.0])
        );
        assert!(reversed[1].is_none());
        assert_eq!(
            reversed[2].map(|profile| profile.weights),
            Some([0.0, 0.0, 0.0, 1.0])
        );
    }

    fn qc(sequence: &str) -> QualityControlResult {
        QualityControlResult {
            per_call: sequence
                .chars()
                .enumerate()
                .map(|(index_0based, _)| CallQuality {
                    index_0based,
                    penalty: 0,
                    relative_quality_score: 60,
                    vendor_quality_applies: false,
                })
                .collect(),
            trim_start_0based: 0,
            trim_end_0based_exclusive: sequence.len(),
            retained_sequence: sequence.into(),
        }
    }

    fn signal(sequence: &str) -> SignalAnalysis {
        let loci = sequence
            .bytes()
            .enumerate()
            .map(|(call_index_0based, base)| {
                let weights = match base {
                    b'A' => [1.0, 0.0, 0.0, 0.0],
                    b'C' => [0.0, 1.0, 0.0, 0.0],
                    b'G' => [0.0, 0.0, 1.0, 0.0],
                    b'T' => [0.0, 0.0, 0.0, 1.0],
                    _ => [0.0; 4],
                };
                LocusEvidence {
                    call_index_0based,
                    locus_position_0based: call_index_0based,
                    window_start_0based: call_index_0based,
                    window_end_0based_exclusive: call_index_0based + 1,
                    context_call_start_0based: call_index_0based,
                    context_call_end_0based_exclusive: call_index_0based + 1,
                    context_sample_start_0based: call_index_0based,
                    context_sample_end_0based_exclusive: call_index_0based + 1,
                    event_position_0based: call_index_0based,
                    channel_heights: [0; 4],
                    channel_baselines: [0.0; 4],
                    channel_noise_sigmas: [1.0; 4],
                    corrected_amplitudes: weights,
                    snrs: weights,
                    profile: Some(EvidenceProfile { weights }),
                }
            })
            .collect();
        SignalAnalysis {
            integrity: SangerIntegrity {
                locus_count: sequence.len(),
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
        }
    }

    fn config() -> AlignmentConfig {
        AlignmentConfig {
            match_score: 3,
            mismatch_score: -5,
            ambiguous_score: 0,
            gap_open_score: -10,
            gap_extension_score: -4,
            minimum_callable_bases: 1,
            minimum_identity: 0.8,
        }
    }

    fn reference() -> Reference {
        Reference {
            name: "ref".into(),
            sequence: "GCCAAAAGTT".into(),
            topology: ReferenceTopology::Linear,
            sequence_sha256: String::new(),
        }
    }

    #[test]
    fn forward_and_reverse_reads_share_canonical_deletion_coordinate() -> Result<()> {
        let forward = align_best(
            &qc("GCCAAAGTT"),
            &ReadCallability::in_phase(9),
            &signal("GCCAAAGTT"),
            &reference(),
            &config(),
        )?;
        let reverse = align_best(
            &qc("AACTTTGGC"),
            &ReadCallability::in_phase(9),
            &signal("AACTTTGGC"),
            &reference(),
            &config(),
        )?;

        assert_eq!(forward.orientation, Orientation::Forward);
        assert_eq!(reverse.orientation, Orientation::Reverse);

        let forward_deleted = forward
            .columns
            .iter()
            .find(|column| column.query_base == '-')
            .and_then(|column| column.reference_index_0based);
        let reverse_deleted = reverse
            .columns
            .iter()
            .find(|column| column.query_base == '-')
            .and_then(|column| column.reference_index_0based);
        assert_eq!(forward_deleted, Some(6));
        assert_eq!(reverse_deleted, forward_deleted);
        Ok(())
    }

    #[test]
    fn upper_bound_orientation_tie_preserves_existing_error() {
        let reference = Reference {
            name: "ref".into(),
            sequence: "ACGT".into(),
            topology: ReferenceTopology::Linear,
            sequence_sha256: String::new(),
        };
        let error = align_best(
            &qc("ACGT"),
            &ReadCallability::in_phase(4),
            &signal("ACGT"),
            &reference,
            &config(),
        )
        .err()
        .map(|error| error.to_string());
        assert_eq!(
            error.as_deref(),
            Some("alignment failed: forward and reverse evidence-profile scores are tied")
        );
    }

    #[test]
    fn repeated_upper_bound_placement_preserves_existing_error() {
        let reference = Reference {
            name: "ref".into(),
            sequence: "AAAAA".into(),
            topology: ReferenceTopology::Linear,
            sequence_sha256: String::new(),
        };
        let error = align_best(
            &qc("AAA"),
            &ReadCallability::in_phase(3),
            &signal("AAA"),
            &reference,
            &config(),
        )
        .err()
        .map(|error| error.to_string());
        assert_eq!(
            error.as_deref(),
            Some("alignment failed: selected orientation has multiple equally scoring placements")
        );
    }

    #[test]
    fn callable_segments_skip_masked_calls_and_deletions_beside_them() {
        let column = |call: Option<usize>, reference: Option<usize>| AlignmentColumn {
            query_base: if call.is_some() { 'A' } else { '-' },
            reference_base: if reference.is_some() { 'A' } else { '-' },
            original_call_index_0based: call,
            reference_index_0based: reference,
        };
        let columns = [
            column(Some(0), Some(10)),
            column(None, Some(11)),
            column(Some(1), Some(12)),
            column(Some(2), Some(13)),
            column(None, Some(14)),
            column(Some(3), Some(15)),
            column(Some(4), None),
            column(Some(5), Some(16)),
        ];
        let masked = |column: &AlignmentColumn| column.original_call_index_0based == Some(3);
        assert_eq!(
            callable_segments(&columns, masked),
            [
                ReferenceSegment {
                    start_0based: 10,
                    end_0based_exclusive: 14,
                },
                ReferenceSegment {
                    start_0based: 16,
                    end_0based_exclusive: 17,
                },
            ]
        );
    }

    #[test]
    fn orientation_comparison_uses_profile_score_only() {
        let left = raw(100, 1, 9, 5);
        let right = raw(100, 10, 0, 0);
        assert_eq!(compare(&left, &right), Ordering::Equal);
        assert_eq!(compare(&raw(101, 0, 10, 10), &right), Ordering::Greater);
    }
}
