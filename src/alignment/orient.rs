//! Forward/reverse profile alignment selection and circular coordinate projection.

use std::cmp::Ordering;

use crate::alignment::gotoh;
use crate::alignment::traceback::RawAlignment;
use crate::config::AlignmentConfig;
use crate::error::{Error, Result};
use crate::model::alignment::{Alignment, AlignmentColumn, Orientation, ReferenceSegment};
use crate::model::locus_evidence::EvidenceProfile;
use crate::model::nucleotide::reverse_complement;
use crate::model::quality::QualityControlResult;
use crate::model::reference::{Reference, ReferenceTopology};
use crate::model::signal::DNAAnalysis;

struct Candidate {
    orientation: Orientation,
    mapping: Vec<usize>,
    placements: Vec<RawAlignment>,
}

/// Aligns both evidence-profile orientations and returns one unique selected result.
pub(crate) fn align_best(
    qc: &QualityControlResult,
    signal: &DNAAnalysis,
    reference: &Reference,
    config: &AlignmentConfig,
) -> Result<Alignment> {
    let forward_query = qc.retained_sequence.clone();
    let reverse_query = reverse_complement(&forward_query);
    let forward_profiles = retained_profiles(qc, signal)?;
    let reverse_profiles = reverse_profiles(&forward_profiles);
    let forward_mapping = (qc.trim_start_0based..qc.trim_end_0based_exclusive).collect();
    let reverse_mapping = (qc.trim_start_0based..qc.trim_end_0based_exclusive)
        .rev()
        .collect();
    let (working_reference, modulo_length) = match reference.topology {
        ReferenceTopology::Linear => (reference.sequence.clone(), None),
        ReferenceTopology::Circular => (
            format!("{}{}", reference.sequence, reference.sequence),
            Some(reference.len()),
        ),
    };
    let forward = Candidate {
        orientation: Orientation::Forward,
        mapping: forward_mapping,
        placements: gotoh::align(
            &forward_query,
            &forward_profiles,
            &working_reference,
            config,
            modulo_length,
        )?,
    };
    let reverse = Candidate {
        orientation: Orientation::Reverse,
        mapping: reverse_mapping,
        placements: gotoh::align(
            &reverse_query,
            &reverse_profiles,
            &working_reference,
            config,
            modulo_length,
        )?,
    };
    let ordering = compare(&forward.placements[0], &reverse.placements[0]);
    let selected = match ordering {
        Ordering::Greater => &forward,
        Ordering::Less => &reverse,
        Ordering::Equal => {
            return Err(Error::Alignment(
                "forward and reverse evidence-profile scores are tied".into(),
            ));
        }
    };
    if selected.placements.len() != 1 {
        return Err(Error::Alignment(
            "selected orientation has multiple equally scoring placements".into(),
        ));
    }
    let raw = &selected.placements[0];
    if raw.metrics.callable_columns < config.minimum_callable_bases {
        return Err(Error::Alignment(format!(
            "alignment has {} callable columns; minimum is {}",
            raw.metrics.callable_columns, config.minimum_callable_bases
        )));
    }
    if raw.metrics.callable_identity < config.minimum_identity {
        return Err(Error::Alignment(format!(
            "alignment callable identity {:.4} is below {:.4}",
            raw.metrics.callable_identity, config.minimum_identity
        )));
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
        .collect();
    Ok(Alignment {
        orientation: selected.orientation,
        score: raw.score,
        reference_segments: segments,
        wraps_origin,
        metrics: raw.metrics.clone(),
        columns,
    })
}

fn retained_profiles(
    qc: &QualityControlResult,
    signal: &DNAAnalysis,
) -> Result<Vec<Option<EvidenceProfile>>> {
    if signal.loci.len() != qc.per_call.len() {
        return Err(Error::Alignment(format!(
            "signal/quality call count mismatch: {} loci, {} quality records",
            signal.loci.len(),
            qc.per_call.len()
        )));
    }
    if qc.trim_start_0based > qc.trim_end_0based_exclusive
        || qc.trim_end_0based_exclusive > signal.loci.len()
    {
        return Err(Error::Alignment(format!(
            "invalid trim interval {}..{} for {} locus profiles",
            qc.trim_start_0based,
            qc.trim_end_0based_exclusive,
            signal.loci.len()
        )));
    }
    let profiles = signal.loci[qc.trim_start_0based..qc.trim_end_0based_exclusive]
        .iter()
        .map(|locus| locus.profile)
        .collect::<Vec<_>>();
    if profiles.len() != qc.retained_sequence.len() {
        return Err(Error::Alignment(format!(
            "retained sequence/profile length mismatch: {} bases, {} profiles",
            qc.retained_sequence.len(),
            profiles.len()
        )));
    }
    Ok(profiles)
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
    use crate::model::signal::TraceIntegrity;

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

    fn signal(sequence: &str) -> DNAAnalysis {
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
                    ploc_0based: call_index_0based,
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
        DNAAnalysis {
            integrity: TraceIntegrity {
                ploc_count: sequence.len(),
                vendor_primary_count: None,
                vendor_quality_count: None,
                minimum_ploc_spacing: None,
                median_ploc_spacing: None,
                maximum_ploc_spacing: None,
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
            &signal("GCCAAAGTT"),
            &reference(),
            &config(),
        )?;
        let reverse = align_best(
            &qc("AACTTTGGC"),
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
    fn orientation_comparison_uses_profile_score_only() {
        let left = raw(100, 1, 9, 5);
        let right = raw(100, 10, 0, 0);
        assert_eq!(compare(&left, &right), Ordering::Equal);
        assert_eq!(compare(&raw(101, 0, 10, 10), &right), Ordering::Greater);
    }
}
