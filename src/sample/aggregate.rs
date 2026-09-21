//! Sample-level validation, deterministic read ordering, and evidence assembly.

use std::collections::BTreeSet;

use crate::config::SampleReconciliationConfig;
use crate::error::{Error, Result};
use crate::model::read_observation::ReadObservation;
use crate::model::sample_evidence::{
    SampleEvidence, SampleReadAlignmentEvidence, SampleReadEvidence,
};

use super::{coverage, loci, overlap, variants};

/// Aggregates independently processed reads without using filenames or pair labels as merge keys.
pub(crate) fn aggregate(
    reads: &[ReadObservation],
    config: &SampleReconciliationConfig,
) -> Result<SampleEvidence> {
    let ordered = validated_ordered_reads(reads)?;
    let first = ordered
        .first()
        .ok_or_else(|| Error::Sample("at least one read observation is required".into()))?;
    let reference_sha256 = first.reference_sha256.clone();
    let configuration_sha256 = first.configuration_sha256.clone();

    let read_evidence: Vec<SampleReadEvidence> = ordered
        .iter()
        .map(|read| SampleReadEvidence {
            input_name: read.input_name.clone(),
            input_sha256: read.input_sha256.clone(),
            integrity: read.signal.integrity.clone(),
            alignment: SampleReadAlignmentEvidence {
                orientation: read.alignment.orientation,
                callable_bases: read.alignment.metrics.callable_columns,
                identity: read.alignment.metrics.callable_identity,
                gap_opens: read.alignment.metrics.gap_opens,
                unresolved_bases: read.alignment.metrics.unresolved_query_bases,
                reference_segments: read.alignment.reference_segments.clone(),
                wraps_origin: read.alignment.wraps_origin,
            },
        })
        .collect();

    let coverage = coverage::summarize(&read_evidence)?;

    Ok(SampleEvidence {
        reference_sha256,
        configuration_sha256,
        reads: read_evidence,
        coverage,
        overlaps: overlap::assess(&ordered, config)?,
        locus_differences: loci::aggregate(&ordered)?,
        variants: variants::aggregate(&ordered)?,
    })
}

pub(crate) fn validated_ordered_reads(reads: &[ReadObservation]) -> Result<Vec<&ReadObservation>> {
    let first = reads
        .first()
        .ok_or_else(|| Error::Sample("at least one read observation is required".into()))?;

    let mut identities = BTreeSet::new();
    for read in reads {
        if read.reference_sha256 != first.reference_sha256 {
            return Err(Error::Sample(
                "all reads must use the same reference identity".into(),
            ));
        }
        if read.configuration_sha256 != first.configuration_sha256 {
            return Err(Error::Sample(
                "all reads must use the same scientific configuration identity".into(),
            ));
        }
        if !identities.insert(read.input_sha256.as_str()) {
            return Err(Error::Sample(
                "duplicate input trace content cannot contribute twice".into(),
            ));
        }
    }

    let mut ordered = reads.iter().collect::<Vec<_>>();
    ordered.sort_by(|left, right| left.input_sha256.cmp(&right.input_sha256));
    Ok(ordered)
}

#[cfg(test)]
mod tests {
    use crate::model::alignment::{
        Alignment, AlignmentColumn, AlignmentMetrics, Orientation, ReferenceSegment,
    };
    use crate::model::basecalls::{
        BaseCall, BaseCalls, ChannelPeak, PeakSource, PrimaryPeakEvidence,
    };
    use crate::model::locus_evidence::{EvidenceProfile, LocusEvidence};
    use crate::model::nucleotide::Nucleotide;
    use crate::model::quality::{CallQuality, QualityControlResult};
    use crate::model::signal::{DNAAnalysis, NoisyRegion};
    use crate::model::variant::{
        ObservedVariant, Variant, VariantCallMapping, VariantCallRole, VariantCallingResult,
        VariantExclusionReason, VariantKind,
    };

    use super::*;

    fn sample_config() -> SampleReconciliationConfig {
        SampleReconciliationConfig {
            minimum_comparable_bases: 1,
            minimum_overlap_agreement: 0.5,
        }
    }

    fn locus(index_0based: usize) -> LocusEvidence {
        LocusEvidence {
            call_index_0based: index_0based,
            ploc_0based: index_0based * 10,
            window_start_0based: index_0based * 10,
            window_end_0based_exclusive: index_0based * 10 + 1,
            context_call_start_0based: index_0based,
            context_call_end_0based_exclusive: index_0based + 1,
            context_sample_start_0based: index_0based * 10,
            context_sample_end_0based_exclusive: index_0based * 10 + 1,
            event_position_0based: index_0based * 10,
            channel_heights: [10, 20, 30, 40],
            channel_baselines: [0.0; 4],
            channel_noise_sigmas: [1.0; 4],
            corrected_amplitudes: [1.0, 2.0, 3.0, 4.0],
            snrs: [1.0, 2.0, 3.0, 4.0],
            profile: Some(EvidenceProfile {
                weights: [0.1, 0.2, 0.3, 0.4],
            }),
        }
    }

    fn observation(
        input_sha256: &str,
        reference_sha256: &str,
        configuration_sha256: &str,
        orientation: Orientation,
        columns: Vec<AlignmentColumn>,
        variants: Vec<Variant>,
    ) -> ReadObservation {
        let call_count = columns
            .iter()
            .filter_map(|column| column.original_call_index_0based)
            .max()
            .map_or(0, |index| index + 1);
        let per_call = (0..call_count)
            .map(|index_0based| CallQuality {
                index_0based,
                penalty: 0,
                relative_quality_score: 50,
                vendor_quality_applies: false,
            })
            .collect();
        let observed = variants
            .iter()
            .cloned()
            .map(|variant| ObservedVariant {
                variant,
                exclusion_reasons: Vec::new(),
            })
            .collect();
        let reference_segments = columns
            .iter()
            .filter_map(|column| column.reference_index_0based)
            .fold(None::<(usize, usize)>, |bounds, index| {
                Some(match bounds {
                    Some((start, end)) => (start.min(index), end.max(index + 1)),
                    None => (index, index + 1),
                })
            })
            .map_or_else(Vec::new, |(start_0based, end_0based_exclusive)| {
                vec![ReferenceSegment {
                    start_0based,
                    end_0based_exclusive,
                }]
            });
        ReadObservation {
            input_name: format!("{input_sha256}.ab1"),
            input_sha256: input_sha256.into(),
            reference_sha256: reference_sha256.into(),
            configuration_sha256: configuration_sha256.into(),
            calls: BaseCalls {
                calls: (0..call_count)
                    .map(|index_0based| BaseCall {
                        index_0based,
                        ploc_0based: index_0based * 10,
                        window_start_0based: index_0based * 10,
                        window_end_0based_exclusive: index_0based * 10 + 1,
                        peaks: Nucleotide::ALL.map(|base| ChannelPeak {
                            base,
                            height: 100,
                            position_0based: index_0based * 10,
                            source: PeakSource::LocalMaximum,
                        }),
                        primary_peak_evidence: Some(PrimaryPeakEvidence {
                            position_0based: index_0based * 10,
                            channel_heights: [10, 20, 100, 30],
                        }),
                        primary: 'G',
                        ambiguity: 'G',
                        qualifying_channels: vec![Nucleotide::G],
                        vendor_agrees: None,
                    })
                    .collect(),
                primary_sequence: "G".repeat(call_count),
            },
            signal: DNAAnalysis {
                integrity: crate::model::signal::TraceIntegrity {
                    ploc_count: call_count,
                    vendor_primary_count: None,
                    vendor_quality_count: None,
                    minimum_ploc_spacing: None,
                    median_ploc_spacing: None,
                    maximum_ploc_spacing: None,
                    clipped_channel_samples: 0,
                    maximum_to_median_event_signal_ratio: None,
                },
                loci: (0..call_count).map(locus).collect(),
                windows: Vec::new(),
                noisy_regions: Vec::new(),
            },
            quality: QualityControlResult {
                per_call,
                trim_start_0based: 0,
                trim_end_0based_exclusive: call_count,
                retained_sequence: String::new(),
            },
            alignment: Alignment {
                orientation,
                score: 1,
                reference_segments,
                wraps_origin: false,
                metrics: AlignmentMetrics {
                    exact_matches: 0,
                    mismatches: 0,
                    gap_opens: 0,
                    callable_columns: 0,
                    callable_identity: 0.0,
                    unresolved_query_bases: 0,
                },
                columns,
            },
            variants: VariantCallingResult {
                reported: variants,
                observed,
                excluded: Vec::new(),
            },
        }
    }

    fn column(
        query_base: char,
        reference_base: char,
        call: Option<usize>,
        reference: usize,
    ) -> AlignmentColumn {
        AlignmentColumn {
            query_base,
            reference_base,
            original_call_index_0based: call,
            reference_index_0based: Some(reference),
        }
    }

    fn snv(position_1based: usize, reference: &str, alternate: &str) -> Variant {
        Variant {
            contig: "reference".into(),
            position_1based,
            reference: reference.into(),
            alternate: alternate.into(),
            kind: VariantKind::Snv,
            calls: vec![VariantCallMapping {
                role: VariantCallRole::Supporting,
                call_index_0based: 0,
                reference_position_0based: Some(position_1based - 1),
            }],
        }
    }

    #[test]
    fn orders_reads_once_and_factors_read_identity_from_evidence() -> Result<()> {
        let mut forward = observation(
            "a",
            "reference",
            "config",
            Orientation::Forward,
            vec![column('G', 'A', Some(0), 72)],
            vec![snv(73, "A", "G")],
        );
        forward.signal.noisy_regions.push(NoisyRegion {
            call_start_0based: 0,
            call_end_0based_exclusive: 1,
            sample_start_0based: 0,
            sample_end_0based_exclusive: 10,
            minimum_primary_snr: 1.0,
        });
        let reverse = observation(
            "b",
            "reference",
            "config",
            Orientation::Reverse,
            vec![column('G', 'A', Some(0), 72)],
            vec![snv(73, "A", "G")],
        );

        let evidence = aggregate(&[reverse, forward], &sample_config())?;

        assert_eq!(evidence.reads[0].input_name, "a.ab1");
        assert_eq!(evidence.reads[1].input_name, "b.ab1");
        assert_eq!(evidence.coverage.len(), 1);
        assert_eq!(evidence.coverage[0].start_0based, 72);
        assert_eq!(evidence.coverage[0].end_0based_exclusive, 73);
        assert_eq!(evidence.coverage[0].read_depth, 2);
        assert_eq!(evidence.coverage[0].forward_depth, 1);
        assert_eq!(evidence.coverage[0].reverse_depth, 1);
        assert_eq!(evidence.overlaps.len(), 1);
        assert_eq!(evidence.overlaps[0].left_read_index, 0);
        assert_eq!(evidence.overlaps[0].right_read_index, 1);
        assert_eq!(evidence.overlaps[0].comparable_bases, 1);
        assert!(evidence.overlaps[0].eligible);
        assert_eq!(evidence.locus_differences.len(), 1);
        assert_eq!(evidence.locus_differences[0].support_topology.reads, 2);
        assert_eq!(
            evidence.locus_differences[0].support_topology.forward_reads,
            1
        );
        assert_eq!(
            evidence.locus_differences[0].support_topology.reverse_reads,
            1
        );
        assert_eq!(
            evidence.locus_differences[0]
                .support_topology
                .alternate_reads,
            2
        );
        assert_eq!(
            evidence.locus_differences[0]
                .support_topology
                .reference_reads,
            0
        );
        assert_eq!(
            evidence.locus_differences[0].support_topology.profile_reads,
            2
        );
        assert_eq!(
            evidence.locus_differences[0]
                .support_topology
                .profile_forward_reads,
            1
        );
        assert_eq!(
            evidence.locus_differences[0]
                .support_topology
                .profile_reverse_reads,
            1
        );
        let nucleotide_support = evidence.locus_differences[0].nucleotide_support;
        assert_eq!(nucleotide_support.contributors, 2);
        assert_eq!(nucleotide_support.forward_contributors, 1);
        assert_eq!(nucleotide_support.reverse_contributors, 1);
        assert_eq!(nucleotide_support.support, [0.5, 0.5, 0.5, 0.5]);
        assert_eq!(nucleotide_support.forward_support, [0.1, 0.2, 0.3, 0.4]);
        assert_eq!(nucleotide_support.reverse_support, [0.4, 0.3, 0.2, 0.1]);
        assert_eq!(
            nucleotide_support
                .mean_profile
                .map(|profile| profile.weights),
            Some([0.25, 0.25, 0.25, 0.25])
        );
        assert_eq!(
            nucleotide_support
                .forward_mean_profile
                .map(|profile| profile.weights),
            Some([0.1, 0.2, 0.3, 0.4])
        );
        assert_eq!(
            nucleotide_support
                .reverse_mean_profile
                .map(|profile| profile.weights),
            Some([0.4, 0.3, 0.2, 0.1])
        );
        let heterogeneity = nucleotide_support
            .heterogeneity
            .ok_or_else(|| Error::Sample("total profile heterogeneity is missing".into()))?;
        assert!((heterogeneity.within_profile_impurity - 0.7).abs() < 1e-12);
        assert!((heterogeneity.between_profile_dispersion - 0.05).abs() < 1e-12);
        assert!((heterogeneity.total_profile_heterogeneity - 0.75).abs() < 1e-12);
        assert_eq!(
            nucleotide_support
                .forward_heterogeneity
                .map(|geometry| geometry.between_profile_dispersion),
            Some(0.0)
        );
        assert_eq!(
            nucleotide_support
                .reverse_heterogeneity
                .map(|geometry| geometry.between_profile_dispersion),
            Some(0.0)
        );
        assert!(
            (nucleotide_support
                .directional_profile_distance
                .ok_or_else(|| Error::Sample("directional profile distance is missing".into()))?
                - 0.4)
                .abs()
                < 1e-12
        );
        assert_eq!(evidence.locus_differences[0].observations.len(), 2);
        assert_eq!(evidence.locus_differences[0].observations[0].read_index, 0);
        let forward_signal = evidence.locus_differences[0].observations[0]
            .signal
            .as_ref()
            .ok_or_else(|| Error::Sample("forward locus signal is missing".into()))?;
        assert_eq!(
            forward_signal.profile.map(|profile| profile.weights),
            Some([0.1, 0.2, 0.3, 0.4])
        );
        assert_eq!(forward_signal.corrected_amplitudes, [1.0, 2.0, 3.0, 4.0]);
        assert_eq!(forward_signal.snrs, [1.0, 2.0, 3.0, 4.0]);
        assert!(forward_signal.in_noisy_region);
        assert_eq!(
            evidence.locus_differences[0].observations[0].nucleotide_contribution,
            crate::model::sample_evidence::NucleotideContribution::Eligible
        );
        assert!(forward_signal.in_noisy_region);
        assert_eq!(evidence.locus_differences[0].observations[1].read_index, 1);
        let reverse_signal = evidence.locus_differences[0].observations[1]
            .signal
            .as_ref()
            .ok_or_else(|| Error::Sample("reverse locus signal is missing".into()))?;
        assert_eq!(
            reverse_signal.profile.map(|profile| profile.weights),
            Some([0.4, 0.3, 0.2, 0.1])
        );
        assert_eq!(reverse_signal.corrected_amplitudes, [4.0, 3.0, 2.0, 1.0]);
        assert_eq!(reverse_signal.snrs, [4.0, 3.0, 2.0, 1.0]);
        assert!(!reverse_signal.in_noisy_region);
        assert_eq!(evidence.variants.len(), 1);
        assert_eq!(evidence.variants[0].support_topology.reads, 2);
        assert_eq!(evidence.variants[0].support_topology.eligible_reads, 2);
        assert_eq!(evidence.variants[0].support_topology.forward_reads, 1);
        assert_eq!(evidence.variants[0].support_topology.reverse_reads, 1);
        assert_eq!(
            evidence.variants[0].support_topology.eligible_forward_reads,
            1
        );
        assert_eq!(
            evidence.variants[0].support_topology.eligible_reverse_reads,
            1
        );
        assert_eq!(evidence.variants[0].support[0].read_index, 0);
        assert_eq!(
            evidence.variants[0].support[0].calls[0]
                .signal
                .profile
                .map(|profile| profile.weights),
            Some([0.1, 0.2, 0.3, 0.4])
        );
        assert_eq!(
            evidence.variants[0].support[0].calls[0]
                .signal
                .corrected_amplitudes,
            [1.0, 2.0, 3.0, 4.0]
        );
        assert!(
            evidence.variants[0].support[0].calls[0]
                .signal
                .in_noisy_region
        );
        assert_eq!(evidence.variants[0].support[1].read_index, 1);
        assert_eq!(
            evidence.variants[0].support[1].calls[0]
                .signal
                .profile
                .map(|profile| profile.weights),
            Some([0.4, 0.3, 0.2, 0.1])
        );
        assert_eq!(
            evidence.variants[0].support[1].calls[0]
                .signal
                .corrected_amplitudes,
            [4.0, 3.0, 2.0, 1.0]
        );
        assert!(
            !evidence.variants[0].support[1].calls[0]
                .signal
                .in_noisy_region
        );
        Ok(())
    }

    #[test]
    fn omits_reference_matches_but_preserves_non_reference_states() -> Result<()> {
        let mut read = observation(
            "a",
            "reference",
            "config",
            Orientation::Forward,
            vec![
                column('A', 'A', Some(0), 10),
                column('G', 'C', Some(1), 11),
                column('N', 'T', Some(2), 12),
                column('-', 'G', None, 13),
                column('C', 'C', Some(3), 14),
            ],
            Vec::new(),
        );
        read.signal.loci[2].profile = None;

        let evidence = aggregate(&[read], &sample_config())?;

        assert_eq!(evidence.locus_differences.len(), 3);
        assert_eq!(evidence.locus_differences[0].position_1based, 12);
        assert_eq!(evidence.locus_differences[0].support_topology.reads, 1);
        assert_eq!(
            evidence.locus_differences[0]
                .support_topology
                .alternate_reads,
            1
        );
        assert_eq!(
            evidence.locus_differences[1]
                .support_topology
                .unresolved_reads,
            1
        );
        assert_eq!(
            evidence.locus_differences[2]
                .support_topology
                .deletion_reads,
            1
        );
        assert_eq!(
            evidence.locus_differences[0].support_topology.profile_reads,
            1
        );
        assert_eq!(
            evidence.locus_differences[1].support_topology.profile_reads,
            0
        );
        assert_eq!(
            evidence.locus_differences[2].support_topology.profile_reads,
            0
        );
        assert_eq!(
            evidence.locus_differences[0]
                .nucleotide_support
                .contributors,
            1
        );
        assert_eq!(
            evidence.locus_differences[1]
                .nucleotide_support
                .contributors,
            0
        );
        assert_eq!(
            evidence.locus_differences[1].nucleotide_support.support,
            [0.0; 4]
        );
        assert_eq!(
            evidence.locus_differences[2]
                .nucleotide_support
                .contributors,
            0
        );
        assert_eq!(
            evidence.locus_differences[2].nucleotide_support.support,
            [0.0; 4]
        );
        assert!(
            evidence.locus_differences[1]
                .nucleotide_support
                .mean_profile
                .is_none()
        );
        assert!(
            evidence.locus_differences[2]
                .nucleotide_support
                .mean_profile
                .is_none()
        );
        assert!(
            evidence.locus_differences[1]
                .nucleotide_support
                .heterogeneity
                .is_none()
        );
        assert!(
            evidence.locus_differences[2]
                .nucleotide_support
                .heterogeneity
                .is_none()
        );
        assert!(
            evidence.locus_differences[1]
                .nucleotide_support
                .directional_profile_distance
                .is_none()
        );
        assert!(
            evidence.locus_differences[2]
                .nucleotide_support
                .directional_profile_distance
                .is_none()
        );
        assert_eq!(
            evidence.locus_differences[0].observations[0].state,
            crate::model::sample_evidence::LocusState::Alternate
        );
        assert_eq!(
            evidence.locus_differences[1].observations[0].state,
            crate::model::sample_evidence::LocusState::Unresolved
        );
        assert!(
            evidence.locus_differences[1].observations[0]
                .signal
                .as_ref()
                .is_some_and(|signal| signal.profile.is_none())
        );
        assert_eq!(
            evidence.locus_differences[1].observations[0].nucleotide_contribution,
            crate::model::sample_evidence::NucleotideContribution::MissingProfile
        );
        assert_eq!(
            evidence.locus_differences[2].observations[0].state,
            crate::model::sample_evidence::LocusState::Deletion
        );
        assert!(
            evidence.locus_differences[2].observations[0]
                .signal
                .is_none()
        );
        assert_eq!(
            evidence.locus_differences[2].observations[0].nucleotide_contribution,
            crate::model::sample_evidence::NucleotideContribution::DeletionEvent
        );
        Ok(())
    }

    #[test]
    fn unresolved_call_with_profile_remains_nucleotide_eligible() -> Result<()> {
        let read = observation(
            "a",
            "reference",
            "config",
            Orientation::Forward,
            vec![column('N', 'T', Some(0), 12)],
            Vec::new(),
        );

        let evidence = aggregate(&[read], &sample_config())?;

        assert_eq!(evidence.locus_differences.len(), 1);
        assert_eq!(
            evidence.locus_differences[0].observations[0].state,
            crate::model::sample_evidence::LocusState::Unresolved
        );
        assert_eq!(
            evidence.locus_differences[0].support_topology.profile_reads,
            1
        );
        assert_eq!(
            evidence.locus_differences[0].observations[0].nucleotide_contribution,
            crate::model::sample_evidence::NucleotideContribution::Eligible
        );
        assert_eq!(
            evidence.locus_differences[0]
                .nucleotide_support
                .contributors,
            1
        );
        assert_eq!(
            evidence.locus_differences[0].nucleotide_support.support,
            [0.1, 0.2, 0.3, 0.4]
        );
        assert_eq!(
            evidence.locus_differences[0]
                .nucleotide_support
                .mean_profile
                .map(|profile| profile.weights),
            Some([0.1, 0.2, 0.3, 0.4])
        );
        let heterogeneity = evidence.locus_differences[0]
            .nucleotide_support
            .heterogeneity
            .ok_or_else(|| Error::Sample("single-read profile heterogeneity is missing".into()))?;
        assert!((heterogeneity.within_profile_impurity - 0.7).abs() < 1e-12);
        assert_eq!(heterogeneity.between_profile_dispersion, 0.0);
        assert!((heterogeneity.total_profile_heterogeneity - 0.7).abs() < 1e-12);
        assert!(
            evidence.locus_differences[0]
                .nucleotide_support
                .directional_profile_distance
                .is_none()
        );
        Ok(())
    }

    #[test]
    fn differential_locus_retains_reference_support_from_overlapping_reads() -> Result<()> {
        let alternate = observation(
            "a",
            "reference",
            "config",
            Orientation::Forward,
            vec![column('G', 'A', Some(0), 72)],
            vec![snv(73, "A", "G")],
        );
        let reference = observation(
            "b",
            "reference",
            "config",
            Orientation::Reverse,
            vec![column('A', 'A', Some(0), 72)],
            Vec::new(),
        );

        let evidence = aggregate(&[reference, alternate], &sample_config())?;

        assert_eq!(evidence.locus_differences.len(), 1);
        let topology = evidence.locus_differences[0].support_topology;
        assert_eq!(topology.reads, 2);
        assert_eq!(topology.forward_reads, 1);
        assert_eq!(topology.reverse_reads, 1);
        assert_eq!(topology.reference_reads, 1);
        assert_eq!(topology.alternate_reads, 1);
        assert_eq!(topology.unresolved_reads, 0);
        assert_eq!(topology.deletion_reads, 0);
        assert_eq!(topology.profile_reads, 2);
        assert_eq!(topology.profile_forward_reads, 1);
        assert_eq!(topology.profile_reverse_reads, 1);
        let observations = &evidence.locus_differences[0].observations;
        assert_eq!(observations.len(), 2);
        assert_eq!(observations[0].read_index, 0);
        assert_eq!(
            observations[0].state,
            crate::model::sample_evidence::LocusState::Alternate
        );
        assert_eq!(observations[1].read_index, 1);
        assert_eq!(
            observations[1].state,
            crate::model::sample_evidence::LocusState::Reference
        );
        assert_eq!(observations[1].quality, Some(50));
        Ok(())
    }

    #[test]
    fn all_reference_overlap_needs_no_per_locus_records() -> Result<()> {
        let first = observation(
            "a",
            "reference",
            "config",
            Orientation::Forward,
            vec![column('A', 'A', Some(0), 10), column('C', 'C', Some(1), 11)],
            Vec::new(),
        );
        let second = observation(
            "b",
            "reference",
            "config",
            Orientation::Reverse,
            vec![column('A', 'A', Some(0), 10), column('C', 'C', Some(1), 11)],
            Vec::new(),
        );

        let evidence = aggregate(&[first, second], &sample_config())?;

        assert!(evidence.locus_differences.is_empty());
        Ok(())
    }

    #[test]
    fn preserves_filtered_variant_observation_without_reporting_it() -> Result<()> {
        let variant = snv(73, "A", "G");
        let forward = observation(
            "a",
            "reference",
            "config",
            Orientation::Forward,
            vec![column('G', 'A', Some(0), 72)],
            vec![variant.clone()],
        );
        let mut reverse = observation(
            "b",
            "reference",
            "config",
            Orientation::Reverse,
            vec![column('G', 'A', Some(0), 72)],
            vec![variant],
        );
        reverse.variants.reported.clear();
        reverse.variants.observed[0].exclusion_reasons =
            vec![VariantExclusionReason::PeakBelowMinimum];

        let evidence = aggregate(&[forward, reverse], &sample_config())?;

        assert_eq!(evidence.variants[0].support.len(), 2);
        assert_eq!(evidence.variants[0].support_topology.reads, 2);
        assert_eq!(evidence.variants[0].support_topology.eligible_reads, 1);
        assert_eq!(evidence.variants[0].support_topology.forward_reads, 1);
        assert_eq!(evidence.variants[0].support_topology.reverse_reads, 1);
        assert_eq!(
            evidence.variants[0].support_topology.eligible_forward_reads,
            1
        );
        assert_eq!(
            evidence.variants[0].support_topology.eligible_reverse_reads,
            0
        );
        assert!(evidence.variants[0].support[0].eligible);
        assert!(!evidence.variants[0].support[1].eligible);
        assert_eq!(
            evidence.variants[0].support[1].exclusion_reasons,
            vec![VariantExclusionReason::PeakBelowMinimum]
        );
        assert_eq!(evidence.variants[0].support[0].calls[0].base, 'G');
        assert_eq!(
            evidence.variants[0].support[0].calls[0].peak_heights,
            [10, 20, 100, 30]
        );
        assert_eq!(evidence.variants[0].support[0].calls[0].quality, 50);
        Ok(())
    }

    #[test]
    fn rejects_misindexed_locus_evidence() {
        let mut read = observation(
            "a",
            "reference",
            "config",
            Orientation::Forward,
            vec![column('G', 'A', Some(0), 72)],
            Vec::new(),
        );
        read.signal.loci[0].call_index_0based = 1;

        assert!(aggregate(&[read], &sample_config()).is_err());
    }

    #[test]
    fn rejects_duplicate_reference_coordinate_within_one_read() {
        let read = observation(
            "a",
            "reference",
            "config",
            Orientation::Forward,
            vec![column('G', 'A', Some(0), 72), column('A', 'A', Some(1), 72)],
            Vec::new(),
        );

        assert!(aggregate(&[read], &sample_config()).is_err());
    }

    #[test]
    fn rejects_duplicate_normalized_variant_identity_within_one_read() {
        let variant = snv(73, "A", "G");
        let read = observation(
            "a",
            "reference",
            "config",
            Orientation::Forward,
            vec![column('G', 'A', Some(0), 72)],
            vec![variant.clone(), variant],
        );

        assert!(aggregate(&[read], &sample_config()).is_err());
    }

    #[test]
    fn rejects_incompatible_or_duplicate_reads_even_when_renamed() {
        let first = observation(
            "a",
            "reference",
            "config",
            Orientation::Forward,
            vec![column('A', 'A', Some(0), 0)],
            Vec::new(),
        );
        let incompatible = observation(
            "b",
            "other-reference",
            "config",
            Orientation::Forward,
            vec![column('A', 'A', Some(0), 0)],
            Vec::new(),
        );
        assert!(aggregate(&[first.clone(), incompatible], &sample_config()).is_err());

        let mut duplicate = observation(
            "a",
            "reference",
            "config",
            Orientation::Reverse,
            vec![column('A', 'A', Some(0), 0)],
            Vec::new(),
        );
        duplicate.input_name = "renamed-copy.ab1".into();
        assert!(aggregate(&[first, duplicate], &sample_config()).is_err());
    }
}
