//! Operational metrics derived from completed sample evidence.

use crate::model::sample_evidence::{NucleotideContribution, SampleEvidence};

pub(super) struct SampleAggregationMetrics {
    pub(super) profiled_locus_observations: usize,
    pub(super) profiled_locus_forward_reads: usize,
    pub(super) profiled_locus_reverse_reads: usize,
    pub(super) profiled_variant_calls: usize,
    pub(super) noisy_locus_observations: usize,
    pub(super) noisy_variant_calls: usize,
    pub(super) eligible_nucleotide_locus_observations: usize,
    pub(super) missing_profile_locus_observations: usize,
    pub(super) deletion_event_locus_observations: usize,
    pub(super) nucleotide_support_loci: usize,
    pub(super) bidirectional_nucleotide_support_loci: usize,
    pub(super) unweighted_nucleotide_profile_mass: f64,
    pub(super) profile_geometry_loci: usize,
    pub(super) within_profile_impurity_sum: f64,
    pub(super) between_profile_dispersion_sum: f64,
    pub(super) total_profile_heterogeneity_sum: f64,
    pub(super) forward_profile_geometry_loci: usize,
    pub(super) reverse_profile_geometry_loci: usize,
    pub(super) directional_profile_distance_loci: usize,
    pub(super) directional_profile_distance_sum: f64,
    pub(super) locus_positive_corrected_channels: usize,
    pub(super) locus_positive_snr_channels: usize,
    pub(super) variant_positive_corrected_channels: usize,
    pub(super) variant_positive_snr_channels: usize,
    pub(super) locus_forward_reads: usize,
    pub(super) locus_reverse_reads: usize,
    pub(super) locus_reference_reads: usize,
    pub(super) locus_alternate_reads: usize,
    pub(super) locus_unresolved_reads: usize,
    pub(super) locus_deletion_reads: usize,
}

pub(super) fn summarize(evidence: &SampleEvidence) -> SampleAggregationMetrics {
    let (
        profile_geometry_loci,
        within_profile_impurity_sum,
        between_profile_dispersion_sum,
        total_profile_heterogeneity_sum,
    ) = evidence
        .locus_differences
        .iter()
        .filter_map(|difference| difference.nucleotide_support.heterogeneity)
        .fold((0usize, 0.0, 0.0, 0.0), |acc, geometry| {
            (
                acc.0 + 1,
                acc.1 + geometry.within_profile_impurity,
                acc.2 + geometry.between_profile_dispersion,
                acc.3 + geometry.total_profile_heterogeneity,
            )
        });
    let forward_profile_geometry_loci = evidence
        .locus_differences
        .iter()
        .filter(|difference| {
            difference
                .nucleotide_support
                .forward_heterogeneity
                .is_some()
        })
        .count();
    let reverse_profile_geometry_loci = evidence
        .locus_differences
        .iter()
        .filter(|difference| {
            difference
                .nucleotide_support
                .reverse_heterogeneity
                .is_some()
        })
        .count();
    let (directional_profile_distance_loci, directional_profile_distance_sum) = evidence
        .locus_differences
        .iter()
        .filter_map(|difference| difference.nucleotide_support.directional_profile_distance)
        .fold((0usize, 0.0), |(count, sum), distance| {
            (count + 1, sum + distance)
        });

    SampleAggregationMetrics {
        profiled_locus_observations: evidence
            .locus_differences
            .iter()
            .map(|difference| difference.support_topology.profile_reads)
            .sum(),
        profiled_locus_forward_reads: evidence
            .locus_differences
            .iter()
            .map(|difference| difference.support_topology.profile_forward_reads)
            .sum(),
        profiled_locus_reverse_reads: evidence
            .locus_differences
            .iter()
            .map(|difference| difference.support_topology.profile_reverse_reads)
            .sum(),
        profiled_variant_calls: evidence
            .variants
            .iter()
            .flat_map(|variant| &variant.support)
            .flat_map(|support| &support.calls)
            .filter(|call| call.signal.profile.is_some())
            .count(),
        noisy_locus_observations: evidence
            .locus_differences
            .iter()
            .flat_map(|difference| &difference.observations)
            .filter(|observation| {
                observation
                    .signal
                    .as_ref()
                    .is_some_and(|signal| signal.in_noisy_region)
            })
            .count(),
        noisy_variant_calls: evidence
            .variants
            .iter()
            .flat_map(|variant| &variant.support)
            .flat_map(|support| &support.calls)
            .filter(|call| call.signal.in_noisy_region)
            .count(),
        eligible_nucleotide_locus_observations: evidence
            .locus_differences
            .iter()
            .flat_map(|difference| &difference.observations)
            .filter(|observation| {
                observation.nucleotide_contribution == NucleotideContribution::Eligible
            })
            .count(),
        missing_profile_locus_observations: evidence
            .locus_differences
            .iter()
            .flat_map(|difference| &difference.observations)
            .filter(|observation| {
                observation.nucleotide_contribution == NucleotideContribution::MissingProfile
            })
            .count(),
        deletion_event_locus_observations: evidence
            .locus_differences
            .iter()
            .flat_map(|difference| &difference.observations)
            .filter(|observation| {
                observation.nucleotide_contribution == NucleotideContribution::DeletionEvent
            })
            .count(),
        nucleotide_support_loci: evidence
            .locus_differences
            .iter()
            .filter(|difference| difference.nucleotide_support.mean_profile.is_some())
            .count(),
        bidirectional_nucleotide_support_loci: evidence
            .locus_differences
            .iter()
            .filter(|difference| {
                difference.nucleotide_support.forward_mean_profile.is_some()
                    && difference.nucleotide_support.reverse_mean_profile.is_some()
            })
            .count(),
        unweighted_nucleotide_profile_mass: evidence
            .locus_differences
            .iter()
            .flat_map(|difference| difference.nucleotide_support.support)
            .sum(),
        profile_geometry_loci,
        within_profile_impurity_sum,
        between_profile_dispersion_sum,
        total_profile_heterogeneity_sum,
        forward_profile_geometry_loci,
        reverse_profile_geometry_loci,
        directional_profile_distance_loci,
        directional_profile_distance_sum,
        locus_positive_corrected_channels: evidence
            .locus_differences
            .iter()
            .flat_map(|difference| &difference.observations)
            .filter_map(|observation| observation.signal.as_ref())
            .flat_map(|signal| signal.corrected_amplitudes)
            .filter(|value| *value > 0.0)
            .count(),
        locus_positive_snr_channels: evidence
            .locus_differences
            .iter()
            .flat_map(|difference| &difference.observations)
            .filter_map(|observation| observation.signal.as_ref())
            .flat_map(|signal| signal.snrs)
            .filter(|value| *value > 0.0)
            .count(),
        variant_positive_corrected_channels: evidence
            .variants
            .iter()
            .flat_map(|variant| &variant.support)
            .flat_map(|support| &support.calls)
            .flat_map(|call| call.signal.corrected_amplitudes)
            .filter(|value| *value > 0.0)
            .count(),
        variant_positive_snr_channels: evidence
            .variants
            .iter()
            .flat_map(|variant| &variant.support)
            .flat_map(|support| &support.calls)
            .flat_map(|call| call.signal.snrs)
            .filter(|value| *value > 0.0)
            .count(),
        locus_forward_reads: evidence
            .locus_differences
            .iter()
            .map(|difference| difference.support_topology.forward_reads)
            .sum(),
        locus_reverse_reads: evidence
            .locus_differences
            .iter()
            .map(|difference| difference.support_topology.reverse_reads)
            .sum(),
        locus_reference_reads: evidence
            .locus_differences
            .iter()
            .map(|difference| difference.support_topology.reference_reads)
            .sum(),
        locus_alternate_reads: evidence
            .locus_differences
            .iter()
            .map(|difference| difference.support_topology.alternate_reads)
            .sum(),
        locus_unresolved_reads: evidence
            .locus_differences
            .iter()
            .map(|difference| difference.support_topology.unresolved_reads)
            .sum(),
        locus_deletion_reads: evidence
            .locus_differences
            .iter()
            .map(|difference| difference.support_topology.deletion_reads)
            .sum(),
    }
}
