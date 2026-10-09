//! Operational metrics derived from completed sample evidence.

use std::fmt;

use crate::model::sample_evidence::{NucleotideContribution, SampleEvidence};

/// Aggregation metrics rendered as the `key=value` tail of the
/// `sample_aggregation_completed` record.
pub(super) struct SampleAggregationMetrics {
    profiled_locus_observations: usize,
    profiled_locus_forward_reads: usize,
    profiled_locus_reverse_reads: usize,
    profiled_variant_calls: usize,
    noisy_locus_observations: usize,
    noisy_variant_calls: usize,
    eligible_nucleotide_locus_observations: usize,
    missing_profile_locus_observations: usize,
    deletion_event_locus_observations: usize,
    nucleotide_support_loci: usize,
    bidirectional_nucleotide_support_loci: usize,
    unweighted_nucleotide_profile_mass: f64,
    profile_geometry_loci: usize,
    within_profile_impurity_sum: f64,
    between_profile_dispersion_sum: f64,
    total_profile_heterogeneity_sum: f64,
    forward_profile_geometry_loci: usize,
    reverse_profile_geometry_loci: usize,
    directional_profile_distance_loci: usize,
    directional_profile_distance_sum: f64,
    locus_positive_corrected_channels: usize,
    locus_positive_snr_channels: usize,
    variant_positive_corrected_channels: usize,
    variant_positive_snr_channels: usize,
    locus_forward_reads: usize,
    locus_reverse_reads: usize,
    locus_reference_reads: usize,
    locus_alternate_reads: usize,
    locus_unresolved_reads: usize,
    locus_deletion_reads: usize,
    variants: usize,
    masked_calls_total: usize,
    callable_calls_total: usize,
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
                acc.3 + geometry.total,
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
        variants: evidence.variants.len(),
        masked_calls_total: evidence
            .reads
            .iter()
            .map(|read| read.callability.masked_count())
            .sum(),
        callable_calls_total: evidence
            .reads
            .iter()
            .map(|read| read.callability.callable_count())
            .sum(),
    }
}

impl fmt::Display for SampleAggregationMetrics {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            concat!(
                "profiled_locus_observations={} profiled_locus_forward_reads={} ",
                "profiled_locus_reverse_reads={} noisy_locus_observations={} ",
                "eligible_nucleotide_locus_observations={} missing_profile_locus_observations={} ",
                "deletion_event_locus_observations={} nucleotide_support_loci={} ",
                "bidirectional_nucleotide_support_loci={} unweighted_nucleotide_profile_mass={:.6} ",
                "profile_geometry_loci={} within_profile_impurity_sum={:.6} ",
                "between_profile_dispersion_sum={:.6} total_profile_heterogeneity_sum={:.6} ",
                "forward_profile_geometry_loci={} reverse_profile_geometry_loci={} ",
                "directional_profile_distance_loci={} directional_profile_distance_sum={:.6} ",
                "locus_positive_corrected_channels={} locus_positive_snr_channels={} ",
                "locus_forward_reads={} locus_reverse_reads={} locus_reference_reads={} ",
                "locus_alternate_reads={} locus_unresolved_reads={} locus_deletion_reads={} variants={} ",
                "profiled_variant_calls={} noisy_variant_calls={} variant_positive_corrected_channels={} ",
                "variant_positive_snr_channels={} masked_calls_total={} callable_calls_total={}"
            ),
            self.profiled_locus_observations,
            self.profiled_locus_forward_reads,
            self.profiled_locus_reverse_reads,
            self.noisy_locus_observations,
            self.eligible_nucleotide_locus_observations,
            self.missing_profile_locus_observations,
            self.deletion_event_locus_observations,
            self.nucleotide_support_loci,
            self.bidirectional_nucleotide_support_loci,
            self.unweighted_nucleotide_profile_mass,
            self.profile_geometry_loci,
            self.within_profile_impurity_sum,
            self.between_profile_dispersion_sum,
            self.total_profile_heterogeneity_sum,
            self.forward_profile_geometry_loci,
            self.reverse_profile_geometry_loci,
            self.directional_profile_distance_loci,
            self.directional_profile_distance_sum,
            self.locus_positive_corrected_channels,
            self.locus_positive_snr_channels,
            self.locus_forward_reads,
            self.locus_reverse_reads,
            self.locus_reference_reads,
            self.locus_alternate_reads,
            self.locus_unresolved_reads,
            self.locus_deletion_reads,
            self.variants,
            self.profiled_variant_calls,
            self.noisy_variant_calls,
            self.variant_positive_corrected_channels,
            self.variant_positive_snr_channels,
            self.masked_calls_total,
            self.callable_calls_total
        )
    }
}
