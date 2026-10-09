//! Operational metrics derived from completed sample evidence.

use std::collections::BTreeMap;
use std::fmt;

use crate::model::attachment::SangerAttachment;
use crate::model::locus_evidence::LocusEvidence;
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

/// Sanger evidence joined to the core sample evidence by read identity and
/// call index.
struct Joined<'a> {
    evidence: &'a SampleEvidence,
    sanger: &'a BTreeMap<String, SangerAttachment>,
}

impl<'a> Joined<'a> {
    fn attachment(&self, read_index: usize) -> Option<&'a SangerAttachment> {
        let read = self.evidence.reads.get(read_index)?;
        self.sanger.get(&read.input_sha256)
    }

    /// Sanger attachment and call index of every locus observation of a call.
    fn locus_calls(&self) -> impl Iterator<Item = (&'a SangerAttachment, usize)> + '_ {
        self.evidence
            .locus_differences
            .iter()
            .flat_map(|difference| &difference.observations)
            .filter_map(|observation| {
                let call = observation.call_index_0based?;
                Some((self.attachment(observation.read_index)?, call))
            })
    }

    /// Sanger attachment and call index of every variant call.
    fn variant_calls(&self) -> impl Iterator<Item = (&'a SangerAttachment, usize)> + '_ {
        self.evidence
            .variants
            .iter()
            .flat_map(|variant| &variant.support)
            .flat_map(|support| {
                support.calls.iter().filter_map(|call| {
                    Some((self.attachment(support.read_index)?, call.call_index_0based))
                })
            })
    }
}

fn loci<'a>(
    calls: impl Iterator<Item = (&'a SangerAttachment, usize)>,
) -> impl Iterator<Item = &'a LocusEvidence> {
    calls.filter_map(|(attachment, call)| attachment.locus(call))
}

fn positive(values: impl Iterator<Item = f64>) -> usize {
    values.filter(|value| *value > 0.0).count()
}

pub(super) fn summarize(
    evidence: &SampleEvidence,
    sanger: &BTreeMap<String, SangerAttachment>,
) -> SampleAggregationMetrics {
    let joined = Joined { evidence, sanger };
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
        profiled_variant_calls: loci(joined.variant_calls())
            .filter(|locus| locus.profile.is_some())
            .count(),
        noisy_locus_observations: joined
            .locus_calls()
            .filter(|(attachment, call)| attachment.in_noisy_region(*call))
            .count(),
        noisy_variant_calls: joined
            .variant_calls()
            .filter(|(attachment, call)| attachment.in_noisy_region(*call))
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
        locus_positive_corrected_channels: positive(
            loci(joined.locus_calls()).flat_map(|locus| locus.corrected_amplitudes),
        ),
        locus_positive_snr_channels: positive(
            loci(joined.locus_calls()).flat_map(|locus| locus.snrs),
        ),
        variant_positive_corrected_channels: positive(
            loci(joined.variant_calls()).flat_map(|locus| locus.corrected_amplitudes),
        ),
        variant_positive_snr_channels: positive(
            loci(joined.variant_calls()).flat_map(|locus| locus.snrs),
        ),
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
        masked_calls_total: (0..evidence.reads.len())
            .filter_map(|index| joined.attachment(index))
            .map(|attachment| attachment.callability.masked_count())
            .sum(),
        callable_calls_total: (0..evidence.reads.len())
            .filter_map(|index| joined.attachment(index))
            .map(|attachment| attachment.callability.callable_count())
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
