//! Unweighted accumulation of eligible nucleotide profiles at retained loci.

use crate::model::alignment::Orientation;
use crate::model::called_read::CalledRead;
use crate::model::sample_evidence::{
    LocusNucleotideSupport, NucleotideContribution, ProfileHeterogeneity, SampleLocusObservation,
};
use dna_kernel::error::{Result, SampleError};
use dna_kernel::read_evidence::EvidenceProfile;

use super::profile_geometry;

/// Adds each eligible normalized profile with unit read mass.
pub(super) fn aggregate(
    observations: &[SampleLocusObservation],
    reads: &[&CalledRead],
) -> Result<LocusNucleotideSupport> {
    let mut result = LocusNucleotideSupport {
        contributors: 0,
        forward_contributors: 0,
        reverse_contributors: 0,
        support: [0.0; 4],
        forward_support: [0.0; 4],
        reverse_support: [0.0; 4],
        mean_profile: None,
        forward_mean_profile: None,
        reverse_mean_profile: None,
        heterogeneity: None,
        forward_heterogeneity: None,
        reverse_heterogeneity: None,
        directional_profile_distance: None,
    };
    let mut impurity_sum = 0.0;
    let mut forward_impurity_sum = 0.0;
    let mut reverse_impurity_sum = 0.0;

    for observation in observations {
        if observation.nucleotide_contribution != NucleotideContribution::Eligible {
            continue;
        }
        let profile = observation.profile.ok_or(SampleError::MissingProfile {
            read: observation.read_index,
        })?;
        let read = reads
            .get(observation.read_index)
            .ok_or(SampleError::MissingRead {
                context: "nucleotide support",
                index: observation.read_index,
            })?;

        let impurity = profile_geometry::impurity(profile);
        impurity_sum += impurity;
        result.contributors += 1;
        match read.alignment.orientation {
            Orientation::Forward => {
                add(&mut result.forward_support, profile.weights);
                forward_impurity_sum += impurity;
                result.forward_contributors += 1;
            }
            Orientation::Reverse => {
                add(&mut result.reverse_support, profile.weights);
                reverse_impurity_sum += impurity;
                result.reverse_contributors += 1;
            }
        }
    }

    result.support = std::array::from_fn(|channel| {
        result.forward_support[channel] + result.reverse_support[channel]
    });
    result.mean_profile = mean_profile(result.support, result.contributors);
    result.forward_mean_profile = mean_profile(result.forward_support, result.forward_contributors);
    result.reverse_mean_profile = mean_profile(result.reverse_support, result.reverse_contributors);
    result.heterogeneity =
        profile_geometry::heterogeneity(result.mean_profile, impurity_sum, result.contributors);
    result.forward_heterogeneity = profile_geometry::heterogeneity(
        result.forward_mean_profile,
        forward_impurity_sum,
        result.forward_contributors,
    );
    result.reverse_heterogeneity = profile_geometry::heterogeneity(
        result.reverse_mean_profile,
        reverse_impurity_sum,
        result.reverse_contributors,
    );
    result.directional_profile_distance =
        profile_geometry::total_variation(result.forward_mean_profile, result.reverse_mean_profile);

    if result.contributors != result.forward_contributors + result.reverse_contributors
        || result.mean_profile.is_some() != (result.contributors > 0)
        || result.forward_mean_profile.is_some() != (result.forward_contributors > 0)
        || result.reverse_mean_profile.is_some() != (result.reverse_contributors > 0)
        || result.heterogeneity.is_some() != result.mean_profile.is_some()
        || result.forward_heterogeneity.is_some() != result.forward_mean_profile.is_some()
        || result.reverse_heterogeneity.is_some() != result.reverse_mean_profile.is_some()
        || result.directional_profile_distance.is_some()
            != (result.forward_mean_profile.is_some() && result.reverse_mean_profile.is_some())
        || !geometry_is_valid(result.heterogeneity)
        || !geometry_is_valid(result.forward_heterogeneity)
        || !geometry_is_valid(result.reverse_heterogeneity)
        || !result
            .directional_profile_distance
            .is_none_or(|distance| distance.is_finite() && (0.0..=1.0).contains(&distance))
        || !result
            .support
            .iter()
            .chain(result.forward_support.iter())
            .chain(result.reverse_support.iter())
            .all(|value| value.is_finite() && *value >= 0.0)
    {
        return Err(
            SampleError::Inconsistent("nucleotide support accumulation is inconsistent").into(),
        );
    }

    Ok(result)
}

fn geometry_is_valid(geometry: Option<ProfileHeterogeneity>) -> bool {
    const EPSILON: f64 = 1e-12;

    geometry.is_none_or(|geometry| {
        [
            geometry.within_profile_impurity,
            geometry.between_profile_dispersion,
            geometry.total,
        ]
        .into_iter()
        .all(|value| value.is_finite() && (0.0..=0.75 + EPSILON).contains(&value))
            && (geometry.within_profile_impurity + geometry.between_profile_dispersion
                - geometry.total)
                .abs()
                <= EPSILON
    })
}

fn mean_profile(support: [f64; 4], contributors: usize) -> Option<EvidenceProfile> {
    (contributors > 0).then(|| EvidenceProfile {
        weights: support.map(|value| value / contributors as f64),
    })
}

fn add(target: &mut [f64; 4], contribution: [f64; 4]) {
    for (target, contribution) in target.iter_mut().zip(contribution) {
        *target += contribution;
    }
}
