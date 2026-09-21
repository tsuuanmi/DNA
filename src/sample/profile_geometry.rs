//! Threshold-free geometry over normalized nucleotide evidence profiles.

use crate::model::locus_evidence::EvidenceProfile;
use crate::model::sample_evidence::ProfileHeterogeneity;

/// Gini impurity of one normalized A/C/G/T profile.
pub(super) fn impurity(profile: EvidenceProfile) -> f64 {
    1.0 - profile
        .weights
        .into_iter()
        .map(|weight| weight * weight)
        .sum::<f64>()
}

/// Decomposes total profile heterogeneity into within-read mixture and between-read dispersion.
pub(super) fn heterogeneity(
    mean_profile: Option<EvidenceProfile>,
    impurity_sum: f64,
    contributors: usize,
) -> Option<ProfileHeterogeneity> {
    const FLOAT_TOLERANCE: f64 = 1e-12;

    let mean_profile = mean_profile.filter(|_| contributors > 0)?;
    let within_profile_impurity = impurity_sum / contributors as f64;
    let total_profile_heterogeneity = impurity(mean_profile);
    let raw_between = total_profile_heterogeneity - within_profile_impurity;
    let between_profile_dispersion = if raw_between < 0.0 && raw_between.abs() <= FLOAT_TOLERANCE {
        0.0
    } else {
        raw_between
    };

    Some(ProfileHeterogeneity {
        within_profile_impurity,
        between_profile_dispersion,
        total_profile_heterogeneity,
    })
}

/// Total Variation distance between two normalized A/C/G/T profile distributions.
pub(super) fn total_variation(
    left: Option<EvidenceProfile>,
    right: Option<EvidenceProfile>,
) -> Option<f64> {
    Some(
        0.5 * left?
            .weights
            .into_iter()
            .zip(right?.weights)
            .map(|(left, right)| (left - right).abs())
            .sum::<f64>(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn profile(weights: [f64; 4]) -> EvidenceProfile {
        EvidenceProfile { weights }
    }

    #[test]
    fn separates_replicated_mixture_from_between_read_disagreement() -> Result<(), &'static str> {
        let mixed = profile([0.5, 0.0, 0.5, 0.0]);
        let replicated = heterogeneity(Some(mixed), impurity(mixed) * 2.0, 2)
            .ok_or("replicated mixture geometry is missing")?;
        assert_eq!(replicated.within_profile_impurity, 0.5);
        assert_eq!(replicated.between_profile_dispersion, 0.0);
        assert_eq!(replicated.total_profile_heterogeneity, 0.5);

        let pure_a = profile([1.0, 0.0, 0.0, 0.0]);
        let pure_g = profile([0.0, 0.0, 1.0, 0.0]);
        let disagreement = heterogeneity(Some(mixed), impurity(pure_a) + impurity(pure_g), 2)
            .ok_or("inter-read disagreement geometry is missing")?;
        assert_eq!(disagreement.within_profile_impurity, 0.0);
        assert_eq!(disagreement.between_profile_dispersion, 0.5);
        assert_eq!(disagreement.total_profile_heterogeneity, 0.5);
        Ok(())
    }

    #[test]
    fn total_variation_is_zero_for_identical_profiles_and_one_for_disjoint_profiles() {
        let pure_a = profile([1.0, 0.0, 0.0, 0.0]);
        let pure_g = profile([0.0, 0.0, 1.0, 0.0]);

        assert_eq!(total_variation(Some(pure_a), Some(pure_a)), Some(0.0));
        assert_eq!(total_variation(Some(pure_a), Some(pure_g)), Some(1.0));
        assert_eq!(total_variation(Some(pure_a), None), None);
    }
}
