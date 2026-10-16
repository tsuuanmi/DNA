//! Stratified estimates of a rate from a weighted panel.
//!
//! Each sample carries the inverse of its stratum's sampling share, so a
//! stratum's weights sum to its estimated size `N_h`. The rate is
//! `p = Σ w·y / Σ w`, which equals `Σ_h W_h·p_h` with `W_h = N_h / N`. Its
//! variance is the stratified variance with the finite-population correction,
//! `Σ_h W_h² · p_h(1 − p_h) / n_h · (1 − n_h / N_h)`, and the 95 % interval is
//! `p ± 1.96·√variance`. Samples of weight zero are left out.

use std::collections::BTreeMap;

/// Normal quantile of a two-sided 95 % interval.
const Z_95: f64 = 1.96;

/// One sample's contribution to an estimate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Observation<'a> {
    pub(crate) stratum: &'a str,
    pub(crate) weight: f64,
    /// Whether the counted event (for example "not exact") happened.
    pub(crate) event: bool,
}

/// An estimated rate over every stratum.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Estimate {
    /// Samples of positive weight.
    pub(crate) samples: usize,
    /// Estimated population, `Σ w`.
    pub(crate) population: f64,
    /// Samples with the event.
    pub(crate) events: usize,
    /// Weighted rate.
    pub(crate) rate: f64,
    /// Half-width of the 95 % interval.
    pub(crate) ci95: f64,
    /// The stratified variance of the rate.
    pub(crate) variance: f64,
    /// Unweighted rate over the same samples.
    pub(crate) unweighted: f64,
    /// Per-stratum rates, by stratum name.
    pub(crate) strata: Vec<Stratum>,
}

/// The estimate within one stratum.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Stratum {
    pub(crate) name: String,
    pub(crate) samples: usize,
    pub(crate) population: f64,
    pub(crate) events: usize,
    pub(crate) rate: f64,
}

/// Estimates the event rate, or `None` when no sample has a positive weight.
pub(crate) fn estimate(observations: &[Observation<'_>]) -> Option<Estimate> {
    #[derive(Default)]
    struct Sums {
        samples: usize,
        events: usize,
        weight: f64,
        weighted_events: f64,
    }
    let mut strata = BTreeMap::<&str, Sums>::new();
    for observation in observations
        .iter()
        .filter(|observation| observation.weight > 0.0)
    {
        let sums = strata.entry(observation.stratum).or_default();
        sums.samples += 1;
        sums.weight += observation.weight;
        if observation.event {
            sums.events += 1;
            sums.weighted_events += observation.weight;
        }
    }
    let population = strata.values().map(|sums| sums.weight).sum::<f64>();
    if population <= 0.0 {
        return None;
    }
    let mut rate = 0.0;
    let mut variance = 0.0;
    let strata = strata
        .into_iter()
        .map(|(name, sums)| {
            let share = sums.weight / population;
            let stratum_rate = sums.weighted_events / sums.weight;
            let samples = sums.samples as f64;
            rate += share * stratum_rate;
            variance += share * share * stratum_rate * (1.0 - stratum_rate) / samples
                * (1.0 - samples / sums.weight).max(0.0);
            Stratum {
                name: name.to_owned(),
                samples: sums.samples,
                population: sums.weight,
                events: sums.events,
                rate: stratum_rate,
            }
        })
        .collect::<Vec<_>>();
    let samples = strata.iter().map(|stratum| stratum.samples).sum::<usize>();
    let events = strata.iter().map(|stratum| stratum.events).sum::<usize>();
    Some(Estimate {
        samples,
        population,
        events,
        rate,
        ci95: Z_95 * variance.sqrt(),
        variance,
        unweighted: events as f64 / samples as f64,
        strata,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn observations(
        stratum: &str,
        weight: f64,
        samples: usize,
        events: usize,
    ) -> Vec<Observation<'_>> {
        (0..samples)
            .map(|index| Observation {
                stratum,
                weight,
                event: index < events,
            })
            .collect()
    }

    fn close(left: f64, right: f64) -> bool {
        (left - right).abs() < 1e-12
    }

    #[test]
    fn weights_strata_by_their_population() -> Result<(), &'static str> {
        // A: N = 100, n = 10, 2 events; B: N = 50, n = 10, 5 events; the
        // zero-weight row is left out.
        let mut all = observations("A", 10.0, 10, 2);
        all.extend(observations("B", 5.0, 10, 5));
        all.extend(observations("C", 0.0, 3, 3));
        let estimate = estimate(&all).ok_or("no estimate")?;
        assert_eq!((estimate.samples, estimate.events), (20, 7));
        assert!(close(estimate.population, 150.0));
        assert!(close(estimate.rate, (100.0 * 0.2 + 50.0 * 0.5) / 150.0));
        let variance = (100.0_f64 / 150.0).powi(2) * 0.2 * 0.8 / 10.0 * 0.9
            + (50.0_f64 / 150.0).powi(2) * 0.5 * 0.5 / 10.0 * 0.8;
        assert!(close(estimate.ci95, 1.96 * variance.sqrt()));
        assert!(close(estimate.unweighted, 0.35));
        assert_eq!(estimate.strata.len(), 2);
        assert_eq!(estimate.strata[1].name, "B");
        assert!(close(estimate.strata[1].rate, 0.5));
        Ok(())
    }

    #[test]
    fn a_fully_sampled_stratum_has_no_variance() -> Result<(), &'static str> {
        let estimate = estimate(&observations("A", 1.0, 4, 1)).ok_or("no estimate")?;
        assert!(close(estimate.rate, 0.25));
        assert!(close(estimate.ci95, 0.0));
        Ok(())
    }

    #[test]
    fn nothing_to_estimate_without_positive_weight() {
        assert_eq!(estimate(&observations("A", 0.0, 4, 1)), None);
        assert_eq!(estimate(&[]), None);
    }
}
