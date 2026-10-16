//! Bounds for the approval endpoints (PROP-0005 approval protocol).
//!
//! - [`beta_upper`]: the one-sided Clopper–Pearson upper bound for a
//!   possibly non-integer event count, through the beta quantile.
//! - [`korn_graubard`]: the Korn–Graubard bound for a stratified proportion,
//!   which applies [`beta_upper`] at the effective sample size of the
//!   stratified variance, so a rate of zero still has a positive bound.
//! - [`bootstrap`]: the 95th percentile of the corpus-weighted leak rate
//!   over batch resamples, with each draw seeded by its index.

use crate::triage::score::ALPHA;

/// The one-sided `1 − ALPHA` upper bound on a proportion with `events` of
/// `trials`, `events` possibly fractional: the `1 − ALPHA` quantile of
/// `Beta(events + 1, trials − events)`; `1` when `events >= trials` or there
/// is no trial.
pub(crate) fn beta_upper(events: f64, trials: f64) -> f64 {
    if trials <= 0.0 || events >= trials {
        return 1.0;
    }
    let (a, b) = (events.max(0.0) + 1.0, trials - events.max(0.0));
    let (mut low, mut high) = (0.0_f64, 1.0_f64);
    for _ in 0..200 {
        let middle = f64::midpoint(low, high);
        if incomplete_beta(a, b, middle) < 1.0 - ALPHA {
            low = middle;
        } else {
            high = middle;
        }
    }
    high
}

/// The Korn–Graubard one-sided upper bound of a stratified proportion `rate`
/// with variance `variance` over `samples` samples: the effective sample size
/// is `rate (1 − rate) / variance`, or `samples` when the rate or the
/// variance is zero, and never more than `samples`.
pub(crate) fn korn_graubard(rate: f64, variance: f64, samples: usize) -> f64 {
    let samples = samples as f64;
    if samples <= 0.0 {
        return 1.0;
    }
    let effective = if rate <= 0.0 || variance <= 0.0 {
        samples
    } else {
        (rate * (1.0 - rate) / variance).min(samples)
    };
    beta_upper(rate * effective, effective)
}

/// `SplitMix64`: a small, fully specified generator, so an outside reviewer
/// can reproduce every draw from the seed.
#[derive(Debug, Clone)]
pub(crate) struct SplitMix64(u64);

impl SplitMix64 {
    /// The generator of draw `draw` under `seed`, independent of the order in
    /// which draws run.
    pub(crate) fn for_draw(seed: u64, draw: u64) -> Self {
        let mut mixer = Self(seed ^ draw.wrapping_mul(0x9E37_79B9_7F4A_7C15));
        let _ = mixer.next_u64();
        mixer
    }

    pub(crate) fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// A uniform index below `bound`, by rejection so it has no modulo bias.
    pub(crate) fn below(&mut self, bound: usize) -> usize {
        let bound = bound as u64;
        if bound == 0 {
            return 0;
        }
        let zone = u64::MAX - u64::MAX % bound;
        loop {
            let value = self.next_u64();
            if value < zone {
                return usize::try_from(value % bound).unwrap_or(0);
            }
        }
    }
}

/// One sample's contribution to the bootstrap: its batch, corpus weight,
/// and whether it passed and leaked.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Draw {
    pub(crate) batch: usize,
    pub(crate) weight: f64,
    pub(crate) passed: bool,
    pub(crate) leaked: bool,
}

/// The bootstrap of the corpus-weighted leak rate among passed samples.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Bootstrap {
    /// The one-sided 95th percentile of the draws' rates.
    pub(crate) upper: f64,
    /// Draws that held no passed sample; they count as rate 1.
    pub(crate) empty_draws: usize,
}

/// Resamples the panel's `batches` batches with replacement `draws` times,
/// seeding draw `d` with `(seed, d)`. Each draw's rate is the combined ratio
/// `Σ w·leak / Σ w·pass` with the samples' fixed weights; a draw without a
/// passed sample counts as rate 1.
pub(crate) fn bootstrap(samples: &[Draw], batches: usize, draws: u64, seed: u64) -> Bootstrap {
    let mut members = vec![Vec::new(); batches];
    for sample in samples {
        if let Some(batch) = members.get_mut(sample.batch) {
            batch.push(*sample);
        }
    }
    // Per batch, its weighted passes and weighted leaks.
    let sums = members
        .iter()
        .map(|batch| {
            batch.iter().fold((0.0, 0.0), |(pass, leak), sample| {
                (
                    pass + if sample.passed { sample.weight } else { 0.0 },
                    leak + if sample.passed && sample.leaked {
                        sample.weight
                    } else {
                        0.0
                    },
                )
            })
        })
        .collect::<Vec<_>>();
    let mut rates = Vec::new();
    let mut empty_draws = 0;
    for draw in 0..draws {
        let mut generator = SplitMix64::for_draw(seed, draw);
        let (mut pass, mut leak) = (0.0, 0.0);
        for _ in 0..batches {
            let (batch_pass, batch_leak) = sums[generator.below(batches)];
            pass += batch_pass;
            leak += batch_leak;
        }
        if pass > 0.0 {
            rates.push(leak / pass);
        } else {
            empty_draws += 1;
            rates.push(1.0);
        }
    }
    rates.sort_by(f64::total_cmp);
    // The `1 − ALPHA` quantile: index `⌈len·(1 − ALPHA)⌉ − 1`, which is
    // `len − ⌊len·ALPHA⌋ − 1`.
    let index = rates
        .len()
        .saturating_sub(largest_count(rates.len(), ALPHA) + 1);
    Bootstrap {
        upper: rates.get(index).copied().unwrap_or(1.0),
        empty_draws,
    }
}

/// The largest count of `size` items whose share does not exceed `rate`:
/// `⌊size·rate⌋`, computed without a float-to-integer cast.
pub(crate) fn largest_count(size: usize, rate: f64) -> usize {
    (0..=size)
        .take_while(|count| *count as f64 <= rate * size as f64)
        .last()
        .unwrap_or(0)
}

/// The regularized incomplete beta function `I_x(a, b)`, by Lentz's
/// continued fraction (Numerical Recipes, `betacf`).
fn incomplete_beta(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let front =
        (ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a * x.ln() + b * (1.0 - x).ln()).exp();
    if x < (a + 1.0) / (a + b + 2.0) {
        front * continued_fraction(a, b, x) / a
    } else {
        1.0 - front * continued_fraction(b, a, 1.0 - x) / b
    }
}

fn continued_fraction(shape: f64, other: f64, point: f64) -> f64 {
    const TINY: f64 = 1e-300;
    let guard = |value: f64| if value.abs() < TINY { TINY } else { value };
    let (sum, above, below) = (shape + other, shape + 1.0, shape - 1.0);
    let mut numerator = 1.0;
    let mut denominator = 1.0 / guard(1.0 - sum * point / above);
    let mut fraction = denominator;
    for step in 1..=300 {
        let step = f64::from(step);
        let twice = 2.0 * step;
        let even = step * (other - step) * point / ((below + twice) * (shape + twice));
        denominator = 1.0 / guard(1.0 + even * denominator);
        numerator = guard(1.0 + even / numerator);
        fraction *= denominator * numerator;
        let odd = -(shape + step) * (sum + step) * point / ((shape + twice) * (above + twice));
        denominator = 1.0 / guard(1.0 + odd * denominator);
        numerator = guard(1.0 + odd / numerator);
        let delta = denominator * numerator;
        fraction *= delta;
        if (delta - 1.0).abs() < 1e-15 {
            break;
        }
    }
    fraction
}

/// `ln Γ(x)` for `x > 0`, by the Lanczos approximation (g = 7, n = 9).
fn ln_gamma(x: f64) -> f64 {
    const COEFFICIENTS: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        let pi = std::f64::consts::PI;
        return (pi / (pi * x).sin()).ln() - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let mut sum = COEFFICIENTS[0];
    for (index, coefficient) in COEFFICIENTS.iter().enumerate().skip(1) {
        sum += coefficient / (x + index as f64);
    }
    let t = x + 7.5;
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + sum.ln()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::triage::score::upper_bound;

    #[test]
    fn beta_upper_matches_clopper_pearson_at_integer_counts() {
        for (events, trials) in [(0_u16, 10_u16), (0, 300), (1, 474), (5, 1057), (3, 20)] {
            let beta = beta_upper(f64::from(events), f64::from(trials));
            let binomial = upper_bound(usize::from(events), usize::from(trials));
            assert!(
                (beta - binomial).abs() < 1e-6,
                "{events}/{trials}: {beta} vs {binomial}"
            );
        }
        assert!((beta_upper(0.0, 300.0) - 0.009_936).abs() < 1e-4);
    }

    #[test]
    fn korn_graubard_at_zero_events_is_clopper_pearson() {
        assert!((korn_graubard(0.0, 0.0, 300) - upper_bound(0, 300)).abs() < 1e-9);
        // A smaller effective size widens the bound.
        let wide = korn_graubard(0.01, 0.000_1, 1000);
        assert!(wide > upper_bound(10, 1000));
    }

    #[test]
    fn draws_depend_only_on_seed_and_index() {
        let first = SplitMix64::for_draw(7, 3).next_u64();
        assert_eq!(first, SplitMix64::for_draw(7, 3).next_u64());
        assert_ne!(first, SplitMix64::for_draw(7, 4).next_u64());
        assert_ne!(first, SplitMix64::for_draw(8, 3).next_u64());
        // A golden value pins the algorithm.
        assert_eq!(SplitMix64(0).next_u64(), 0xE220_A839_7B1D_CDAF);
    }

    #[test]
    fn the_quantile_index_is_the_ceiling_minus_one() {
        // 10,000 draws: index 9,499, the 9,500th value.
        assert_eq!(10_000 - largest_count(10_000, ALPHA) - 1, 9_499);
        assert_eq!(largest_count(0, ALPHA), 0);
    }

    #[test]
    fn bootstrap_counts_empty_draws_as_rate_one() {
        let samples = [
            Draw {
                batch: 0,
                weight: 1.0,
                passed: true,
                leaked: false,
            },
            Draw {
                batch: 1,
                weight: 1.0,
                passed: false,
                leaked: false,
            },
        ];
        let result = bootstrap(&samples, 2, 1000, 1);
        assert!(result.empty_draws > 0);
        assert!((result.upper - 1.0).abs() < f64::EPSILON);
        // Without a leak and without an empty draw the bound is 0: it cannot
        // bound a zero-event rate alone.
        let clean = [Draw {
            batch: 0,
            weight: 2.0,
            passed: true,
            leaked: false,
        }];
        let result = bootstrap(&clean, 1, 100, 1);
        assert_eq!((result.upper, result.empty_draws), (0.0, 0));
    }
}
