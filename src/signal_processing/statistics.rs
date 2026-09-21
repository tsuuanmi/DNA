//! Shared local baseline, noise, amplitude, and SNR primitives.

use crate::error::{Error, Result};

const NORMAL_MAD_SCALE: f64 = 0.674_489_75;
const FIRST_DIFFERENCE_SCALE: f64 = std::f64::consts::SQRT_2;
const MINIMUM_NOISE_SIGMA: f64 = 1.0;
const OUTPUT_PRECISION: f64 = 1_000_000.0;

/// Local statistics for one analyzed channel sample span.
#[derive(Debug, Clone, Copy)]
pub(super) struct ChannelStatistics {
    pub(super) baseline: f64,
    pub(super) noise_sigma: f64,
}

/// Estimates the median baseline and first-difference-MAD noise scale.
pub(super) fn estimate(samples: &[i32]) -> Result<ChannelStatistics> {
    Ok(ChannelStatistics {
        baseline: median_i32(samples)?,
        noise_sigma: noise_sigma(samples)?,
    })
}

/// Returns non-negative amplitude above the local baseline.
pub(super) fn corrected_amplitude(raw_height: i32, baseline: f64) -> f64 {
    (f64::from(raw_height) - baseline).max(0.0)
}

/// Returns local signal-to-noise ratio from a corrected amplitude.
pub(super) fn snr(corrected_amplitude: f64, noise_sigma: f64) -> f64 {
    corrected_amplitude / noise_sigma
}

/// Rounds an emitted or thresholded signal metric to six decimal places.
pub(super) fn round_metric(value: f64) -> f64 {
    (value * OUTPUT_PRECISION).round() / OUTPUT_PRECISION
}

pub(super) fn median_usize(values: &[usize]) -> Result<f64> {
    if values.is_empty() {
        return Err(Error::DNAProcessing(
            "cannot calculate a median from empty integer values".into(),
        ));
    }
    let mut sorted = values.to_vec();
    sorted.sort_unstable();
    let middle = sorted.len() / 2;
    Ok(if sorted.len().is_multiple_of(2) {
        (sorted[middle - 1] as f64 + sorted[middle] as f64) / 2.0
    } else {
        sorted[middle] as f64
    })
}

pub(super) fn median_f64(values: &[f64]) -> Result<f64> {
    if values.is_empty() || values.iter().any(|value| !value.is_finite()) {
        return Err(Error::DNAProcessing(
            "floating-point median requires finite non-empty values".into(),
        ));
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    Ok(median_sorted_f64(&sorted))
}

fn median_i32(values: &[i32]) -> Result<f64> {
    if values.is_empty() {
        return Err(Error::DNAProcessing(
            "cannot calculate a median from an empty sample interval".into(),
        ));
    }
    let mut sorted = values.to_vec();
    sorted.sort_unstable();
    Ok(median_sorted_i32(&sorted))
}

fn median_sorted_i32(sorted: &[i32]) -> f64 {
    let middle = sorted.len() / 2;
    if sorted.len().is_multiple_of(2) {
        (f64::from(sorted[middle - 1]) + f64::from(sorted[middle])) / 2.0
    } else {
        f64::from(sorted[middle])
    }
}

fn noise_sigma(samples: &[i32]) -> Result<f64> {
    if samples.len() < 2 {
        return Err(Error::DNAProcessing(
            "noise estimation requires at least two channel samples".into(),
        ));
    }
    let mut differences = samples
        .windows(2)
        .map(|pair| f64::from(pair[1]) - f64::from(pair[0]))
        .collect::<Vec<_>>();
    differences.sort_by(f64::total_cmp);
    let center = median_sorted_f64(&differences);
    let mut deviations = differences
        .into_iter()
        .map(|difference| (difference - center).abs())
        .collect::<Vec<_>>();
    deviations.sort_by(f64::total_cmp);
    let mad = median_sorted_f64(&deviations);
    Ok((mad / (NORMAL_MAD_SCALE * FIRST_DIFFERENCE_SCALE)).max(MINIMUM_NOISE_SIGMA))
}

fn median_sorted_f64(sorted: &[f64]) -> f64 {
    let middle = sorted.len() / 2;
    if sorted.len().is_multiple_of(2) {
        (sorted[middle - 1] + sorted[middle]) / 2.0
    } else {
        sorted[middle]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn median_handles_even_and_odd_sample_counts() -> Result<()> {
        assert_eq!(median_i32(&[4, 1, 3])?, 3.0);
        assert_eq!(median_i32(&[4, 1, 3, 2])?, 2.5);
        Ok(())
    }

    #[test]
    fn shared_medians_handle_integer_and_float_inputs() -> Result<()> {
        assert_eq!(median_usize(&[2, 8, 4])?, 4.0);
        assert_eq!(median_usize(&[2, 8, 4, 6])?, 5.0);
        assert_eq!(median_f64(&[1.0, 9.0, 3.0])?, 3.0);
        assert_eq!(median_f64(&[1.0, 9.0, 3.0, 5.0])?, 4.0);
        assert!(median_f64(&[f64::NAN]).is_err());
        Ok(())
    }

    #[test]
    fn flat_samples_use_the_quantization_floor() -> Result<()> {
        assert_eq!(estimate(&[4, 4, 4, 4])?.noise_sigma, 1.0);
        Ok(())
    }

    #[test]
    fn first_difference_mad_uses_the_existing_scale() -> Result<()> {
        let statistics = estimate(&[0, 2, 0, 2, 0])?;
        assert_eq!(statistics.baseline, 0.0);
        assert_eq!(
            statistics.noise_sigma,
            2.0 / (NORMAL_MAD_SCALE * FIRST_DIFFERENCE_SCALE)
        );
        Ok(())
    }

    #[test]
    fn corrects_and_rounds_finite_metrics() {
        assert_eq!(corrected_amplitude(-7, 0.0), 0.0);
        assert_eq!(snr(3.0, 1.0), 3.0);
        assert_eq!(round_metric(1.234_567_89), 1.234_568);
    }
}
