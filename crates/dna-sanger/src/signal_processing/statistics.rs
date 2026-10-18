//! Shared local baseline, noise, amplitude, and SNR primitives.

use std::collections::BTreeMap;

use crate::error::{Result, SignalError};
use crate::model::sanger::Chromatogram;

const NORMAL_MAD_SCALE: f64 = 0.674_489_75;
const FIRST_DIFFERENCE_SCALE: f64 = std::f64::consts::SQRT_2;
const OUTPUT_PRECISION: f64 = 1_000_000.0;

/// Baseline and noise of each A/C/G/T channel over one sample span.
#[derive(Debug, Clone, Copy)]
pub(super) struct SpanStatistics {
    pub(super) baselines: [f64; 4],
    pub(super) noise_sigmas: [f64; 4],
}

/// One trace and the channel statistics of its sample spans, each estimated
/// once.
///
/// The windowed features and the locus contexts both span runs of whole
/// locus windows, so most spans recur; a recurring span reuses the first
/// estimate, which is exactly what a fresh estimate would return. The
/// calculators reach the trace only through this record, so its statistics
/// always belong to the trace they read.
pub(super) struct SpanEstimates<'a> {
    trace: &'a Chromatogram,
    /// Floor of every channel's noise scale.
    minimum_noise_sigma: f64,
    estimated: BTreeMap<(usize, usize), SpanStatistics>,
}

impl<'a> SpanEstimates<'a> {
    pub(super) fn new(trace: &'a Chromatogram, minimum_noise_sigma: f64) -> Self {
        Self {
            trace,
            minimum_noise_sigma,
            estimated: BTreeMap::new(),
        }
    }

    /// The trace whose spans this record estimates.
    pub(super) fn trace(&self) -> &'a Chromatogram {
        self.trace
    }

    /// Statistics of the samples `start..end`.
    pub(super) fn span(&mut self, start: usize, end: usize) -> Result<SpanStatistics> {
        if let Some(statistics) = self.estimated.get(&(start, end)) {
            return Ok(*statistics);
        }
        let mut statistics = SpanStatistics {
            baselines: [0.0; 4],
            noise_sigmas: [0.0; 4],
        };
        for (channel, samples) in self.trace.channels.iter().enumerate() {
            let samples = samples
                .get(start..end)
                .ok_or(SignalError::InvalidContextInterval { start, end })?;
            let local = estimate(samples, self.minimum_noise_sigma)?;
            statistics.baselines[channel] = local.baseline;
            statistics.noise_sigmas[channel] = local.noise_sigma;
        }
        self.estimated.insert((start, end), statistics);
        Ok(statistics)
    }
}

/// Local statistics for one analyzed channel sample span.
#[derive(Debug, Clone, Copy)]
struct ChannelStatistics {
    baseline: f64,
    noise_sigma: f64,
}

/// Estimates the median baseline and first-difference-MAD noise scale,
/// floored at `minimum_noise_sigma`.
fn estimate(samples: &[i32], minimum_noise_sigma: f64) -> Result<ChannelStatistics> {
    Ok(ChannelStatistics {
        baseline: median_i32(samples)?,
        noise_sigma: noise_sigma(samples, minimum_noise_sigma)?,
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
pub(crate) fn round_metric(value: f64) -> f64 {
    (value * OUTPUT_PRECISION).round() / OUTPUT_PRECISION
}

pub(super) fn median_usize(values: &[usize]) -> Result<f64> {
    if values.is_empty() {
        return Err(SignalError::InsufficientData(
            "cannot calculate a median from empty integer values",
        )
        .into());
    }
    let mut sorted = values.to_vec();
    sorted.sort_unstable();
    let middle = sorted.len() / 2;
    Ok(if sorted.len().is_multiple_of(2) {
        f64::midpoint(sorted[middle - 1] as f64, sorted[middle] as f64)
    } else {
        sorted[middle] as f64
    })
}

pub(super) fn median_f64(values: &[f64]) -> Result<f64> {
    if values.is_empty() || values.iter().any(|value| !value.is_finite()) {
        return Err(SignalError::InsufficientData(
            "floating-point median requires finite non-empty values",
        )
        .into());
    }
    Ok(median_selected_f64(&mut values.to_vec()))
}

fn median_i32(values: &[i32]) -> Result<f64> {
    if values.is_empty() {
        return Err(SignalError::InsufficientData(
            "cannot calculate a median from an empty sample interval",
        )
        .into());
    }
    // Every `i32` is exact in `f64`, so the median is the integer one.
    let mut values = values.iter().copied().map(f64::from).collect::<Vec<_>>();
    Ok(median_selected_f64(&mut values))
}

fn noise_sigma(samples: &[i32], minimum_noise_sigma: f64) -> Result<f64> {
    if samples.len() < 2 {
        return Err(SignalError::InsufficientData(
            "noise estimation requires at least two channel samples",
        )
        .into());
    }
    let mut differences = samples
        .windows(2)
        .map(|pair| f64::from(pair[1]) - f64::from(pair[0]))
        .collect::<Vec<_>>();
    let center = median_selected_f64(&mut differences);
    for difference in &mut differences {
        *difference = (*difference - center).abs();
    }
    let mad = median_selected_f64(&mut differences);
    Ok((mad / (NORMAL_MAD_SCALE * FIRST_DIFFERENCE_SCALE)).max(minimum_noise_sigma))
}

/// The median of non-empty `values` under `f64::total_cmp`, found by
/// selection instead of a full sort; reorders `values`.
fn median_selected_f64(values: &mut [f64]) -> f64 {
    let even = values.len().is_multiple_of(2);
    let middle = values.len() / 2;
    let (lower, upper, _) = values.select_nth_unstable_by(middle, f64::total_cmp);
    let upper = *upper;
    match lower.iter().copied().max_by(f64::total_cmp) {
        Some(lower) if even => f64::midpoint(lower, upper),
        _ => upper,
    }
}

#[cfg(test)]
mod tests {
    use crate::error::{SangerError as Error, SignalError};
    use crate::model::sanger::VendorEvidence;

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
        assert!(matches!(
            median_f64(&[f64::NAN]),
            Err(Error::Signal(SignalError::InsufficientData(_)))
        ));
        Ok(())
    }

    #[test]
    fn flat_samples_use_the_quantization_floor() -> Result<()> {
        assert_eq!(estimate(&[4, 4, 4, 4], 1.0)?.noise_sigma, 1.0);
        Ok(())
    }

    #[test]
    fn first_difference_mad_uses_the_existing_scale() -> Result<()> {
        let statistics = estimate(&[0, 2, 0, 2, 0], 1.0)?;
        assert_eq!(statistics.baseline, 0.0);
        assert_eq!(
            statistics.noise_sigma,
            2.0 / (NORMAL_MAD_SCALE * FIRST_DIFFERENCE_SCALE)
        );
        Ok(())
    }

    #[test]
    fn a_recurring_span_reuses_the_estimate_a_fresh_one_would_give() -> Result<()> {
        let channels = [
            vec![0, 2, 0, 2, 0, 7],
            vec![1, 1, 1, 1, 1, 1],
            vec![4, 1, 3, 2, 5, 9],
            vec![0, 0, 9, 0, 0, 0],
        ];
        let trace = Chromatogram {
            source_name: "synthetic.ab1".into(),
            source_sha256: String::new(),
            channels: channels.clone(),
            locus_positions: vec![1, 3, 5],
            duplicate_loci: 0,
            vendor: VendorEvidence::default(),
        };
        let mut estimates = SpanEstimates::new(&trace, 1.0);
        let first = estimates.span(0, 5)?;
        let other = estimates.span(1, 6)?;
        let again = estimates.span(0, 5)?;
        assert_eq!(estimates.estimated.len(), 2);
        assert!(matches!(
            estimates.span(4, 7),
            Err(Error::Signal(SignalError::InvalidContextInterval {
                start: 4,
                end: 7
            }))
        ));
        for (channel, samples) in channels.iter().enumerate() {
            for (span, start, end) in [(first, 0, 5), (again, 0, 5), (other, 1, 6)] {
                let fresh = estimate(&samples[start..end], 1.0)?;
                assert_eq!(span.baselines[channel], fresh.baseline);
                assert_eq!(span.noise_sigmas[channel], fresh.noise_sigma);
            }
        }
        Ok(())
    }

    #[test]
    fn corrects_and_rounds_finite_metrics() {
        assert_eq!(corrected_amplitude(-7, 0.0), 0.0);
        assert_eq!(snr(3.0, 1.0), 3.0);
        assert_eq!(round_metric(1.234_567_89), 1.234_568);
    }
}
