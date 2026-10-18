//! The `[signal_processing]` configuration section.

use serde::Deserialize;

use dna_kernel::bounds;
use dna_kernel::error::{ConfigError, Result};

/// Observation-only rolling signal-quality settings.
#[derive(Debug, Clone)]
pub struct SignalProcessingConfig {
    /// Calls per rolling SNR window.
    pub window_size_bases: usize,
    /// Primary-channel SNR below which a window is noisy.
    pub minimum_primary_snr: f64,
    /// Fewest adjacent noisy windows that form a noisy region.
    pub minimum_noisy_windows: usize,
    /// Floor of a channel's noise scale, in signal units.
    pub minimum_noise_sigma: f64,
    /// Lowest second-channel threshold, on the peak scale (PROP-0009).
    pub second_floor: f64,
    /// How a key's run level is taken from its run's signal (PROP-0009).
    pub run_signal: RunSignalParameters,
}

/// How a key's run level is taken from its run's signal (PROP-0009):
/// `median + max(margin, mad_factor * MAD)` of the key's shares, when the run
/// has at least `minimum_reads` reads and the key at least `minimum_calls`
/// calls once the sample's own are taken out.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RunSignalParameters {
    /// Fewest reads a run signal needs to be used at all.
    pub minimum_reads: usize,
    /// Fewest calls a key needs to set a threshold.
    pub minimum_calls: u64,
    /// Smallest margin above the key's median.
    pub margin: f64,
    /// Median absolute deviations above the key's median.
    pub mad_factor: f64,
}

/// The `[signal_processing]` section as written in the configuration.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RawSignalProcessingConfig {
    pub(crate) window_size_bases: usize,
    pub(crate) minimum_primary_snr: f64,
    pub(crate) minimum_noisy_windows: usize,
    pub(crate) minimum_noise_sigma: f64,
    pub(crate) second_floor: f64,
    pub(crate) run_signal_minimum_reads: usize,
    pub(crate) run_signal_minimum_calls: u64,
    pub(crate) run_signal_margin: f64,
    pub(crate) run_signal_mad_factor: f64,
}

impl RawSignalProcessingConfig {
    /// Validates the section.
    pub(crate) fn validate(self) -> Result<SignalProcessingConfig> {
        if !(5..=10).contains(&self.window_size_bases) {
            return Err(ConfigError::Constraint(
                "signal_processing.window_size_bases must be in 5..=10",
            )
            .into());
        }
        bounds::positive_finite(
            "signal_processing.minimum_primary_snr",
            self.minimum_primary_snr,
        )?;
        if self.minimum_noisy_windows < 2 {
            return Err(ConfigError::Constraint(
                "signal_processing.minimum_noisy_windows must be at least 2",
            )
            .into());
        }
        if self.run_signal_minimum_reads == 0 || self.run_signal_minimum_calls == 0 {
            return Err(ConfigError::Constraint(
                "signal_processing.run_signal_minimum_reads and run_signal_minimum_calls must be positive",
            )
            .into());
        }
        bounds::positive_finite(
            "signal_processing.minimum_noise_sigma",
            self.minimum_noise_sigma,
        )?;
        bounds::fraction("signal_processing.second_floor", self.second_floor)?;
        bounds::fraction(
            "signal_processing.run_signal_margin",
            self.run_signal_margin,
        )?;
        bounds::positive_finite(
            "signal_processing.run_signal_mad_factor",
            self.run_signal_mad_factor,
        )?;
        Ok(SignalProcessingConfig {
            window_size_bases: self.window_size_bases,
            minimum_primary_snr: self.minimum_primary_snr,
            minimum_noisy_windows: self.minimum_noisy_windows,
            minimum_noise_sigma: self.minimum_noise_sigma,
            second_floor: self.second_floor,
            run_signal: RunSignalParameters {
                minimum_reads: self.run_signal_minimum_reads,
                minimum_calls: self.run_signal_minimum_calls,
                margin: self.run_signal_margin,
                mad_factor: self.run_signal_mad_factor,
            },
        })
    }
}
