//! The `[signal_processing]` configuration section.

use serde::Deserialize;

use dna_kernel::bounds;
use dna_kernel::error::{ConfigError, Result};

/// Widest run signal context flank: a key's count grows fourfold per call.
const MAX_FLANK: usize = 4;
/// Coarsest and finest run signal share resolutions.
const MIN_RESOLUTION: u16 = 100;
const MAX_RESOLUTION: u16 = 10_000;

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
    /// How a call's second peak is labelled with its cause (SRS-CALL-021).
    pub second_cause: SecondCauseParameters,
    /// Nearest-rank quantile of a read's in-phase residuals above which a
    /// call of a phase mixture is unexplained (SRS-CALL-023).
    pub decompose_unexplained_quantile: f64,
    /// How an electrophoretic spike is found (SRS-SIG-012).
    pub spike: SpikeParameters,
}

/// How an electrophoretic spike is found (SRS-SIG-012): every channel peaks
/// within `span` median call spacings, each at least `height` times the
/// local level, the median primary height of the `level_calls` calls on each
/// side.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpikeParameters {
    /// Multiple of the local level every channel's apex reaches.
    pub height: f64,
    /// Median call spacings within which the four apexes lie.
    pub span: f64,
    /// Calls on each side whose primary heights give the local level.
    pub level_calls: usize,
}

/// How a call's second peak is labelled with its cause (SRS-CALL-021).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SecondCauseParameters {
    /// Below this multiple of the call's lowest channel a second peak is
    /// noise.
    pub noise_ratio: f64,
    /// Offset from the call's centre, in call spacings, from which a second
    /// peak is shifted.
    pub tail_offset: f64,
    /// Multiple of the neighbouring calls' median primary width from which a
    /// centred second peak is a dye blob.
    pub blob_width: f64,
    /// Calls on each side whose primary widths give that median.
    pub width_calls: usize,
}

/// How a run signal is made and how a key's run level is taken from it
/// (PROP-0009): keys of `flank` calls on each side, shares counted in
/// `resolution` steps, and `median + max(margin, mad_factor * MAD)` of the
/// key's shares, when the run has at least `minimum_reads` reads and the key
/// at least `minimum_calls` calls once the sample's own are taken out.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RunSignalParameters {
    /// Calls on each side of a call that form its key's context.
    pub flank: usize,
    /// Steps of a share from 0 to 1 that a run signal counts.
    pub resolution: u16,
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
    pub(crate) run_signal_flank: usize,
    pub(crate) run_signal_resolution: u16,
    pub(crate) run_signal_minimum_reads: usize,
    pub(crate) run_signal_minimum_calls: u64,
    pub(crate) run_signal_margin: f64,
    pub(crate) run_signal_mad_factor: f64,
    pub(crate) second_cause_noise_ratio: f64,
    pub(crate) second_cause_tail_offset: f64,
    pub(crate) second_cause_blob_width: f64,
    pub(crate) second_cause_width_calls: usize,
    pub(crate) decompose_unexplained_quantile: f64,
    pub(crate) spike_height: f64,
    pub(crate) spike_span: f64,
    pub(crate) spike_level_calls: usize,
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
        if !(1..=MAX_FLANK).contains(&self.run_signal_flank) {
            return Err(ConfigError::Constraint(
                "signal_processing.run_signal_flank must be in 1..=4",
            )
            .into());
        }
        if !(MIN_RESOLUTION..=MAX_RESOLUTION).contains(&self.run_signal_resolution) {
            return Err(ConfigError::Constraint(
                "signal_processing.run_signal_resolution must be in 100..=10000",
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
        for (key, value) in [
            (
                "signal_processing.second_cause_noise_ratio",
                self.second_cause_noise_ratio,
            ),
            (
                "signal_processing.second_cause_tail_offset",
                self.second_cause_tail_offset,
            ),
            (
                "signal_processing.second_cause_blob_width",
                self.second_cause_blob_width,
            ),
        ] {
            bounds::positive_finite(key, value)?;
        }
        if !(self.decompose_unexplained_quantile > 0.0 && self.decompose_unexplained_quantile < 1.0)
        {
            return Err(ConfigError::Constraint(
                "signal_processing.decompose_unexplained_quantile must be in (0, 1)",
            )
            .into());
        }
        if self.second_cause_width_calls == 0 {
            return Err(ConfigError::Constraint(
                "signal_processing.second_cause_width_calls must be positive",
            )
            .into());
        }
        if !(self.spike_height.is_finite() && self.spike_height >= 1.0) {
            return Err(ConfigError::Constraint(
                "signal_processing.spike_height must be finite and at least 1",
            )
            .into());
        }
        bounds::positive_finite("signal_processing.spike_span", self.spike_span)?;
        if self.spike_level_calls == 0 {
            return Err(ConfigError::Constraint(
                "signal_processing.spike_level_calls must be positive",
            )
            .into());
        }
        Ok(SignalProcessingConfig {
            window_size_bases: self.window_size_bases,
            minimum_primary_snr: self.minimum_primary_snr,
            minimum_noisy_windows: self.minimum_noisy_windows,
            minimum_noise_sigma: self.minimum_noise_sigma,
            second_floor: self.second_floor,
            run_signal: RunSignalParameters {
                flank: self.run_signal_flank,
                resolution: self.run_signal_resolution,
                minimum_reads: self.run_signal_minimum_reads,
                minimum_calls: self.run_signal_minimum_calls,
                margin: self.run_signal_margin,
                mad_factor: self.run_signal_mad_factor,
            },
            second_cause: SecondCauseParameters {
                noise_ratio: self.second_cause_noise_ratio,
                tail_offset: self.second_cause_tail_offset,
                blob_width: self.second_cause_blob_width,
                width_calls: self.second_cause_width_calls,
            },
            decompose_unexplained_quantile: self.decompose_unexplained_quantile,
            spike: SpikeParameters {
                height: self.spike_height,
                span: self.spike_span,
                level_calls: self.spike_level_calls,
            },
        })
    }
}
