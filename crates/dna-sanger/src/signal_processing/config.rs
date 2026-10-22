//! The `[signal_processing]` configuration section.

use serde::Deserialize;

use dna_kernel::bounds;
use dna_kernel::error::{ConfigError, Result};

/// Calls per rolling SNR window.
const WINDOW_SIZE_BASES: std::ops::RangeInclusive<usize> = 5..=10;
/// Fewest adjacent noisy windows that can form a noisy region.
const MIN_NOISY_WINDOWS: usize = 2;
/// Run signal context flanks: a key's count grows fourfold per call.
const FLANK: std::ops::RangeInclusive<usize> = 1..=4;
/// Coarsest to finest run signal share resolutions.
const RESOLUTION: std::ops::RangeInclusive<usize> = 100..=10_000;

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
    pub(crate) minimum_noise_sigma: f64,
    /// Lowest second-channel threshold, on the peak scale (PROP-0009).
    pub(crate) second_floor: f64,
    /// How a key's run level is taken from its run's signal (PROP-0009).
    pub run_signal: RunSignalParams,
    /// How a call's second peak is labelled with its cause (SRS-CALL-021).
    pub(crate) second_cause: SecondCauseParams,
    /// Nearest-rank quantile of a read's in-phase residuals above which a
    /// call of a phase mixture is unexplained (SRS-CALL-023).
    pub(crate) decompose_unexplained_quantile: f64,
    /// Fewest in-phase single-channel calls a read needs to measure its own
    /// noise before any phase mixture is fitted (SRS-CALL-023).
    pub(crate) decompose_minimum_calls: usize,
    /// Calls from a dephased segment's start within which the run that
    /// opens its phase mixture must end (SRS-CALL-023).
    pub(crate) decompose_onset_window: usize,
    /// How an electrophoretic spike is found (SRS-SIG-012).
    pub(crate) spike: SpikeParams,
}

/// How an electrophoretic spike is found (SRS-SIG-012): every channel peaks
/// within `span` median call spacings, each at least `height` times the
/// local level, the median primary height of the `level_calls` calls on each
/// side.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpikeParams {
    /// Multiple of the local level every channel's apex reaches.
    pub(crate) height: f64,
    /// Median call spacings within which the four apexes lie.
    pub(crate) span: f64,
    /// Calls on each side whose primary heights give the local level.
    pub(crate) level_calls: usize,
}

/// How a call's second peak is labelled with its cause (SRS-CALL-021).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SecondCauseParams {
    /// Below this multiple of the call's lowest channel a second peak is
    /// noise.
    pub(crate) noise_ratio: f64,
    /// Offset from the call's centre, in call spacings, from which a second
    /// peak is shifted.
    pub(crate) tail_offset: f64,
    /// Multiple of the neighbouring calls' median primary width from which a
    /// centred second peak is a dye blob.
    pub(crate) blob_width: f64,
    /// Calls on each side whose primary widths give that median.
    pub(crate) width_calls: usize,
}

/// How a run signal is made and how a key's run level is taken from it
/// (PROP-0009): keys of `flank` calls on each side, shares counted in
/// `resolution` steps, and `median + max(margin, mad_factor * MAD)` of the
/// key's shares, when the run has at least `minimum_reads` reads and the key
/// at least `minimum_calls` calls once the sample's own are taken out.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RunSignalParams {
    /// Calls on each side of a call that form its key's context.
    pub flank: usize,
    /// Steps of a share from 0 to 1 that a run signal counts.
    pub(crate) resolution: u16,
    /// Fewest reads a run signal needs to be used at all.
    pub(crate) minimum_reads: usize,
    /// Fewest calls a key needs to set a threshold.
    pub(crate) minimum_calls: u64,
    /// Smallest margin above the key's median.
    pub(crate) margin: f64,
    /// Median absolute deviations above the key's median.
    pub(crate) mad_factor: f64,
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
    pub(crate) decompose_minimum_calls: usize,
    pub(crate) decompose_onset_window: usize,
    pub(crate) spike_height: f64,
    pub(crate) spike_span: f64,
    pub(crate) spike_level_calls: usize,
}

impl RawSignalProcessingConfig {
    /// Validates the section.
    pub(crate) fn validate(self) -> Result<SignalProcessingConfig> {
        bounds::count_range(
            "signal_processing.window_size_bases",
            self.window_size_bases,
            WINDOW_SIZE_BASES,
        )?;
        bounds::positive_finite(
            "signal_processing.minimum_primary_snr",
            self.minimum_primary_snr,
        )?;
        bounds::count_at_least(
            "signal_processing.minimum_noisy_windows",
            self.minimum_noisy_windows,
            MIN_NOISY_WINDOWS,
        )?;
        bounds::count_range(
            "signal_processing.run_signal_flank",
            self.run_signal_flank,
            FLANK,
        )?;
        bounds::count_range(
            "signal_processing.run_signal_resolution",
            usize::from(self.run_signal_resolution),
            RESOLUTION,
        )?;
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
        bounds::positive_count(
            "signal_processing.second_cause_width_calls",
            self.second_cause_width_calls,
        )?;
        bounds::positive_count(
            "signal_processing.decompose_minimum_calls",
            self.decompose_minimum_calls,
        )?;
        bounds::positive_count(
            "signal_processing.decompose_onset_window",
            self.decompose_onset_window,
        )?;
        if !(self.spike_height.is_finite() && self.spike_height >= 1.0) {
            return Err(ConfigError::Constraint(
                "signal_processing.spike_height must be finite and at least 1",
            )
            .into());
        }
        bounds::positive_finite("signal_processing.spike_span", self.spike_span)?;
        bounds::positive_count(
            "signal_processing.spike_level_calls",
            self.spike_level_calls,
        )?;
        Ok(SignalProcessingConfig {
            window_size_bases: self.window_size_bases,
            minimum_primary_snr: self.minimum_primary_snr,
            minimum_noisy_windows: self.minimum_noisy_windows,
            minimum_noise_sigma: self.minimum_noise_sigma,
            second_floor: self.second_floor,
            run_signal: RunSignalParams {
                flank: self.run_signal_flank,
                resolution: self.run_signal_resolution,
                minimum_reads: self.run_signal_minimum_reads,
                minimum_calls: self.run_signal_minimum_calls,
                margin: self.run_signal_margin,
                mad_factor: self.run_signal_mad_factor,
            },
            second_cause: SecondCauseParams {
                noise_ratio: self.second_cause_noise_ratio,
                tail_offset: self.second_cause_tail_offset,
                blob_width: self.second_cause_blob_width,
                width_calls: self.second_cause_width_calls,
            },
            decompose_unexplained_quantile: self.decompose_unexplained_quantile,
            decompose_minimum_calls: self.decompose_minimum_calls,
            decompose_onset_window: self.decompose_onset_window,
            spike: SpikeParams {
                height: self.spike_height,
                span: self.spike_span,
                level_calls: self.spike_level_calls,
            },
        })
    }
}

/// The shipped `[signal_processing]` values with `window_size_bases` calls
/// per window and a second floor of 0.21: the tests' configuration.
#[cfg(test)]
pub(crate) fn fixture(window_size_bases: usize) -> SignalProcessingConfig {
    SignalProcessingConfig {
        window_size_bases,
        minimum_primary_snr: 3.0,
        minimum_noisy_windows: 2,
        minimum_noise_sigma: 1.0,
        second_floor: 0.21,
        run_signal: RunSignalParams {
            flank: 2,
            resolution: 1000,
            minimum_reads: 40,
            minimum_calls: 40,
            margin: 0.05,
            mad_factor: 4.0,
        },
        second_cause: SecondCauseParams {
            noise_ratio: 3.0,
            tail_offset: 0.4,
            blob_width: 1.8,
            width_calls: 6,
        },
        decompose_unexplained_quantile: 0.99,
        decompose_minimum_calls: 20,
        decompose_onset_window: 16,
        spike: SpikeParams {
            height: 3.0,
            span: 1.5,
            level_calls: 6,
        },
    }
}

#[cfg(test)]
mod tests {
    use dna_kernel::error::Error;

    use super::RawSignalProcessingConfig;

    /// A change to the section and the message it gives, or `ok`.
    type Case = (fn(&mut RawSignalProcessingConfig), String);

    /// The shipped `[signal_processing]` values.
    fn raw() -> RawSignalProcessingConfig {
        RawSignalProcessingConfig {
            window_size_bases: 10,
            minimum_primary_snr: 3.0,
            minimum_noisy_windows: 2,
            minimum_noise_sigma: 1.0,
            second_floor: 0.20,
            run_signal_flank: 2,
            run_signal_resolution: 1000,
            run_signal_minimum_reads: 40,
            run_signal_minimum_calls: 40,
            run_signal_margin: 0.05,
            run_signal_mad_factor: 4.0,
            second_cause_noise_ratio: 3.0,
            second_cause_tail_offset: 0.4,
            second_cause_blob_width: 1.8,
            second_cause_width_calls: 6,
            decompose_unexplained_quantile: 0.99,
            decompose_minimum_calls: 20,
            decompose_onset_window: 16,
            spike_height: 3.0,
            spike_span: 1.5,
            spike_level_calls: 6,
        }
    }

    /// The message of the section with `change` applied, or `ok`.
    fn message(change: fn(&mut RawSignalProcessingConfig)) -> String {
        let mut config = raw();
        change(&mut config);
        match config.validate() {
            Ok(_) => "ok".to_owned(),
            Err(Error::Config(error)) => error.to_string(),
            Err(other) => format!("not a configuration error: {other}"),
        }
    }

    // The messages are pinned: a refactor of the checks keeps them.
    #[test]
    fn validates_the_signal_processing_section() {
        assert_eq!(message(|_| {}), "ok");
        let positive = |key: &str| format!("signal_processing.{key} must be finite and positive");
        let fraction = |key: &str| {
            format!(
                "signal_processing.{key} must be finite and in [{}, 1]",
                f64::MIN_POSITIVE
            )
        };
        let window = "signal_processing.window_size_bases must be in 5..=10";
        let flank = "signal_processing.run_signal_flank must be in 1..=4";
        let resolution = "signal_processing.run_signal_resolution must be in 100..=10000";
        let minimum = "signal_processing.run_signal_minimum_reads and run_signal_minimum_calls \
                       must be positive";
        let quantile = "signal_processing.decompose_unexplained_quantile must be in (0, 1)";
        let spike = "signal_processing.spike_height must be finite and at least 1";
        let cases: [Case; 30] = [
            (|raw| raw.window_size_bases = 4, window.to_owned()),
            (|raw| raw.window_size_bases = 11, window.to_owned()),
            (|raw| raw.window_size_bases = 5, "ok".to_owned()),
            (
                |raw| raw.minimum_primary_snr = 0.0,
                positive("minimum_primary_snr"),
            ),
            (
                |raw| raw.minimum_noisy_windows = 1,
                "signal_processing.minimum_noisy_windows must be at least 2".to_owned(),
            ),
            (|raw| raw.run_signal_flank = 0, flank.to_owned()),
            (|raw| raw.run_signal_flank = 5, flank.to_owned()),
            (|raw| raw.run_signal_flank = 4, "ok".to_owned()),
            (|raw| raw.run_signal_resolution = 99, resolution.to_owned()),
            (
                |raw| raw.run_signal_resolution = 10_001,
                resolution.to_owned(),
            ),
            (|raw| raw.run_signal_resolution = 10_000, "ok".to_owned()),
            (|raw| raw.run_signal_minimum_reads = 0, minimum.to_owned()),
            (|raw| raw.run_signal_minimum_calls = 0, minimum.to_owned()),
            (
                |raw| raw.minimum_noise_sigma = f64::NAN,
                positive("minimum_noise_sigma"),
            ),
            (|raw| raw.second_floor = 0.0, fraction("second_floor")),
            (|raw| raw.second_floor = 1.01, fraction("second_floor")),
            (
                |raw| raw.run_signal_margin = -0.1,
                fraction("run_signal_margin"),
            ),
            (
                |raw| raw.run_signal_mad_factor = 0.0,
                positive("run_signal_mad_factor"),
            ),
            (
                |raw| raw.second_cause_noise_ratio = 0.0,
                positive("second_cause_noise_ratio"),
            ),
            (
                |raw| raw.second_cause_tail_offset = f64::INFINITY,
                positive("second_cause_tail_offset"),
            ),
            (
                |raw| raw.second_cause_blob_width = -1.0,
                positive("second_cause_blob_width"),
            ),
            (
                |raw| raw.decompose_unexplained_quantile = 0.0,
                quantile.to_owned(),
            ),
            (
                |raw| raw.decompose_unexplained_quantile = 1.0,
                quantile.to_owned(),
            ),
            (
                |raw| raw.decompose_unexplained_quantile = f64::NAN,
                quantile.to_owned(),
            ),
            (
                |raw| raw.second_cause_width_calls = 0,
                "signal_processing.second_cause_width_calls must be positive".to_owned(),
            ),
            (|raw| raw.spike_height = 0.99, spike.to_owned()),
            (|raw| raw.spike_height = f64::INFINITY, spike.to_owned()),
            (|raw| raw.spike_height = 1.0, "ok".to_owned()),
            (|raw| raw.spike_span = 0.0, positive("spike_span")),
            (
                |raw| raw.spike_level_calls = 0,
                "signal_processing.spike_level_calls must be positive".to_owned(),
            ),
        ];
        for (change, expected) in cases {
            assert_eq!(message(change), expected);
        }
    }
}
