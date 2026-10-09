//! The `[signal_processing]` configuration section.

use serde::Deserialize;

use crate::bounds;
use crate::error::{ConfigError, Result};

/// Observation-only rolling signal-quality settings.
#[derive(Debug, Clone)]
pub(crate) struct SignalProcessingConfig {
    pub(crate) window_size_bases: usize,
    pub(crate) minimum_primary_snr: f64,
    pub(crate) minimum_noisy_windows: usize,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct RawSignalProcessingConfig {
    pub(crate) window_size_bases: usize,
    pub(crate) minimum_primary_snr: f64,
    pub(crate) minimum_noisy_windows: usize,
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
        Ok(SignalProcessingConfig {
            window_size_bases: self.window_size_bases,
            minimum_primary_snr: self.minimum_primary_snr,
            minimum_noisy_windows: self.minimum_noisy_windows,
        })
    }
}
