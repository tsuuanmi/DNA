//! The `[callability]` configuration section.

use serde::Deserialize;

use dna_kernel::bounds;
use dna_kernel::error::{ConfigError, Result};

/// Signal-derived read-callability settings.
#[derive(Debug, Clone)]
pub struct CallabilityConfig {
    /// Calls per rolling window of the phase statistics.
    pub(crate) window_calls: usize,
    /// Defect fraction at or above which an in-phase stretch ends.
    pub(crate) onset_defect_fraction: f64,
    /// Defect fraction at or below which a masked stretch ends.
    pub(crate) exit_defect_fraction: f64,
    /// Smallest main-ladder share of a segment's shadow fit for the segment to
    /// count as dephased rather than mixed.
    pub(crate) minimum_main_share: f64,
    /// Largest far-shadow (offsets of two or three calls) share of a dephased
    /// segment's shadow fit.
    pub(crate) maximum_far_share: f64,
    /// Smallest share for a shadow offset to be reported; a dephased segment
    /// needs a one-call shadow at or above it.
    pub(crate) minimum_shadow_share: f64,
    /// Fraction of the read's median primary amplitude below which a position is weak.
    pub(crate) weak_amplitude_fraction: f64,
    /// Shortest run of identical or alternating primary calls treated as a repeat.
    pub(crate) repeat_min_length: usize,
    /// Fewest callable calls a read needs to be analyzed.
    pub(crate) minimum_callable_calls: usize,
}

/// The `[callability]` section as written in the configuration.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RawCallabilityConfig {
    pub(crate) window_calls: usize,
    pub(crate) onset_defect_fraction: f64,
    pub(crate) exit_defect_fraction: f64,
    pub(crate) minimum_main_share: f64,
    pub(crate) maximum_far_share: f64,
    pub(crate) minimum_shadow_share: f64,
    pub(crate) weak_amplitude_fraction: f64,
    pub(crate) repeat_min_length: usize,
    pub(crate) minimum_callable_calls: usize,
}

impl RawCallabilityConfig {
    /// Validates the section.
    pub(crate) fn validate(self) -> Result<CallabilityConfig> {
        if !(8..=64).contains(&self.window_calls) {
            return Err(
                ConfigError::Constraint("callability.window_calls must be in 8..=64").into(),
            );
        }
        bounds::fraction(
            "callability.onset_defect_fraction",
            self.onset_defect_fraction,
        )?;
        bounds::finite_range(
            "callability.exit_defect_fraction",
            self.exit_defect_fraction,
            0.0,
            1.0,
        )?;
        if self.exit_defect_fraction >= self.onset_defect_fraction {
            return Err(ConfigError::Constraint(
                "callability.exit_defect_fraction must be less than callability.onset_defect_fraction",
            )
            .into());
        }
        bounds::fraction("callability.minimum_main_share", self.minimum_main_share)?;
        bounds::finite_range(
            "callability.maximum_far_share",
            self.maximum_far_share,
            0.0,
            1.0,
        )?;
        bounds::fraction(
            "callability.minimum_shadow_share",
            self.minimum_shadow_share,
        )?;
        bounds::finite_range(
            "callability.weak_amplitude_fraction",
            self.weak_amplitude_fraction,
            0.0,
            0.5,
        )?;
        if self.repeat_min_length < 2 {
            return Err(ConfigError::Constraint(
                "callability.repeat_min_length must be at least 2",
            )
            .into());
        }
        if self.minimum_callable_calls == 0 {
            return Err(ConfigError::Constraint(
                "callability.minimum_callable_calls must be positive",
            )
            .into());
        }
        Ok(CallabilityConfig {
            window_calls: self.window_calls,
            onset_defect_fraction: self.onset_defect_fraction,
            exit_defect_fraction: self.exit_defect_fraction,
            minimum_main_share: self.minimum_main_share,
            maximum_far_share: self.maximum_far_share,
            minimum_shadow_share: self.minimum_shadow_share,
            weak_amplitude_fraction: self.weak_amplitude_fraction,
            repeat_min_length: self.repeat_min_length,
            minimum_callable_calls: self.minimum_callable_calls,
        })
    }
}
