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
    /// Calls after a repeat run that form its prior window.
    pub(crate) repeat_window_calls: usize,
    /// Secondary ratio at or above which a position carries a double peak.
    pub(crate) double_peak_ratio: f64,
    /// Relative spacing deviation at or above which a position is a spacing
    /// defect.
    pub(crate) spacing_defect: f64,
    /// Spacings on each side of a position whose median anchors its spacing
    /// deviation.
    pub(crate) spacing_half_window: usize,
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
    pub(crate) repeat_window_calls: usize,
    pub(crate) double_peak_ratio: f64,
    pub(crate) spacing_defect: f64,
    pub(crate) spacing_half_window: usize,
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
        if self.repeat_window_calls == 0 || self.spacing_half_window == 0 {
            return Err(ConfigError::Constraint(
                "callability.repeat_window_calls and spacing_half_window must be positive",
            )
            .into());
        }
        bounds::fraction("callability.double_peak_ratio", self.double_peak_ratio)?;
        bounds::positive_finite("callability.spacing_defect", self.spacing_defect)?;
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
            repeat_window_calls: self.repeat_window_calls,
            double_peak_ratio: self.double_peak_ratio,
            spacing_defect: self.spacing_defect,
            spacing_half_window: self.spacing_half_window,
        })
    }
}

#[cfg(test)]
mod tests {
    use dna_kernel::error::Error;

    use super::RawCallabilityConfig;

    /// A change to the section and the message it gives, or `ok`.
    type Case = (fn(&mut RawCallabilityConfig), String);

    /// The shipped `[callability]` values.
    fn raw() -> RawCallabilityConfig {
        RawCallabilityConfig {
            window_calls: 16,
            onset_defect_fraction: 0.375,
            exit_defect_fraction: 0.125,
            minimum_main_share: 0.35,
            maximum_far_share: 0.12,
            minimum_shadow_share: 0.1,
            weak_amplitude_fraction: 0.2,
            repeat_min_length: 8,
            minimum_callable_calls: 20,
            repeat_window_calls: 16,
            double_peak_ratio: 0.33,
            spacing_defect: 0.5,
            spacing_half_window: 8,
        }
    }

    /// The message of the section with `change` applied, or `ok`.
    fn message(change: fn(&mut RawCallabilityConfig)) -> String {
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
    fn validates_the_callability_section() {
        assert_eq!(message(|_| {}), "ok");
        let fraction = |key: &str| {
            format!(
                "callability.{key} must be finite and in [{}, 1]",
                f64::MIN_POSITIVE
            )
        };
        let closed = |key: &str, maximum: &str| {
            format!("callability.{key} must be finite and in [0, {maximum}]")
        };
        let window = "callability.window_calls must be in 8..=64";
        let positive = "callability.repeat_window_calls and spacing_half_window must be positive";
        let spacing = "callability.spacing_defect must be finite and positive";
        let cases: [Case; 24] = [
            (|raw| raw.window_calls = 7, window.to_owned()),
            (|raw| raw.window_calls = 65, window.to_owned()),
            (|raw| raw.window_calls = 8, "ok".to_owned()),
            (|raw| raw.window_calls = 64, "ok".to_owned()),
            (
                |raw| raw.onset_defect_fraction = 0.0,
                fraction("onset_defect_fraction"),
            ),
            (
                |raw| raw.exit_defect_fraction = -0.1,
                closed("exit_defect_fraction", "1"),
            ),
            (|raw| raw.exit_defect_fraction = 0.0, "ok".to_owned()),
            (
                |raw| raw.exit_defect_fraction = 0.375,
                "callability.exit_defect_fraction must be less than \
                 callability.onset_defect_fraction"
                    .to_owned(),
            ),
            (
                |raw| raw.minimum_main_share = 0.0,
                fraction("minimum_main_share"),
            ),
            (
                |raw| raw.maximum_far_share = 1.5,
                closed("maximum_far_share", "1"),
            ),
            (|raw| raw.maximum_far_share = 0.0, "ok".to_owned()),
            (
                |raw| raw.minimum_shadow_share = f64::NAN,
                fraction("minimum_shadow_share"),
            ),
            (
                |raw| raw.weak_amplitude_fraction = 0.51,
                closed("weak_amplitude_fraction", "0.5"),
            ),
            (|raw| raw.weak_amplitude_fraction = 0.5, "ok".to_owned()),
            (
                |raw| raw.repeat_min_length = 1,
                "callability.repeat_min_length must be at least 2".to_owned(),
            ),
            (|raw| raw.repeat_min_length = 2, "ok".to_owned()),
            (
                |raw| raw.minimum_callable_calls = 0,
                "callability.minimum_callable_calls must be positive".to_owned(),
            ),
            (|raw| raw.repeat_window_calls = 0, positive.to_owned()),
            (|raw| raw.spacing_half_window = 0, positive.to_owned()),
            (
                |raw| raw.double_peak_ratio = 0.0,
                fraction("double_peak_ratio"),
            ),
            (|raw| raw.double_peak_ratio = 1.0, "ok".to_owned()),
            (|raw| raw.spacing_defect = 0.0, spacing.to_owned()),
            (|raw| raw.spacing_defect = f64::INFINITY, spacing.to_owned()),
            (|raw| raw.spacing_defect = 2.0, "ok".to_owned()),
        ];
        for (change, expected) in cases {
            assert_eq!(message(change), expected);
        }
    }
}
