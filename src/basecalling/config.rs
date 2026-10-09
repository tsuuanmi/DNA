//! The `[basecalling]` configuration section.

use serde::Deserialize;

use crate::bounds;
use crate::error::Result;

/// DNA re-calling settings.
#[derive(Debug, Clone)]
pub(crate) struct BasecallingConfig {
    pub(crate) secondary_peak_ratio: f64,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct RawBasecallingConfig {
    pub(crate) secondary_peak_ratio: f64,
}

impl RawBasecallingConfig {
    /// Validates the section.
    pub(crate) fn validate(self) -> Result<BasecallingConfig> {
        bounds::fraction(
            "basecalling.secondary_peak_ratio",
            self.secondary_peak_ratio,
        )?;
        Ok(BasecallingConfig {
            secondary_peak_ratio: self.secondary_peak_ratio,
        })
    }
}
