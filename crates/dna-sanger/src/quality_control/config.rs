//! The `[quality_control]` configuration section.

use serde::Deserialize;

use dna_kernel::bounds;
use dna_kernel::error::Result;

/// Relative quality settings.
#[derive(Debug, Clone)]
pub struct QualityControlConfig {
    pub(crate) penalty_window_size: usize,
    pub(crate) max_relative_quality_score: u8,
    /// Most calls of an adjacent dephased segment the trim keeps beyond the
    /// callable span (SRS-QC-004).
    pub(crate) trim_context_margin: usize,
}

/// The `[quality_control]` section as written in the configuration.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RawQualityControlConfig {
    pub(crate) penalty_window_size: usize,
    pub(crate) max_relative_quality_score: u8,
    pub(crate) trim_context_margin: usize,
}

impl RawQualityControlConfig {
    /// Validates the section.
    pub(crate) fn validate(self) -> Result<QualityControlConfig> {
        bounds::positive_count(
            "quality_control.penalty_window_size",
            self.penalty_window_size,
        )?;
        bounds::positive_count(
            "quality_control.max_relative_quality_score",
            usize::from(self.max_relative_quality_score),
        )?;
        Ok(QualityControlConfig {
            penalty_window_size: self.penalty_window_size,
            max_relative_quality_score: self.max_relative_quality_score,
            trim_context_margin: self.trim_context_margin,
        })
    }
}
