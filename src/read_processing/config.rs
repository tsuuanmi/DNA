//! The Sanger plugin's configuration: its `[sanger_evidence]` section and
//! every section it owns, validated together (ADR-0069).

use serde::Deserialize;

use crate::error::{ConfigError, Result};

use crate::basecalling::{BasecallingConfig, RawBasecallingConfig};
use crate::callability::{CallabilityConfig, RawCallabilityConfig};
use crate::quality_control::{QualityControlConfig, RawQualityControlConfig};
use crate::signal_processing::{RawSignalProcessingConfig, SignalProcessingConfig};

/// Largest peak height representable by an ABIF signed short.
pub(crate) const MAX_PEAK_HEIGHT: i32 = i16::MAX as i32;

/// Every configuration section the Sanger plugin owns.
#[derive(Debug, Clone)]
pub(crate) struct SangerConfig {
    pub(crate) basecalling: BasecallingConfig,
    pub(crate) signal_processing: SignalProcessingConfig,
    pub(crate) callability: CallabilityConfig,
    pub(crate) quality_control: QualityControlConfig,
    pub(crate) sanger_evidence: SangerEvidenceConfig,
}

/// The Sanger plugin's sections as written in the configuration.
pub(crate) struct RawSangerConfig {
    pub(crate) basecalling: RawBasecallingConfig,
    pub(crate) signal_processing: RawSignalProcessingConfig,
    pub(crate) callability: RawCallabilityConfig,
    pub(crate) quality_control: RawQualityControlConfig,
    pub(crate) sanger_evidence: RawSangerEvidenceConfig,
}

impl RawSangerConfig {
    /// Validates every section, then the checks that span sections.
    pub(crate) fn validate(self) -> Result<SangerConfig> {
        let basecalling = self.basecalling.validate()?;
        let signal_processing = self.signal_processing.validate()?;
        let callability = self.callability.validate()?;
        let quality_control = self.quality_control.validate()?;
        let sanger_evidence = self.sanger_evidence.validate()?;
        if sanger_evidence.relative_quality_threshold >= quality_control.max_relative_quality_score
        {
            return Err(ConfigError::Constraint(
                "sanger_evidence.relative_quality_threshold must be less than quality_control.max_relative_quality_score",
            )
            .into());
        }
        Ok(SangerConfig {
            basecalling,
            signal_processing,
            callability,
            quality_control,
            sanger_evidence,
        })
    }
}

/// Sanger support thresholds that the Sanger evidence adapter turns into
/// support vetoes (ADR-0069).
#[derive(Debug, Clone)]
pub(crate) struct SangerEvidenceConfig {
    pub(crate) minimum_peak_height: i32,
    pub(crate) relative_quality_threshold: u8,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct RawSangerEvidenceConfig {
    pub(crate) minimum_peak_height: i32,
    pub(crate) relative_quality_threshold: u8,
}

impl RawSangerEvidenceConfig {
    /// Validates the section.
    pub(crate) fn validate(self) -> Result<SangerEvidenceConfig> {
        if self.minimum_peak_height <= 0 || self.minimum_peak_height > MAX_PEAK_HEIGHT {
            return Err(ConfigError::MinimumPeakHeight {
                maximum: MAX_PEAK_HEIGHT,
            }
            .into());
        }
        Ok(SangerEvidenceConfig {
            minimum_peak_height: self.minimum_peak_height,
            relative_quality_threshold: self.relative_quality_threshold,
        })
    }
}
