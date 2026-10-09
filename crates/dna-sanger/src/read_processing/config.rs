//! The Sanger plugin's configuration: its `[sanger_evidence]` section and
//! every section it owns, validated together (ADR-0069).

use serde::Deserialize;

use dna_kernel::error::{ConfigError, Result};

use crate::basecalling::{BasecallingConfig, RawBasecallingConfig};
use crate::callability::{CallabilityConfig, RawCallabilityConfig};
use crate::quality_control::{QualityControlConfig, RawQualityControlConfig};
use crate::signal_processing::{RawSignalProcessingConfig, SignalProcessingConfig};

/// Largest peak height representable by an ABIF signed short.
pub(crate) const MAX_PEAK_HEIGHT: i32 = i16::MAX as i32;

/// Every configuration section the Sanger plugin owns.
#[derive(Debug, Clone)]
pub struct SangerConfig {
    pub(crate) basecalling: BasecallingConfig,
    /// The `[signal_processing]` section.
    pub signal_processing: SignalProcessingConfig,
    pub(crate) callability: CallabilityConfig,
    pub(crate) quality_control: QualityControlConfig,
    /// The `[sanger_evidence]` section.
    pub sanger_evidence: SangerEvidenceConfig,
}

/// The Sanger plugin's sections as written in the configuration.
pub struct RawSangerConfig {
    /// The `[basecalling]` section as written.
    pub basecalling: RawBasecallingConfig,
    /// The `[signal_processing]` section as written.
    pub signal_processing: RawSignalProcessingConfig,
    /// The `[callability]` section as written.
    pub callability: RawCallabilityConfig,
    /// The `[quality_control]` section as written.
    pub quality_control: RawQualityControlConfig,
    /// The `[sanger_evidence]` section as written.
    pub sanger_evidence: RawSangerEvidenceConfig,
}

impl RawSangerConfig {
    /// Validates every section, then the checks that span sections.
    ///
    /// # Errors
    ///
    /// Returns `ConfigError` naming the first section value or cross-section
    /// check that fails.
    pub fn validate(self) -> Result<SangerConfig> {
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
pub struct SangerEvidenceConfig {
    /// Smallest highest A/C/G/T peak of a supporting call.
    pub minimum_peak_height: i32,
    /// Relative quality a supporting call must exceed.
    pub relative_quality_threshold: u8,
}

/// The `[sanger_evidence]` section as written in the configuration.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RawSangerEvidenceConfig {
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
