//! The Sanger plugin's configuration: its `[sanger_evidence]` section and
//! every section it owns, validated together (ADR-0069).

use serde::Deserialize;

use dna_kernel::bounds;
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

impl SangerConfig {
    /// Fewest callable calls a read needs to be analyzed
    /// (`callability.minimum_callable_calls`).
    #[must_use]
    pub const fn minimum_callable_calls(&self) -> usize {
        self.callability.minimum_callable_calls
    }
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
/// support vetoes (ADR-0069), and the parameters of its read-edge rules
/// (SRS-CALL-013, SRS-CALL-014, SRS-CALL-016 to SRS-CALL-018).
#[derive(Debug, Clone)]
pub struct SangerEvidenceConfig {
    /// Smallest highest A/C/G/T peak of a supporting call.
    pub minimum_peak_height: i32,
    /// Relative quality a supporting call must exceed.
    pub relative_quality_threshold: u8,
    /// Relative quality at or above which a call below the peak floor still
    /// shows its base; its indel stays withheld.
    pub peak_rescue_quality: u8,
    /// A rescued call's second base must hold less than this peak share.
    pub peak_rescue_second_share: f64,
    /// Calls from an end of the usable signal that lie at a signal edge
    /// (SRS-CALL-013).
    pub signal_edge_margin: usize,
    /// Smallest primary peak of a leading call vouched for by its base, as a
    /// fraction of the read's median primary peak (SRS-CALL-014).
    pub leading_height_fraction: f64,
    /// Calls from an informative end within which the run before a
    /// post-repeat phase loss counts as bounded (SRS-CALL-016).
    pub post_repeat_margin: usize,
    /// Fewest informative calls between a run's bound and a post-repeat edge
    /// (SRS-CALL-016).
    pub post_repeat_minimum_calls: usize,
    /// Calls from an informative end within which an edge call may be
    /// vouched for by its own signal (SRS-CALL-017).
    pub edge_identity_margin: usize,
    /// Smallest primary peak of such a call, as a fraction of the read's
    /// median (SRS-CALL-017).
    pub edge_height_fraction: f64,
    /// Second-channel share below which such a call's signal is clean
    /// (SRS-CALL-017, SRS-CALL-018).
    pub edge_second_share: f64,
    /// Widest primary peak of such a call, as a multiple of the median width
    /// (SRS-CALL-017, SRS-CALL-018).
    pub edge_width_factor: f64,
    /// Narrowest spacing to a neighbouring call, as a multiple of the median
    /// spacing (SRS-CALL-017, SRS-CALL-018).
    pub edge_spacing_minimum: f64,
    /// Widest such spacing (SRS-CALL-017, SRS-CALL-018).
    pub edge_spacing_maximum: f64,
    /// Calls past the callable span's start within which a start call may be
    /// vouched for by its own signal (SRS-CALL-018).
    pub start_identity_margin: usize,
    /// First in-phase single-channel calls of the callable span whose medians
    /// are the start's norms (SRS-CALL-018).
    pub start_norm_calls: usize,
    /// Smallest primary peak of a start call, as a fraction of the start's
    /// median (SRS-CALL-018).
    pub start_height_fraction: f64,
    /// Second-channel share below which a start call's signal is clean
    /// (SRS-CALL-018).
    pub start_second_share: f64,
}

/// The `[sanger_evidence]` section as written in the configuration.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RawSangerEvidenceConfig {
    pub(crate) minimum_peak_height: i32,
    pub(crate) relative_quality_threshold: u8,
    pub(crate) peak_rescue_quality: u8,
    pub(crate) peak_rescue_second_share: f64,
    pub(crate) signal_edge_margin: usize,
    pub(crate) leading_height_fraction: f64,
    pub(crate) post_repeat_margin: usize,
    pub(crate) post_repeat_minimum_calls: usize,
    pub(crate) edge_identity_margin: usize,
    pub(crate) edge_height_fraction: f64,
    pub(crate) edge_second_share: f64,
    pub(crate) edge_width_factor: f64,
    pub(crate) edge_spacing_minimum: f64,
    pub(crate) edge_spacing_maximum: f64,
    pub(crate) start_identity_margin: usize,
    pub(crate) start_norm_calls: usize,
    pub(crate) start_height_fraction: f64,
    pub(crate) start_second_share: f64,
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
        for (key, value) in [
            (
                "sanger_evidence.leading_height_fraction",
                self.leading_height_fraction,
            ),
            (
                "sanger_evidence.edge_height_fraction",
                self.edge_height_fraction,
            ),
            ("sanger_evidence.edge_width_factor", self.edge_width_factor),
            (
                "sanger_evidence.edge_spacing_minimum",
                self.edge_spacing_minimum,
            ),
            (
                "sanger_evidence.edge_spacing_maximum",
                self.edge_spacing_maximum,
            ),
            (
                "sanger_evidence.start_height_fraction",
                self.start_height_fraction,
            ),
        ] {
            bounds::positive_finite(key, value)?;
        }
        bounds::fraction("sanger_evidence.edge_second_share", self.edge_second_share)?;
        bounds::fraction(
            "sanger_evidence.peak_rescue_second_share",
            self.peak_rescue_second_share,
        )?;
        bounds::fraction(
            "sanger_evidence.start_second_share",
            self.start_second_share,
        )?;
        if self.edge_spacing_minimum > self.edge_spacing_maximum {
            return Err(ConfigError::Constraint(
                "sanger_evidence.edge_spacing_minimum must not exceed edge_spacing_maximum",
            )
            .into());
        }
        if self.peak_rescue_quality <= self.relative_quality_threshold {
            return Err(ConfigError::Constraint(
                "sanger_evidence.peak_rescue_quality must exceed relative_quality_threshold",
            )
            .into());
        }
        if self.post_repeat_minimum_calls == 0 || self.start_norm_calls == 0 {
            return Err(ConfigError::Constraint(
                "sanger_evidence.post_repeat_minimum_calls and start_norm_calls must be positive",
            )
            .into());
        }
        Ok(SangerEvidenceConfig {
            minimum_peak_height: self.minimum_peak_height,
            relative_quality_threshold: self.relative_quality_threshold,
            peak_rescue_quality: self.peak_rescue_quality,
            peak_rescue_second_share: self.peak_rescue_second_share,
            signal_edge_margin: self.signal_edge_margin,
            leading_height_fraction: self.leading_height_fraction,
            post_repeat_margin: self.post_repeat_margin,
            post_repeat_minimum_calls: self.post_repeat_minimum_calls,
            edge_identity_margin: self.edge_identity_margin,
            edge_height_fraction: self.edge_height_fraction,
            edge_second_share: self.edge_second_share,
            edge_width_factor: self.edge_width_factor,
            edge_spacing_minimum: self.edge_spacing_minimum,
            edge_spacing_maximum: self.edge_spacing_maximum,
            start_identity_margin: self.start_identity_margin,
            start_norm_calls: self.start_norm_calls,
            start_height_fraction: self.start_height_fraction,
            start_second_share: self.start_second_share,
        })
    }
}
