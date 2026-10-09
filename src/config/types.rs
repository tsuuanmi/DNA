//! Typed and validated configuration records.

use std::path::{Path, PathBuf};

use serde::Deserialize;

use crate::config::defaults::{MAX_INDEL_LENGTH, MAX_PEAK_HEIGHT};
use crate::error::{ConfigError, Result};

/// Configuration schema version this build accepts.
const SCHEMA_VERSION: u32 = 7;

/// Complete effective configuration and source identity.
#[derive(Debug, Clone)]
pub(crate) struct Config {
    /// Target profile path, resolved against the configuration file's directory.
    pub(crate) profile_path: PathBuf,
    pub(crate) basecalling: BasecallingConfig,
    pub(crate) signal_processing: SignalProcessingConfig,
    pub(crate) callability: CallabilityConfig,
    pub(crate) quality_control: QualityControlConfig,
    pub(crate) alignment: AlignmentConfig,
    pub(crate) sample_reconciliation: SampleReconciliationConfig,
    pub(crate) sanger_evidence: SangerEvidenceConfig,
    pub(crate) variant_calling: VariantCallingConfig,
    pub(crate) source_path: PathBuf,
    pub(crate) source_sha256: String,
}

/// DNA re-calling settings.
#[derive(Debug, Clone)]
pub(crate) struct BasecallingConfig {
    pub(crate) secondary_peak_ratio: f64,
}

/// Observation-only rolling signal-quality settings.
#[derive(Debug, Clone)]
pub(crate) struct SignalProcessingConfig {
    pub(crate) window_size_bases: usize,
    pub(crate) minimum_primary_snr: f64,
    pub(crate) minimum_noisy_windows: usize,
}

/// Signal-derived read-callability settings.
#[derive(Debug, Clone)]
pub(crate) struct CallabilityConfig {
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

/// Relative quality settings.
#[derive(Debug, Clone)]
pub(crate) struct QualityControlConfig {
    pub(crate) penalty_window_size: usize,
    pub(crate) max_relative_quality_score: u8,
}

/// Pairwise alignment settings.
#[derive(Debug, Clone)]
pub(crate) struct AlignmentConfig {
    pub(crate) match_score: i32,
    pub(crate) mismatch_score: i32,
    pub(crate) ambiguous_score: i32,
    pub(crate) gap_open_score: i32,
    pub(crate) gap_extension_score: i32,
    pub(crate) minimum_callable_bases: usize,
    pub(crate) minimum_identity: f64,
}

/// Cross-read overlap admission settings used before sample consensus.
#[derive(Debug, Clone)]
pub(crate) struct SampleReconciliationConfig {
    pub(crate) minimum_comparable_bases: usize,
    pub(crate) minimum_overlap_agreement: f64,
}

/// Sanger support thresholds that the Sanger evidence adapter turns into
/// support vetoes (ADR-0069).
#[derive(Debug, Clone)]
pub(crate) struct SangerEvidenceConfig {
    pub(crate) minimum_peak_height: i32,
    pub(crate) relative_quality_threshold: u8,
}

/// Modality-neutral primary-difference calling settings.
#[derive(Debug, Clone)]
pub(crate) struct VariantCallingConfig {
    pub(crate) max_indel_length: usize,
    /// Calls this close to an uninformative call (beyond the trim interval or
    /// masked as unresolved) cannot support a variant.
    pub(crate) read_end_margin: usize,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct RawConfig {
    schema_version: u32,
    profile: PathBuf,
    basecalling: RawBasecallingConfig,
    signal_processing: RawSignalProcessingConfig,
    callability: RawCallabilityConfig,
    quality_control: RawQualityControlConfig,
    alignment: RawAlignmentConfig,
    sample_reconciliation: RawSampleReconciliationConfig,
    sanger_evidence: RawSangerEvidenceConfig,
    variant_calling: RawVariantCallingConfig,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawBasecallingConfig {
    secondary_peak_ratio: f64,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawSignalProcessingConfig {
    window_size_bases: usize,
    minimum_primary_snr: f64,
    minimum_noisy_windows: usize,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawCallabilityConfig {
    window_calls: usize,
    onset_defect_fraction: f64,
    exit_defect_fraction: f64,
    minimum_main_share: f64,
    maximum_far_share: f64,
    minimum_shadow_share: f64,
    weak_amplitude_fraction: f64,
    repeat_min_length: usize,
    minimum_callable_calls: usize,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawQualityControlConfig {
    penalty_window_size: usize,
    max_relative_quality_score: u8,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawAlignmentConfig {
    match_score: i32,
    mismatch_score: i32,
    ambiguous_score: i32,
    gap_open_score: i32,
    gap_extension_score: i32,
    minimum_callable_bases: usize,
    minimum_identity: f64,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawSampleReconciliationConfig {
    minimum_comparable_bases: usize,
    minimum_overlap_agreement: f64,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawSangerEvidenceConfig {
    minimum_peak_height: i32,
    relative_quality_threshold: u8,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawVariantCallingConfig {
    max_indel_length: usize,
    read_end_margin: usize,
}

impl RawConfig {
    pub(super) fn validate(self, source_path: PathBuf, source_sha256: String) -> Result<Config> {
        if self.schema_version != SCHEMA_VERSION {
            return Err(ConfigError::UnsupportedSchemaVersion {
                found: self.schema_version,
                expected: SCHEMA_VERSION,
            }
            .into());
        }
        require_fraction(
            "basecalling.secondary_peak_ratio",
            self.basecalling.secondary_peak_ratio,
        )?;
        if !(5..=10).contains(&self.signal_processing.window_size_bases) {
            return Err(ConfigError::Constraint(
                "signal_processing.window_size_bases must be in 5..=10",
            )
            .into());
        }
        require_positive_finite(
            "signal_processing.minimum_primary_snr",
            self.signal_processing.minimum_primary_snr,
        )?;
        if self.signal_processing.minimum_noisy_windows < 2 {
            return Err(ConfigError::Constraint(
                "signal_processing.minimum_noisy_windows must be at least 2",
            )
            .into());
        }
        if !(8..=64).contains(&self.callability.window_calls) {
            return Err(
                ConfigError::Constraint("callability.window_calls must be in 8..=64").into(),
            );
        }
        require_fraction(
            "callability.onset_defect_fraction",
            self.callability.onset_defect_fraction,
        )?;
        require_finite_range(
            "callability.exit_defect_fraction",
            self.callability.exit_defect_fraction,
            0.0,
            1.0,
        )?;
        if self.callability.exit_defect_fraction >= self.callability.onset_defect_fraction {
            return Err(ConfigError::Constraint(
                "callability.exit_defect_fraction must be less than callability.onset_defect_fraction",
            )
            .into());
        }
        require_fraction(
            "callability.minimum_main_share",
            self.callability.minimum_main_share,
        )?;
        require_finite_range(
            "callability.maximum_far_share",
            self.callability.maximum_far_share,
            0.0,
            1.0,
        )?;
        require_fraction(
            "callability.minimum_shadow_share",
            self.callability.minimum_shadow_share,
        )?;
        require_finite_range(
            "callability.weak_amplitude_fraction",
            self.callability.weak_amplitude_fraction,
            0.0,
            0.5,
        )?;
        if self.quality_control.penalty_window_size == 0 {
            return Err(ConfigError::Constraint(
                "quality_control.penalty_window_size must be positive",
            )
            .into());
        }
        if self.quality_control.max_relative_quality_score == 0 {
            return Err(ConfigError::Constraint(
                "quality_control.max_relative_quality_score must be positive",
            )
            .into());
        }
        if self.alignment.match_score <= 0 {
            return Err(ConfigError::Constraint("alignment.match_score must be positive").into());
        }
        if self.alignment.mismatch_score >= 0
            || self.alignment.gap_open_score >= 0
            || self.alignment.gap_extension_score >= 0
        {
            return Err(ConfigError::Constraint(
                "alignment mismatch and gap scores must be negative",
            )
            .into());
        }
        if self.alignment.minimum_callable_bases == 0 {
            return Err(ConfigError::Constraint(
                "alignment.minimum_callable_bases must be positive",
            )
            .into());
        }
        require_fraction(
            "alignment.minimum_identity",
            self.alignment.minimum_identity,
        )?;
        if self.sample_reconciliation.minimum_comparable_bases == 0 {
            return Err(ConfigError::Constraint(
                "sample_reconciliation.minimum_comparable_bases must be positive",
            )
            .into());
        }
        require_fraction(
            "sample_reconciliation.minimum_overlap_agreement",
            self.sample_reconciliation.minimum_overlap_agreement,
        )?;
        if self.variant_calling.max_indel_length == 0
            || self.variant_calling.max_indel_length > MAX_INDEL_LENGTH
        {
            return Err(ConfigError::MaxIndelLength {
                maximum: MAX_INDEL_LENGTH,
            }
            .into());
        }
        if self.sanger_evidence.minimum_peak_height <= 0
            || self.sanger_evidence.minimum_peak_height > MAX_PEAK_HEIGHT
        {
            return Err(ConfigError::MinimumPeakHeight {
                maximum: MAX_PEAK_HEIGHT,
            }
            .into());
        }
        if self.sanger_evidence.relative_quality_threshold
            >= self.quality_control.max_relative_quality_score
        {
            return Err(ConfigError::Constraint(
                "sanger_evidence.relative_quality_threshold must be less than quality_control.max_relative_quality_score",
            ).into());
        }
        if self.profile.as_os_str().is_empty() {
            return Err(ConfigError::Constraint("profile must name a target profile file").into());
        }
        if self.callability.repeat_min_length < 2 {
            return Err(ConfigError::Constraint(
                "callability.repeat_min_length must be at least 2",
            )
            .into());
        }
        if self.callability.minimum_callable_calls == 0 {
            return Err(ConfigError::Constraint(
                "callability.minimum_callable_calls must be positive",
            )
            .into());
        }
        let profile_path = source_path
            .parent()
            .unwrap_or_else(|| Path::new(""))
            .join(&self.profile);
        Ok(Config {
            profile_path,
            basecalling: BasecallingConfig {
                secondary_peak_ratio: self.basecalling.secondary_peak_ratio,
            },
            signal_processing: SignalProcessingConfig {
                window_size_bases: self.signal_processing.window_size_bases,
                minimum_primary_snr: self.signal_processing.minimum_primary_snr,
                minimum_noisy_windows: self.signal_processing.minimum_noisy_windows,
            },
            callability: CallabilityConfig {
                window_calls: self.callability.window_calls,
                onset_defect_fraction: self.callability.onset_defect_fraction,
                exit_defect_fraction: self.callability.exit_defect_fraction,
                minimum_main_share: self.callability.minimum_main_share,
                maximum_far_share: self.callability.maximum_far_share,
                minimum_shadow_share: self.callability.minimum_shadow_share,
                weak_amplitude_fraction: self.callability.weak_amplitude_fraction,
                repeat_min_length: self.callability.repeat_min_length,
                minimum_callable_calls: self.callability.minimum_callable_calls,
            },
            quality_control: QualityControlConfig {
                penalty_window_size: self.quality_control.penalty_window_size,
                max_relative_quality_score: self.quality_control.max_relative_quality_score,
            },
            alignment: AlignmentConfig {
                match_score: self.alignment.match_score,
                mismatch_score: self.alignment.mismatch_score,
                ambiguous_score: self.alignment.ambiguous_score,
                gap_open_score: self.alignment.gap_open_score,
                gap_extension_score: self.alignment.gap_extension_score,
                minimum_callable_bases: self.alignment.minimum_callable_bases,
                minimum_identity: self.alignment.minimum_identity,
            },
            sample_reconciliation: SampleReconciliationConfig {
                minimum_comparable_bases: self.sample_reconciliation.minimum_comparable_bases,
                minimum_overlap_agreement: self.sample_reconciliation.minimum_overlap_agreement,
            },
            sanger_evidence: SangerEvidenceConfig {
                minimum_peak_height: self.sanger_evidence.minimum_peak_height,
                relative_quality_threshold: self.sanger_evidence.relative_quality_threshold,
            },
            variant_calling: VariantCallingConfig {
                max_indel_length: self.variant_calling.max_indel_length,
                read_end_margin: self.variant_calling.read_end_margin,
            },
            source_path,
            source_sha256,
        })
    }
}

fn require_fraction(key: &'static str, value: f64) -> Result<()> {
    require_finite_range(key, value, f64::MIN_POSITIVE, 1.0)
}

fn require_positive_finite(key: &'static str, value: f64) -> Result<()> {
    if !value.is_finite() || value <= 0.0 {
        return Err(ConfigError::NotFinitePositive { key }.into());
    }
    Ok(())
}

fn require_finite_range(key: &'static str, value: f64, minimum: f64, maximum: f64) -> Result<()> {
    if !value.is_finite() || value < minimum || value > maximum {
        return Err(ConfigError::NotFiniteInRange {
            key,
            minimum,
            maximum,
        }
        .into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::error::{ConfigError, Error};

    use super::*;

    const VALID: &str = "schema_version=7\nprofile='profiles/target.toml'\n[basecalling]\nsecondary_peak_ratio=0.33\n[signal_processing]\nwindow_size_bases=10\nminimum_primary_snr=3.0\nminimum_noisy_windows=2\n[quality_control]\npenalty_window_size=10\nmax_relative_quality_score=60\n[callability]\nwindow_calls=16\nonset_defect_fraction=0.375\nexit_defect_fraction=0.125\nminimum_main_share=0.35\nmaximum_far_share=0.12\nminimum_shadow_share=0.1\nweak_amplitude_fraction=0.1\nrepeat_min_length=8\nminimum_callable_calls=20\n[alignment]\nmatch_score=3\nmismatch_score=-5\nambiguous_score=0\ngap_open_score=-10\ngap_extension_score=-4\nminimum_callable_bases=20\nminimum_identity=0.8\n[sample_reconciliation]\nminimum_comparable_bases=25\nminimum_overlap_agreement=0.5\n[sanger_evidence]\nminimum_peak_height=150\nrelative_quality_threshold=30\n[variant_calling]\nmax_indel_length=50\nread_end_margin=10\n";

    /// Parses TOML, then returns the typed scientific validation outcome.
    fn validate_raw(text: &str) -> std::result::Result<Result<Config>, toml::de::Error> {
        let raw: RawConfig = toml::from_str(text)?;
        Ok(raw.validate("dna.toml".into(), String::new()))
    }

    #[test]
    fn rejects_out_of_range_values_with_the_offending_key()
    -> std::result::Result<(), toml::de::Error> {
        assert!(matches!(
            validate_raw(&VALID.replace("window_size_bases=10", "window_size_bases=11"))?,
            Err(Error::Config(ConfigError::Constraint(rule)))
                if rule.starts_with("signal_processing.window_size_bases")
        ));
        assert!(matches!(
            validate_raw(&VALID.replace("secondary_peak_ratio=0.33", "secondary_peak_ratio=1.5"))?,
            Err(Error::Config(ConfigError::NotFiniteInRange {
                key: "basecalling.secondary_peak_ratio",
                ..
            }))
        ));
        Ok(())
    }

    fn validate(text: &str) -> std::result::Result<Config, Box<dyn std::error::Error>> {
        let raw: RawConfig = toml::from_str(text)?;
        Ok(raw.validate(PathBuf::from("dna.toml"), String::new())?)
    }

    #[test]
    fn accepts_signal_settings_and_variant_filter_list_of_lists()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        let config = validate(VALID)?;

        assert_eq!(config.signal_processing.window_size_bases, 10);
        assert_eq!(config.signal_processing.minimum_primary_snr, 3.0);
        assert_eq!(config.signal_processing.minimum_noisy_windows, 2);
        assert_eq!(config.sample_reconciliation.minimum_comparable_bases, 25);
        assert_eq!(config.sample_reconciliation.minimum_overlap_agreement, 0.5);
        assert_eq!(config.sanger_evidence.minimum_peak_height, 150);
        assert_eq!(config.sanger_evidence.relative_quality_threshold, 30);
        Ok(())
    }

    #[test]
    fn resolves_the_profile_against_the_configuration_directory()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        let raw: RawConfig = toml::from_str(VALID)?;
        let config = raw.validate(PathBuf::from("deploy/config/dna.toml"), String::new())?;
        assert_eq!(
            config.profile_path,
            PathBuf::from("deploy/config/profiles/target.toml")
        );
        assert_eq!(
            validate(&VALID.replace("'profiles/target.toml'", "'/opt/target.toml'"))?.profile_path,
            PathBuf::from("/opt/target.toml")
        );
        assert!(matches!(
            validate_raw(&VALID.replace("'profiles/target.toml'", "''"))?,
            Err(Error::Config(ConfigError::Constraint(rule))) if rule.starts_with("profile")
        ));
        assert!(
            toml::from_str::<RawConfig>(&VALID.replace("profile='profiles/target.toml'\n", ""))
                .is_err()
        );
        Ok(())
    }

    #[test]
    fn rejects_old_schema_and_missing_required_fields() -> std::result::Result<(), toml::de::Error>
    {
        assert!(matches!(
            validate_raw(&VALID.replace("schema_version=7", "schema_version=6"))?,
            Err(Error::Config(ConfigError::UnsupportedSchemaVersion {
                found: 6,
                expected: 7
            }))
        ));
        assert!(
            toml::from_str::<RawConfig>(&VALID.replace(
                "max_indel_length=50\n",
                "max_indel_length=50\nminimum_peak_height=150\n"
            ))
            .is_err()
        );
        assert!(
            toml::from_str::<RawConfig>(&VALID.replace("minimum_primary_snr=3.0\n", "")).is_err()
        );
        assert!(
            toml::from_str::<RawConfig>(&VALID.replace("minimum_peak_height=150\n", "")).is_err()
        );
        assert!(
            toml::from_str::<RawConfig>(&VALID.replace("minimum_comparable_bases=25\n", ""))
                .is_err()
        );
        Ok(())
    }

    #[test]
    fn rejects_invalid_signal_settings() {
        for invalid in [
            VALID.replace("window_size_bases=10", "window_size_bases=4"),
            VALID.replace("window_size_bases=10", "window_size_bases=11"),
            VALID.replace("minimum_noisy_windows=2", "minimum_noisy_windows=1"),
            VALID.replace("minimum_primary_snr=3.0", "minimum_primary_snr=0.0"),
            VALID.replace("minimum_primary_snr=3.0", "minimum_primary_snr=nan"),
        ] {
            assert!(
                validate(&invalid).is_err(),
                "unexpectedly accepted {invalid}"
            );
        }
    }

    #[test]
    fn rejects_invalid_callability_settings() {
        for invalid in [
            VALID.replace("window_calls=16", "window_calls=7"),
            VALID.replace("exit_defect_fraction=0.125", "exit_defect_fraction=0.5"),
            VALID.replace("minimum_main_share=0.35", "minimum_main_share=0.0"),
            VALID.replace("minimum_main_share=0.35", "minimum_main_share=1.5"),
            VALID.replace("maximum_far_share=0.12", "maximum_far_share=-0.1"),
            VALID.replace("maximum_far_share=0.12", "maximum_far_share=nan"),
            VALID.replace("minimum_shadow_share=0.1", "minimum_shadow_share=0.0"),
            VALID.replace("weak_amplitude_fraction=0.1", "weak_amplitude_fraction=0.6"),
        ] {
            assert!(
                validate(&invalid).is_err(),
                "unexpectedly accepted {invalid}"
            );
        }
        assert!(
            toml::from_str::<RawConfig>(&VALID.replace(
                "minimum_shadow_share=0.1\n",
                "minimum_shadow_share=0.1\nshift_coherence=0.75\n"
            ))
            .is_err()
        );
    }

    #[test]
    fn rejects_invalid_sample_reconciliation_settings() {
        for invalid in [
            VALID.replace("minimum_comparable_bases=25", "minimum_comparable_bases=0"),
            VALID.replace(
                "minimum_overlap_agreement=0.5",
                "minimum_overlap_agreement=0.0",
            ),
            VALID.replace(
                "minimum_overlap_agreement=0.5",
                "minimum_overlap_agreement=1.1",
            ),
            VALID.replace(
                "minimum_overlap_agreement=0.5",
                "minimum_overlap_agreement=nan",
            ),
        ] {
            assert!(
                validate(&invalid).is_err(),
                "unexpectedly accepted {invalid}"
            );
        }
    }

    #[test]
    fn rejects_invalid_filter_thresholds() {
        for invalid in [
            VALID.replace("minimum_peak_height=150", "minimum_peak_height=0"),
            VALID.replace("minimum_peak_height=150", "minimum_peak_height=32768"),
            VALID.replace(
                "relative_quality_threshold=30",
                "relative_quality_threshold=60",
            ),
        ] {
            assert!(
                validate(&invalid).is_err(),
                "unexpectedly accepted {invalid}"
            );
        }
        assert!(
            toml::from_str::<RawConfig>(
                &VALID.replace("read_end_margin=10", "read_end_margin=10\nregions=[[1,10]]")
            )
            .is_err(),
            "regions belong to the target profile"
        );
    }
}
