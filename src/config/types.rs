//! The validated configuration envelope: target profile, source identity, and
//! the sections each plugin owns (ADR-0069).

use std::path::{Path, PathBuf};

use serde::Deserialize;

use dna_core::alignment::RawAlignmentConfig;
use dna_core::read_call::{CoreConfig, RawCoreConfig};
use dna_core::sample::RawSampleReconciliationConfig;
use dna_core::variant_calling::RawVariantCallingConfig;
use dna_kernel::error::{ConfigError, Result};
use dna_sanger::basecalling::RawBasecallingConfig;
use dna_sanger::callability::RawCallabilityConfig;
use dna_sanger::quality_control::RawQualityControlConfig;
use dna_sanger::read_processing::{RawSangerConfig, RawSangerEvidenceConfig, SangerConfig};
use dna_sanger::signal_processing::RawSignalProcessingConfig;

/// Configuration schema version this build accepts.
const SCHEMA_VERSION: u32 = 7;

/// Complete effective configuration and source identity.
#[derive(Debug, Clone)]
pub(crate) struct Config {
    /// Target profile path, resolved against the configuration file's directory.
    pub(crate) profile_path: PathBuf,
    /// Sections owned by the Sanger modality plugin.
    pub(crate) sanger: SangerConfig,
    /// Sections owned by the core plugin.
    pub(crate) core: CoreConfig,
    pub(crate) source_path: PathBuf,
    pub(crate) source_sha256: String,
}

/// The configuration as written: one table per plugin-owned section.
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

impl RawConfig {
    pub(super) fn validate(self, source_path: PathBuf, source_sha256: String) -> Result<Config> {
        if self.schema_version != SCHEMA_VERSION {
            return Err(ConfigError::UnsupportedSchemaVersion {
                found: self.schema_version,
                expected: SCHEMA_VERSION,
            }
            .into());
        }
        if self.profile.as_os_str().is_empty() {
            return Err(ConfigError::Constraint("profile must name a target profile file").into());
        }
        let sanger = RawSangerConfig {
            basecalling: self.basecalling,
            signal_processing: self.signal_processing,
            callability: self.callability,
            quality_control: self.quality_control,
            sanger_evidence: self.sanger_evidence,
        }
        .validate()?;
        let core = RawCoreConfig {
            alignment: self.alignment,
            variant_calling: self.variant_calling,
            sample_reconciliation: self.sample_reconciliation,
        }
        .validate()?;
        let profile_path = source_path
            .parent()
            .unwrap_or_else(|| Path::new(""))
            .join(&self.profile);
        Ok(Config {
            profile_path,
            sanger,
            core,
            source_path,
            source_sha256,
        })
    }
}

#[cfg(test)]
mod tests {
    use dna_kernel::error::{ConfigError, Error};

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

        assert_eq!(config.sanger.signal_processing.window_size_bases, 10);
        assert_eq!(config.sanger.signal_processing.minimum_primary_snr, 3.0);
        assert_eq!(config.sanger.signal_processing.minimum_noisy_windows, 2);
        assert_eq!(
            config.core.sample_reconciliation.minimum_comparable_bases,
            25
        );
        assert_eq!(
            config.core.sample_reconciliation.minimum_overlap_agreement,
            0.5
        );
        assert_eq!(config.sanger.sanger_evidence.minimum_peak_height, 150);
        assert_eq!(config.sanger.sanger_evidence.relative_quality_threshold, 30);
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
