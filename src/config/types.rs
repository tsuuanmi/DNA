//! The validated configuration envelope: target profile, source identity, and
//! the sections each plugin owns (ADR-0069).

use std::path::{Path, PathBuf};

use serde::Deserialize;

use crate::error::{ConfigError, Result};
use dna_align::RawAlignmentConfig;
use dna_call::RawVariantCallingConfig;
use dna_call::read_call::{CallConfig, RawCallConfig};
use dna_consensus::{ConsensusConfig, RawConsensusConfig};
use dna_sample::{RawSampleReconciliationConfig, SampleReconciliationConfig};
use dna_sanger::basecalling::RawBasecallingConfig;
use dna_sanger::callability::RawCallabilityConfig;
use dna_sanger::quality_control::RawQualityControlConfig;
use dna_sanger::read_processing::{RawSangerConfig, RawSangerEvidenceConfig, SangerConfig};
use dna_sanger::signal_processing::RawSignalProcessingConfig;
use dna_sanger::tract::RawTraceConfig;

/// Configuration schema version this build accepts.
const SCHEMA_VERSION: u32 = 20;

/// Complete effective configuration and source identity.
#[derive(Debug, Clone)]
pub(crate) struct Config {
    /// Target profile path, resolved against the configuration file's directory.
    pub(crate) profile_path: PathBuf,
    /// Sections owned by the Sanger modality plugin.
    pub(crate) sanger: SangerConfig,
    /// The core plugin's sections a called read needs.
    pub(crate) core: CallConfig,
    /// The consensus plugin's `[consensus]` section.
    pub(crate) consensus: ConsensusConfig,
    /// The core plugin's `[sample_reconciliation]` section.
    pub(crate) sample: SampleReconciliationConfig,
    /// The review plugin's `[review]` section.
    pub(crate) review: ReviewConfig,
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
    trace: RawTraceConfig,
    variant_calling: RawVariantCallingConfig,
    consensus: RawConsensusConfig,
    review: RawReviewConfig,
}

/// The review plugin's thresholds for a decision that one strand or an
/// anchored end alone cannot assure (SRS-API-014).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ReviewConfig {
    /// Fewest bases a decided run may exceed its reference run by before a
    /// length that rests on an anchored end or a single read is a check.
    pub(crate) run_growth: usize,
    /// Positions from a reportable region's end that the call does not
    /// analyse within which any difference is a check.
    pub(crate) edge_zone: usize,
    /// Positions from a reportable region's first or last position within
    /// which a deletion that rests on a single read is a check.
    pub(crate) deletion_zone: usize,
    /// Equal bases in a row after which a read's in-phase end of an item no
    /// longer shows bases the review compares with the decision.
    pub(crate) phase_run: usize,
}

/// The `[review]` section as written.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct RawReviewConfig {
    run_growth: usize,
    edge_zone: usize,
    deletion_zone: usize,
    phase_run: usize,
}

impl RawReviewConfig {
    fn validate(self) -> Result<ReviewConfig> {
        if self.run_growth == 0 {
            return Err(ConfigError::Constraint("review.run_growth must be positive").into());
        }
        if self.phase_run < 2 {
            return Err(ConfigError::Constraint("review.phase_run must be at least 2").into());
        }
        Ok(ReviewConfig {
            run_growth: self.run_growth,
            edge_zone: self.edge_zone,
            deletion_zone: self.deletion_zone,
            phase_run: self.phase_run,
        })
    }
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
            trace: self.trace,
        }
        .validate()?;
        let core = RawCallConfig {
            alignment: self.alignment,
            variant_calling: self.variant_calling,
        }
        .validate()?;
        let sample = self.sample_reconciliation.validate()?;
        let consensus = self.consensus.validate()?;
        let review = self.review.validate()?;
        let profile_path = source_path
            .parent()
            .unwrap_or_else(|| Path::new(""))
            .join(&self.profile);
        Ok(Config {
            profile_path,
            sanger,
            core,
            consensus,
            sample,
            review,
            source_path,
            source_sha256,
        })
    }
}

#[cfg(test)]
mod tests {
    use crate::error::{ConfigError, Error, KernelError};

    use super::*;

    const VALID: &str = "schema_version=20\nprofile='profiles/target.toml'\n[basecalling]\nsecondary_peak_ratio=0.33\n[signal_processing]\nwindow_size_bases=10\nminimum_primary_snr=3.0\nminimum_noisy_windows=2\nminimum_noise_sigma=1.0\nsecond_floor=0.21\nrun_signal_flank=2\nrun_signal_resolution=1000\nrun_signal_minimum_reads=40\nrun_signal_minimum_calls=40\nrun_signal_margin=0.05\nrun_signal_mad_factor=4.0\nsecond_cause_noise_ratio=3.0\nsecond_cause_tail_offset=0.4\nsecond_cause_blob_width=1.8\nsecond_cause_width_calls=6\ndecompose_unexplained_quantile=0.99\n[quality_control]\npenalty_window_size=10\nmax_relative_quality_score=60\ntrim_context_margin=10\n[callability]\nwindow_calls=16\nonset_defect_fraction=0.375\nexit_defect_fraction=0.125\nminimum_main_share=0.35\nmaximum_far_share=0.12\nminimum_shadow_share=0.1\nweak_amplitude_fraction=0.1\nrepeat_min_length=8\nminimum_callable_calls=20\nrepeat_window_calls=16\ndouble_peak_ratio=0.33\nspacing_defect=0.5\nspacing_half_window=8\n[alignment]\nmatch_score=3\nmismatch_score=-5\nambiguous_score=0\ngap_open_score=-10\ngap_extension_score=-4\nminimum_callable_bases=20\nminimum_identity=0.8\nmax_stretch=64\n[sample_reconciliation]\nminimum_comparable_bases=25\nminimum_overlap_agreement=0.5\n[sanger_evidence]\nminimum_peak_height=150\nrelative_quality_threshold=30\npeak_rescue_quality=50\npeak_rescue_second_share=0.20\nsignal_edge_margin=10\nleading_height_fraction=1.0\npost_repeat_margin=10\npost_repeat_minimum_calls=2\nedge_identity_margin=10\nedge_height_fraction=1.0\nedge_second_share=0.10\nedge_width_factor=1.5\nedge_spacing_minimum=0.5\nedge_spacing_maximum=1.5\nstart_identity_margin=10\nstart_norm_calls=30\nstart_height_fraction=0.5\nstart_second_share=0.10\n[trace]\nlocate_min_run=4\nlocate_min_before=2\nspan_below=10\nspan_above=3\nmin_length=3\ndownstream_calls=60\njoint_lengths=6\niterations=6\nsignature_calls=200\nsignature_gap=3\nleak_min_run=3\nleak_pad=2\nleak_min_calls=10\nnear_equal_z=2.0\n[variant_calling]\nmax_indel_length=50\nread_end_margin=10\n[consensus]\nseeded_run=4\nmaximum_undecided_gap=20\nminimum_resolved_bases=20\nsite_mismatch_cost=2\nsite_gap_open_cost=3\nsite_gap_extension_cost=2\ndissent_edge=40\n[review]\nrun_growth=3\nedge_zone=20\ndeletion_zone=20\nphase_run=4\n";

    /// Parses TOML, then returns the typed scientific validation outcome.
    fn validate_raw(text: &str) -> std::result::Result<Result<Config>, toml::de::Error> {
        let raw: RawConfig = toml::from_str(text)?;
        Ok(raw.validate("dna.toml".into(), String::new()))
    }

    // SRS-CFG-003
    #[test]
    fn rejects_out_of_range_values_with_the_offending_key()
    -> std::result::Result<(), toml::de::Error> {
        assert!(matches!(
            validate_raw(&VALID.replace("window_size_bases=10", "window_size_bases=11"))?,
            Err(Error::Kernel(KernelError::Config(ConfigError::Constraint(rule))))
                if rule.starts_with("signal_processing.window_size_bases")
        ));
        assert!(matches!(
            validate_raw(&VALID.replace("secondary_peak_ratio=0.33", "secondary_peak_ratio=1.5"))?,
            Err(Error::Kernel(KernelError::Config(
                ConfigError::NotFiniteInRange {
                    key: "basecalling.secondary_peak_ratio",
                    ..
                }
            )))
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
        assert_eq!(config.sample.minimum_comparable_bases, 25);
        assert_eq!(config.sample.minimum_overlap_agreement, 0.5);
        assert_eq!(config.sanger.sanger_evidence.minimum_peak_height, 150);
        assert_eq!(config.sanger.sanger_evidence.relative_quality_threshold, 30);
        Ok(())
    }

    // SRS-CFG-002
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
            Err(Error::Kernel(KernelError::Config(ConfigError::Constraint(rule)))) if rule.starts_with("profile")
        ));
        assert!(
            toml::from_str::<RawConfig>(&VALID.replace("profile='profiles/target.toml'\n", ""))
                .is_err()
        );
        Ok(())
    }

    // SRS-CFG-002 SRS-CFG-003
    #[test]
    fn rejects_old_schema_and_missing_required_fields() -> std::result::Result<(), toml::de::Error>
    {
        assert!(matches!(
            validate_raw(&VALID.replace("schema_version=20", "schema_version=10"))?,
            Err(Error::Kernel(KernelError::Config(
                ConfigError::UnsupportedSchemaVersion {
                    found: 10,
                    expected: 20
                }
            )))
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
        assert!(toml::from_str::<RawConfig>(&VALID.replace("deletion_zone=20\n", "")).is_err());
        Ok(())
    }

    // SRS-SIG-002 SRS-CFG-003 SRS-CFG-009
    #[test]
    fn rejects_invalid_signal_settings() {
        for invalid in [
            VALID.replace("window_size_bases=10", "window_size_bases=4"),
            VALID.replace("window_size_bases=10", "window_size_bases=11"),
            VALID.replace("minimum_noisy_windows=2", "minimum_noisy_windows=1"),
            VALID.replace("minimum_primary_snr=3.0", "minimum_primary_snr=0.0"),
            VALID.replace("minimum_primary_snr=3.0", "minimum_primary_snr=nan"),
            VALID.replace("minimum_noise_sigma=1.0", "minimum_noise_sigma=0.0"),
        ] {
            assert!(
                validate(&invalid).is_err(),
                "unexpectedly accepted {invalid}"
            );
        }
    }

    // SRS-CALL-002 SRS-CFG-003 SRS-CFG-008
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
            VALID.replace("repeat_window_calls=16", "repeat_window_calls=0"),
            VALID.replace("double_peak_ratio=0.33", "double_peak_ratio=1.5"),
            VALID.replace("spacing_defect=0.5", "spacing_defect=0.0"),
            VALID.replace("spacing_half_window=8", "spacing_half_window=0"),
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

    // SRS-CFG-003 SRS-CFG-009
    #[test]
    fn rejects_invalid_read_edge_settings() {
        for invalid in [
            VALID.replace("leading_height_fraction=1.0", "leading_height_fraction=0.0"),
            VALID.replace("post_repeat_minimum_calls=2", "post_repeat_minimum_calls=0"),
            VALID.replace("edge_height_fraction=1.0", "edge_height_fraction=nan"),
            VALID.replace("edge_second_share=0.10", "edge_second_share=1.5"),
            VALID.replace("edge_width_factor=1.5", "edge_width_factor=-1.0"),
            VALID.replace("edge_spacing_minimum=0.5", "edge_spacing_minimum=2.0"),
            VALID.replace("start_norm_calls=30", "start_norm_calls=0"),
            VALID.replace("start_height_fraction=0.5", "start_height_fraction=0.0"),
        ] {
            assert!(
                validate(&invalid).is_err(),
                "unexpectedly accepted {invalid}"
            );
        }
    }

    // SRS-CFG-003 SRS-CFG-009
    #[test]
    fn rejects_invalid_consensus_and_stretch_settings() {
        for invalid in [
            VALID.replace("max_stretch=64", "max_stretch=0"),
            VALID.replace("seeded_run=4", "seeded_run=1"),
            VALID.replace("minimum_resolved_bases=20", "minimum_resolved_bases=0"),
            VALID.replace("site_mismatch_cost=2", "site_mismatch_cost=6"),
            VALID.replace("site_gap_extension_cost=2", "site_gap_extension_cost=0"),
        ] {
            assert!(
                validate(&invalid).is_err(),
                "unexpectedly accepted {invalid}"
            );
        }
    }

    // SRS-CFG-003 SRS-CFG-007
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

    // SRS-CFG-003 SRS-CFG-006 SRS-PRF-001
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
