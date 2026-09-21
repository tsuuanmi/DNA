//! Typed and validated configuration records.

use std::path::PathBuf;

use serde::Deserialize;

use crate::config::defaults::{MAX_INDEL_LENGTH, MAX_PEAK_HEIGHT, MAX_REFERENCE_LENGTH};
use crate::error::{Error, Result};
use crate::model::reference::ReferenceTopology;

/// Complete effective configuration and source identity.
#[derive(Debug, Clone)]
pub struct Config {
    pub(crate) reference: ReferenceConfig,
    pub(crate) basecalling: BasecallingConfig,
    pub(crate) signal_processing: DNAProcessingConfig,
    pub(crate) quality_control: QualityControlConfig,
    pub(crate) alignment: AlignmentConfig,
    pub(crate) sample_reconciliation: SampleReconciliationConfig,
    pub(crate) variant_calling: VariantCallingConfig,
    pub(crate) source_path: PathBuf,
    pub(crate) source_sha256: String,
}

/// Reference interpretation settings.
#[derive(Debug, Clone)]
pub struct ReferenceConfig {
    pub(crate) topology: ReferenceTopology,
}

/// DNA re-calling settings.
#[derive(Debug, Clone)]
pub struct BasecallingConfig {
    pub(crate) secondary_peak_ratio: f64,
}

/// Observation-only rolling signal-quality settings.
#[derive(Debug, Clone)]
pub struct DNAProcessingConfig {
    pub(crate) window_size_bases: usize,
    pub(crate) minimum_primary_snr: f64,
    pub(crate) minimum_noisy_windows: usize,
}

/// Relative score and trimming settings.
#[derive(Debug, Clone)]
pub struct QualityControlConfig {
    pub(crate) trim_window_size: usize,
    pub(crate) best_section_fraction: f64,
    pub(crate) max_relative_quality_score: u8,
    pub(crate) trim_stringency: f64,
    pub(crate) minimum_retained_bases: usize,
}

/// Pairwise alignment settings.
#[derive(Debug, Clone)]
pub struct AlignmentConfig {
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
pub struct SampleReconciliationConfig {
    pub(crate) minimum_comparable_bases: usize,
    pub(crate) minimum_overlap_agreement: f64,
}

/// Primary-difference calling settings.
#[derive(Debug, Clone)]
pub struct VariantCallingConfig {
    pub(crate) max_indel_length: usize,
    pub(crate) minimum_peak_height: i32,
    pub(crate) relative_quality_threshold: u8,
    pub(crate) regions: Vec<[usize; 2]>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct RawConfig {
    schema_version: u32,
    reference: RawReferenceConfig,
    basecalling: RawBasecallingConfig,
    signal_processing: RawDNAProcessingConfig,
    quality_control: RawQualityControlConfig,
    alignment: RawAlignmentConfig,
    sample_reconciliation: RawSampleReconciliationConfig,
    variant_calling: RawVariantCallingConfig,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawReferenceConfig {
    topology: ReferenceTopology,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawBasecallingConfig {
    secondary_peak_ratio: f64,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawDNAProcessingConfig {
    window_size_bases: usize,
    minimum_primary_snr: f64,
    minimum_noisy_windows: usize,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawQualityControlConfig {
    trim_window_size: usize,
    best_section_fraction: f64,
    max_relative_quality_score: u8,
    trim_stringency: f64,
    minimum_retained_bases: usize,
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
struct RawVariantCallingConfig {
    max_indel_length: usize,
    minimum_peak_height: i32,
    relative_quality_threshold: u8,
    regions: Vec<[usize; 2]>,
}

impl RawConfig {
    pub(super) fn validate(self, source_path: PathBuf, source_sha256: String) -> Result<Config> {
        if self.schema_version != 5 {
            return Err(Error::Config(format!(
                "unsupported schema_version {}; expected 5",
                self.schema_version
            )));
        }
        require_fraction(
            "basecalling.secondary_peak_ratio",
            self.basecalling.secondary_peak_ratio,
        )?;
        if !(5..=10).contains(&self.signal_processing.window_size_bases) {
            return Err(Error::Config(
                "signal_processing.window_size_bases must be in 5..=10".into(),
            ));
        }
        require_positive_finite(
            "signal_processing.minimum_primary_snr",
            self.signal_processing.minimum_primary_snr,
        )?;
        if self.signal_processing.minimum_noisy_windows < 2 {
            return Err(Error::Config(
                "signal_processing.minimum_noisy_windows must be at least 2".into(),
            ));
        }
        if self.quality_control.trim_window_size == 0 {
            return Err(Error::Config(
                "quality_control.trim_window_size must be positive".into(),
            ));
        }
        require_fraction(
            "quality_control.best_section_fraction",
            self.quality_control.best_section_fraction,
        )?;
        if self.quality_control.max_relative_quality_score == 0 {
            return Err(Error::Config(
                "quality_control.max_relative_quality_score must be positive".into(),
            ));
        }
        require_finite_range(
            "quality_control.trim_stringency",
            self.quality_control.trim_stringency,
            0.0,
            9.0,
        )?;
        if self.quality_control.minimum_retained_bases == 0 {
            return Err(Error::Config(
                "quality_control.minimum_retained_bases must be positive".into(),
            ));
        }
        if self.alignment.match_score <= 0 {
            return Err(Error::Config(
                "alignment.match_score must be positive".into(),
            ));
        }
        if self.alignment.mismatch_score >= 0
            || self.alignment.gap_open_score >= 0
            || self.alignment.gap_extension_score >= 0
        {
            return Err(Error::Config(
                "alignment mismatch and gap scores must be negative".into(),
            ));
        }
        if self.alignment.minimum_callable_bases == 0 {
            return Err(Error::Config(
                "alignment.minimum_callable_bases must be positive".into(),
            ));
        }
        require_fraction(
            "alignment.minimum_identity",
            self.alignment.minimum_identity,
        )?;
        if self.sample_reconciliation.minimum_comparable_bases == 0 {
            return Err(Error::Config(
                "sample_reconciliation.minimum_comparable_bases must be positive".into(),
            ));
        }
        require_fraction(
            "sample_reconciliation.minimum_overlap_agreement",
            self.sample_reconciliation.minimum_overlap_agreement,
        )?;
        if self.variant_calling.max_indel_length == 0
            || self.variant_calling.max_indel_length > MAX_INDEL_LENGTH
        {
            return Err(Error::Config(format!(
                "variant_calling.max_indel_length must be in 1..={MAX_INDEL_LENGTH}"
            )));
        }
        if self.variant_calling.minimum_peak_height <= 0
            || self.variant_calling.minimum_peak_height > MAX_PEAK_HEIGHT
        {
            return Err(Error::Config(format!(
                "variant_calling.minimum_peak_height must be in 1..={MAX_PEAK_HEIGHT}"
            )));
        }
        if self.variant_calling.relative_quality_threshold
            >= self.quality_control.max_relative_quality_score
        {
            return Err(Error::Config(
                "variant_calling.relative_quality_threshold must be less than quality_control.max_relative_quality_score"
                    .into(),
            ));
        }
        if self.variant_calling.regions.is_empty() {
            return Err(Error::Config(
                "variant_calling.regions must contain at least one inclusive range".into(),
            ));
        }
        for (index, region) in self.variant_calling.regions.iter().enumerate() {
            let [start, end] = *region;
            if start == 0 || start > end || end > MAX_REFERENCE_LENGTH {
                return Err(Error::Config(format!(
                    "variant_calling.regions[{index}] must satisfy 1 <= start <= end <= {MAX_REFERENCE_LENGTH}"
                )));
            }
        }
        Ok(Config {
            reference: ReferenceConfig {
                topology: self.reference.topology,
            },
            basecalling: BasecallingConfig {
                secondary_peak_ratio: self.basecalling.secondary_peak_ratio,
            },
            signal_processing: DNAProcessingConfig {
                window_size_bases: self.signal_processing.window_size_bases,
                minimum_primary_snr: self.signal_processing.minimum_primary_snr,
                minimum_noisy_windows: self.signal_processing.minimum_noisy_windows,
            },
            quality_control: QualityControlConfig {
                trim_window_size: self.quality_control.trim_window_size,
                best_section_fraction: self.quality_control.best_section_fraction,
                max_relative_quality_score: self.quality_control.max_relative_quality_score,
                trim_stringency: self.quality_control.trim_stringency,
                minimum_retained_bases: self.quality_control.minimum_retained_bases,
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
            variant_calling: VariantCallingConfig {
                max_indel_length: self.variant_calling.max_indel_length,
                minimum_peak_height: self.variant_calling.minimum_peak_height,
                relative_quality_threshold: self.variant_calling.relative_quality_threshold,
                regions: self.variant_calling.regions,
            },
            source_path,
            source_sha256,
        })
    }
}

fn require_fraction(name: &str, value: f64) -> Result<()> {
    require_finite_range(name, value, f64::MIN_POSITIVE, 1.0)
}

fn require_positive_finite(name: &str, value: f64) -> Result<()> {
    if !value.is_finite() || value <= 0.0 {
        return Err(Error::Config(format!("{name} must be finite and positive")));
    }
    Ok(())
}

fn require_finite_range(name: &str, value: f64, minimum: f64, maximum: f64) -> Result<()> {
    if !value.is_finite() || value < minimum || value > maximum {
        return Err(Error::Config(format!(
            "{name} must be finite and in [{minimum}, {maximum}]"
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    const VALID: &str = "schema_version=5\n[reference]\ntopology='circular'\n[basecalling]\nsecondary_peak_ratio=0.33\n[signal_processing]\nwindow_size_bases=10\nminimum_primary_snr=3.0\nminimum_noisy_windows=2\n[quality_control]\ntrim_window_size=10\nbest_section_fraction=0.1\nmax_relative_quality_score=60\ntrim_stringency=7.0\nminimum_retained_bases=20\n[alignment]\nmatch_score=3\nmismatch_score=-5\nambiguous_score=0\ngap_open_score=-10\ngap_extension_score=-4\nminimum_callable_bases=20\nminimum_identity=0.8\n[sample_reconciliation]\nminimum_comparable_bases=25\nminimum_overlap_agreement=0.5\n[variant_calling]\nmax_indel_length=50\nminimum_peak_height=150\nrelative_quality_threshold=30\nregions=[[16024,16365],[73,340],[438,576]]\n";

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
        assert_eq!(config.variant_calling.minimum_peak_height, 150);
        assert_eq!(config.variant_calling.relative_quality_threshold, 30);
        assert_eq!(
            config.variant_calling.regions,
            vec![[16024, 16365], [73, 340], [438, 576]]
        );
        Ok(())
    }

    #[test]
    fn rejects_old_schema_and_missing_required_fields() {
        assert!(validate(&VALID.replace("schema_version=5", "schema_version=4")).is_err());
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
    fn rejects_invalid_filter_thresholds_and_regions() {
        for invalid in [
            VALID.replace("minimum_peak_height=150", "minimum_peak_height=0"),
            VALID.replace("minimum_peak_height=150", "minimum_peak_height=32768"),
            VALID.replace(
                "relative_quality_threshold=30",
                "relative_quality_threshold=60",
            ),
            VALID.replace("regions=[[16024,16365],[73,340],[438,576]]", "regions=[]"),
            VALID.replace(
                "regions=[[16024,16365],[73,340],[438,576]]",
                "regions=[[0,1]]",
            ),
            VALID.replace(
                "regions=[[16024,16365],[73,340],[438,576]]",
                "regions=[[2,1]]",
            ),
            VALID.replace(
                "regions=[[16024,16365],[73,340],[438,576]]",
                "regions=[[1,50001]]",
            ),
        ] {
            assert!(
                validate(&invalid).is_err(),
                "unexpectedly accepted {invalid}"
            );
        }
        assert!(
            toml::from_str::<RawConfig>(&VALID.replace(
                "regions=[[16024,16365],[73,340],[438,576]]",
                "regions=[[1]]"
            ))
            .is_err()
        );
    }
}
