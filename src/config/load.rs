//! Configuration path resolution, parsing, and identity.

use std::env;
use std::fs;
use std::path::{Path, PathBuf};

use crate::config::defaults::{DEFAULT_CONFIG_PATH, MAX_CONFIG_BYTES};
use crate::config::types::{Config, RawConfig};
use crate::error::{ConfigError, KernelError, Result};
use dna_kernel::checksum::hex_sha256;

/// Resolves the authoritative configuration path without performing I/O.
pub(crate) fn resolve_path() -> PathBuf {
    env::var_os("DNA_CONFIG").map_or_else(|| PathBuf::from(DEFAULT_CONFIG_PATH), PathBuf::from)
}

/// Loads and validates one explicit configuration file.
pub(crate) fn load_path(path: &Path) -> Result<Config> {
    let metadata = fs::metadata(path).map_err(|source| KernelError::Read {
        kind: "configuration",
        path: path.to_path_buf(),
        source,
    })?;
    if metadata.len() > MAX_CONFIG_BYTES as u64 {
        return Err(ConfigError::TooLarge {
            maximum: MAX_CONFIG_BYTES,
        }
        .into());
    }
    let bytes = fs::read(path).map_err(|source| KernelError::Read {
        kind: "configuration",
        path: path.to_path_buf(),
        source,
    })?;
    if bytes.len() > MAX_CONFIG_BYTES {
        return Err(ConfigError::TooLarge {
            maximum: MAX_CONFIG_BYTES,
        }
        .into());
    }
    let text = std::str::from_utf8(&bytes).map_err(ConfigError::NotUtf8)?;
    let raw: RawConfig = toml::from_str(text).map_err(|source| KernelError::ConfigParse {
        path: path.to_path_buf(),
        source: Box::new(source),
    })?;
    raw.validate(path.to_path_buf(), hex_sha256(&bytes))
}

#[cfg(test)]
mod tests {
    use std::fs;

    use tempfile::tempdir;

    use super::*;

    #[test]
    fn rejects_unknown_keys() -> Result<()> {
        let directory = tempdir().map_err(|source| KernelError::Output {
            path: PathBuf::from("temporary directory"),
            source,
        })?;
        let path = directory.path().join("dna.toml");
        fs::write(
            &path,
            "schema_version=13\nunknown=true\nprofile='profiles/target.toml'\n[basecalling]\nsecondary_peak_ratio=0.33\n[signal_processing]\nwindow_size_bases=10\nminimum_primary_snr=3.0\nminimum_noisy_windows=2\nminimum_noise_sigma=1.0\nsecond_floor=0.21\nrun_signal_minimum_reads=40\nrun_signal_minimum_calls=40\nrun_signal_margin=0.05\nrun_signal_mad_factor=4.0\n[quality_control]\npenalty_window_size=10\nmax_relative_quality_score=60\ntrim_context_margin=10\n[callability]\nwindow_calls=16\nonset_defect_fraction=0.375\nexit_defect_fraction=0.125\nminimum_main_share=0.35\nmaximum_far_share=0.12\nminimum_shadow_share=0.1\nweak_amplitude_fraction=0.1\nrepeat_min_length=8\nminimum_callable_calls=20\nrepeat_window_calls=16\ndouble_peak_ratio=0.33\nspacing_defect=0.5\nspacing_half_window=8\n[alignment]\nmatch_score=3\nmismatch_score=-5\nambiguous_score=0\ngap_open_score=-10\ngap_extension_score=-4\nminimum_callable_bases=20\nminimum_identity=0.8\nmax_stretch=64\n[sample_reconciliation]\nminimum_comparable_bases=25\nminimum_overlap_agreement=0.5\n[sanger_evidence]\nminimum_peak_height=150\nrelative_quality_threshold=30\npeak_rescue_quality=50\nsignal_edge_margin=10\nleading_height_fraction=1.0\npost_repeat_margin=10\npost_repeat_minimum_calls=2\nedge_identity_margin=10\nedge_height_fraction=1.0\nedge_second_share=0.10\nedge_width_factor=1.5\nedge_spacing_minimum=0.5\nedge_spacing_maximum=1.5\nstart_identity_margin=10\nstart_norm_calls=30\nstart_height_fraction=0.5\n[variant_calling]\nmax_indel_length=50\nread_end_margin=10\n[consensus]\nseeded_run=4\nmaximum_undecided_gap=20\nminimum_resolved_bases=20\nsite_mismatch_cost=2\nsite_gap_open_cost=3\nsite_gap_extension_cost=2\n[review]\nrun_growth=3\nedge_zone=20\nphase_run=4\n",
        )
        .map_err(|source| KernelError::Output {
            path: path.clone(),
            source,
        })?;
        assert!(load_path(&path).is_err());
        Ok(())
    }

    #[test]
    fn loads_the_shipped_configuration_and_its_profile() -> Result<()> {
        let config = load_path(&Path::new(env!("CARGO_MANIFEST_DIR")).join(DEFAULT_CONFIG_PATH))?;
        let profile = dna_kernel::profile::Profile::load(&config.profile_path)?;
        assert_eq!(profile.identity().id, "human-mtdna-rcrs");
        Ok(())
    }
}
