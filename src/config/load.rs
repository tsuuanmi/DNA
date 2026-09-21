//! Configuration path resolution, parsing, and identity.

use std::env;
use std::fs;
use std::path::{Path, PathBuf};

use crate::checksum::hex_sha256;
use crate::config::defaults::{DEFAULT_CONFIG_PATH, MAX_CONFIG_BYTES};
use crate::config::types::{Config, RawConfig};
use crate::error::{Error, Result};

/// Loads and validates the one authoritative configuration file.
pub(crate) fn load() -> Result<Config> {
    let path = env::var_os("DNA_CONFIG")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(DEFAULT_CONFIG_PATH));
    load_path(&path)
}

fn load_path(path: &Path) -> Result<Config> {
    let metadata = fs::metadata(path).map_err(|source| Error::Read {
        kind: "configuration",
        path: path.to_path_buf(),
        source,
    })?;
    if metadata.len() > MAX_CONFIG_BYTES as u64 {
        return Err(Error::Config(format!(
            "configuration file exceeds {MAX_CONFIG_BYTES} bytes"
        )));
    }
    let bytes = fs::read(path).map_err(|source| Error::Read {
        kind: "configuration",
        path: path.to_path_buf(),
        source,
    })?;
    if bytes.len() > MAX_CONFIG_BYTES {
        return Err(Error::Config(format!(
            "configuration file exceeds {MAX_CONFIG_BYTES} bytes"
        )));
    }
    let text = std::str::from_utf8(&bytes)
        .map_err(|error| Error::Config(format!("configuration must be UTF-8: {error}")))?;
    let raw: RawConfig = toml::from_str(text).map_err(|source| Error::ConfigParse {
        path: path.to_path_buf(),
        source,
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
        let directory = tempdir().map_err(|source| Error::Output {
            path: PathBuf::from("temporary directory"),
            source,
        })?;
        let path = directory.path().join("dna.toml");
        fs::write(
            &path,
            "schema_version=5\nunknown=true\n[reference]\ntopology='circular'\n[basecalling]\nsecondary_peak_ratio=0.33\n[signal_processing]\nwindow_size_bases=10\nminimum_primary_snr=3.0\nminimum_noisy_windows=2\n[quality_control]\ntrim_window_size=10\nbest_section_fraction=0.1\nmax_relative_quality_score=60\ntrim_stringency=7.0\nminimum_retained_bases=20\n[alignment]\nmatch_score=3\nmismatch_score=-5\nambiguous_score=0\ngap_open_score=-10\ngap_extension_score=-4\nminimum_callable_bases=20\nminimum_identity=0.8\n[sample_reconciliation]\nminimum_comparable_bases=25\nminimum_overlap_agreement=0.5\n[variant_calling]\nmax_indel_length=50\nminimum_peak_height=150\nrelative_quality_threshold=30\nregions=[[1, 50000]]\n",
        )
        .map_err(|source| Error::Output {
            path: path.clone(),
            source,
        })?;
        assert!(load_path(&path).is_err());
        Ok(())
    }
}
