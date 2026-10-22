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
    use crate::error::Error;

    // SRS-CFG-003
    /// An unknown key is refused at the root and in every section of the
    /// shipped configuration, as an unknown field naming it, while the
    /// unchanged configuration loads.
    #[test]
    fn rejects_unknown_keys() -> Result<()> {
        const UNKNOWN: &str = "unknown_key_for_test";
        let shipped = Path::new(env!("CARGO_MANIFEST_DIR")).join(DEFAULT_CONFIG_PATH);
        let text = fs::read_to_string(&shipped).map_err(|source| KernelError::Read {
            kind: "configuration",
            path: shipped.clone(),
            source,
        })?;
        let directory = tempdir().map_err(|source| KernelError::Output {
            path: PathBuf::from("temporary directory"),
            source,
        })?;
        let path = directory.path().join("dna.toml");
        let write = |text: &str| {
            fs::write(&path, text).map_err(|source| KernelError::Output {
                path: path.clone(),
                source,
            })
        };
        write(&text)?;
        load_path(&path)?;
        let lines = text.lines().collect::<Vec<_>>();
        let mut placements = vec![0];
        placements.extend(
            lines
                .iter()
                .enumerate()
                .filter(|(_, line)| line.starts_with('['))
                .map(|(index, _)| index + 1),
        );
        assert!(placements.len() > 10, "every section is exercised");
        for at in placements {
            let mut edited = lines.clone();
            let key = format!("{UNKNOWN} = 1");
            edited.insert(at, key.as_str());
            write(&(edited.join("\n") + "\n"))?;
            let result = load_path(&path);
            let Err(Error::Kernel(KernelError::ConfigParse { source, .. })) = &result else {
                unreachable!("an unknown key at line {at} must be refused: {result:?}");
            };
            let message = source.to_string();
            assert!(
                message.contains("unknown field") && message.contains(UNKNOWN),
                "line {at}: {message}"
            );
        }
        Ok(())
    }

    // SRS-CFG-001 SRS-CFG-002
    #[test]
    fn loads_the_shipped_configuration_and_its_profile() -> Result<()> {
        let config = load_path(&Path::new(env!("CARGO_MANIFEST_DIR")).join(DEFAULT_CONFIG_PATH))?;
        let profile = dna_kernel::profile::Profile::load(&config.profile_path)?;
        assert_eq!(profile.identity().id, "human-mtdna-rcrs");
        Ok(())
    }
}
