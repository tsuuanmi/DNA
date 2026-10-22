//! Configuration path resolution, parsing, and identity.

use std::collections::BTreeMap;
use std::env;
use std::fmt::Write as _;
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
    raw.validate(
        path.to_path_buf(),
        hex_sha256(&bytes),
        section_sha256(path, text)?,
    )
}

/// Each top-level table's digest: the SHA-256 of its canonical form
/// ([`canonical`]), so comments, layout, and key order do not count.
fn section_sha256(path: &Path, text: &str) -> Result<BTreeMap<String, String>> {
    let table: toml::Table = toml::from_str(text).map_err(|source| KernelError::ConfigParse {
        path: path.to_path_buf(),
        source: Box::new(source),
    })?;
    Ok(table
        .iter()
        .filter(|(_, value)| value.is_table())
        .map(|(name, value)| {
            let mut written = String::new();
            canonical(value, &mut written);
            (name.clone(), hex_sha256(written.as_bytes()))
        })
        .collect())
}

/// One value written in a canonical form: tables with their keys in order,
/// strings quoted and escaped, floats in Rust's shortest round-trip form, so
/// equal values are written alike and a float never reads as an integer.
fn canonical(value: &toml::Value, out: &mut String) {
    match value {
        toml::Value::String(text) => {
            let _ = write!(out, "{text:?}");
        }
        toml::Value::Integer(number) => out.push_str(&number.to_string()),
        toml::Value::Float(number) => {
            let _ = write!(out, "{number:?}");
        }
        toml::Value::Boolean(flag) => out.push_str(if *flag { "true" } else { "false" }),
        toml::Value::Datetime(time) => out.push_str(&time.to_string()),
        toml::Value::Array(items) => {
            out.push('[');
            for (index, item) in items.iter().enumerate() {
                if index > 0 {
                    out.push(',');
                }
                canonical(item, out);
            }
            out.push(']');
        }
        toml::Value::Table(table) => {
            out.push('{');
            let mut keys = table.keys().collect::<Vec<_>>();
            keys.sort();
            for (index, key) in keys.into_iter().enumerate() {
                if index > 0 {
                    out.push(',');
                }
                let _ = write!(out, "{key:?}=");
                if let Some(item) = table.get(key) {
                    canonical(item, out);
                }
            }
            out.push('}');
        }
    }
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

    /// A plugin's configuration digest follows its own sections' values
    /// alone: comments, layout, the way a number is written, and other
    /// plugins' sections leave it unchanged; any value of its sections
    /// changes it.
    #[test]
    fn plugin_digests_follow_their_sections_values() -> Result<()> {
        use dna_sanger::read_processing::PLUGIN as SANGER;

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
        let load = |name: &str, text: &str| {
            let path = directory.path().join(name);
            fs::write(&path, text).map_err(|source| KernelError::Output {
                path: path.clone(),
                source,
            })?;
            load_path(&path)
        };
        let base = load("base.toml", &text)?;
        let edited = |from: &str, to: &str| {
            assert!(text.contains(from), "{from}");
            text.replacen(from, to, 1)
        };
        let alike = [
            edited("[basecalling]", "# a comment\n\n[basecalling]  # another"),
            edited(
                "secondary_peak_ratio = 0.33",
                "secondary_peak_ratio = 3.3e-1",
            ),
            edited("edge_zone = 40", "edge_zone = 41"),
        ];
        for (index, text) in alike.iter().enumerate() {
            let config = load(&format!("alike{index}.toml"), text)?;
            assert_ne!(config.source_sha256, base.source_sha256);
            assert_eq!(
                config.plugin_sha256(&SANGER),
                base.plugin_sha256(&SANGER),
                "{index}"
            );
        }
        let review = load("review.toml", &alike[2])?;
        assert_ne!(
            review.section_sha256.get("review"),
            base.section_sha256.get("review")
        );
        let sanger = load(
            "sanger.toml",
            &edited("secondary_peak_ratio = 0.33", "secondary_peak_ratio = 0.34"),
        )?;
        assert_ne!(sanger.plugin_sha256(&SANGER), base.plugin_sha256(&SANGER));
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
