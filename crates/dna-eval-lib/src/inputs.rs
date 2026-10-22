//! `run --inputs <file>`: a run's inputs from one TOML file, every input
//! file pinned by its SHA-256, so a probe names its panel, truth,
//! configuration, reference, and lists once and a changed or replaced file
//! is refused before anything runs.
//!
//! ```toml
//! current_since = "2026-01"
//! evidence_cache = "/data/evidence-cache"
//! manifest = { path = "test-panel/dev-v5.tsv", sha256 = "…" }
//! truth = { path = "pins/statistics.sqlite", sha256 = "…" }
//! config = { path = "dna.toml", sha256 = "…" }
//! profile = { path = "profiles/human-mtdna-rcrs.toml", sha256 = "…" }
//! reference = { path = "rCRS.fasta", sha256 = "…" }
//! recheck = { path = "recheck.tsv", sha256 = "…" }
//! ```
//!
//! `manifest`, `truth`, `config`, `profile`, and `reference` are required;
//! `profile` must be the file the configuration's `profile` key names. A pin
//! is 64 lowercase hexadecimal digits, and every pinned file is checked
//! again before the run writes `run.json` (`Inputs::verify`). `recheck`
//! and `forced_windows` are optional pinned files; `current_since`,
//! `evidence_cache`, and `run_context` are optional and not pinned (the
//! caches key themselves, and a run context directory is read by name). A
//! relative path is resolved against the file's directory. The file supplies
//! these inputs alone: the command line gives none of them beside it.

use std::fs;
use std::path::{Path, PathBuf};

use serde::Deserialize;

use crate::error::{Error, Result};
use crate::manifest::Month;
use crate::profile_path;
use crate::runner::sha256_file;

/// Largest inputs file read.
const MAX_INPUTS_BYTES: u64 = 64 * 1024;

/// A run's inputs, from the command line or an inputs file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Inputs {
    pub(crate) manifest: PathBuf,
    pub(crate) truth: PathBuf,
    pub(crate) config: PathBuf,
    pub(crate) reference: PathBuf,
    pub(crate) recheck: Option<PathBuf>,
    pub(crate) forced_windows: Option<PathBuf>,
    pub(crate) current_since: Option<Month>,
    pub(crate) evidence_cache: Option<PathBuf>,
    pub(crate) run_context: Option<PathBuf>,
    /// The inputs file and its SHA-256, when the inputs came from one.
    pub(crate) file: Option<(PathBuf, String)>,
    /// Every pinned file: its name in the inputs file, path, and SHA-256.
    pub(crate) pins: Vec<(&'static str, PathBuf, String)>,
}

/// The inputs file as written.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawInputs {
    manifest: Pinned,
    truth: Pinned,
    config: Pinned,
    profile: Pinned,
    reference: Pinned,
    recheck: Option<Pinned>,
    forced_windows: Option<Pinned>,
    current_since: Option<String>,
    evidence_cache: Option<PathBuf>,
    run_context: Option<PathBuf>,
}

/// One input file and the SHA-256 its content must have.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct Pinned {
    path: PathBuf,
    sha256: String,
}

impl Inputs {
    /// The inputs the file at `path` names, each pinned file checked
    /// against its SHA-256.
    ///
    /// # Errors
    ///
    /// [`Error::Inputs`] when the file cannot be read or parsed, names an
    /// invalid month, or pins a file whose content has another SHA-256
    /// (or that cannot be read).
    pub(crate) fn load(path: &Path) -> Result<Self> {
        let invalid = |message: String| Error::Inputs {
            path: path.to_path_buf(),
            message,
        };
        let bytes = fs::metadata(path)
            .map_err(|error| invalid(error.to_string()))
            .and_then(|metadata| {
                if metadata.len() > MAX_INPUTS_BYTES {
                    Err(invalid(format!(
                        "{} bytes, more than the maximum of {MAX_INPUTS_BYTES}",
                        metadata.len()
                    )))
                } else {
                    fs::read(path).map_err(|error| invalid(error.to_string()))
                }
            })?;
        let text = std::str::from_utf8(&bytes).map_err(|error| invalid(error.to_string()))?;
        let raw: RawInputs = toml::from_str(text).map_err(|error| invalid(error.to_string()))?;
        let directory = path.parent().unwrap_or_else(|| Path::new(""));
        let resolve = |path: PathBuf| directory.join(path);
        let mut pins = Vec::new();
        let mut pinned = |name: &'static str, pinned: Pinned| -> Result<PathBuf> {
            if !is_sha256(&pinned.sha256) {
                return Err(invalid(format!(
                    "{name}: sha256 {:?} is not 64 lowercase hexadecimal digits",
                    pinned.sha256
                )));
            }
            let file = resolve(pinned.path);
            check(name, &file, &pinned.sha256).map_err(invalid)?;
            pins.push((name, file.clone(), pinned.sha256));
            Ok(file)
        };
        let manifest = pinned("manifest", raw.manifest)?;
        let truth = pinned("truth", raw.truth)?;
        let config = pinned("config", raw.config)?;
        let profile = pinned("profile", raw.profile)?;
        let reference = pinned("reference", raw.reference)?;
        let recheck = raw
            .recheck
            .map(|file| pinned("recheck", file))
            .transpose()?;
        let forced_windows = raw
            .forced_windows
            .map(|file| pinned("forced_windows", file))
            .transpose()?;
        let named = profile_path(&config).map_err(|error| invalid(error.to_string()))?;
        let same = |a: &Path, b: &Path| match (fs::canonicalize(a), fs::canonicalize(b)) {
            (Ok(a), Ok(b)) => a == b,
            _ => false,
        };
        if !same(&profile, &named) {
            return Err(invalid(format!(
                "profile {} is not the configuration's profile {}",
                profile.display(),
                named.display()
            )));
        }
        let current_since = raw
            .current_since
            .map(|month| month.parse::<Month>())
            .transpose()
            .map_err(invalid)?;
        let sha256 =
            sha256_file(path).map_err(|error| invalid(format!("{}: {error}", path.display())))?;
        Ok(Self {
            manifest,
            truth,
            config,
            reference,
            recheck,
            forced_windows,
            current_since,
            evidence_cache: raw.evidence_cache.map(resolve),
            run_context: raw.run_context.map(resolve),
            file: Some((path.to_path_buf(), sha256)),
            pins,
        })
    }

    /// Checks every pinned file again, so `run.json`'s claim that each had
    /// its pinned SHA-256 holds for the whole run.
    ///
    /// # Errors
    ///
    /// [`Error::Inputs`] when a pinned file changed or cannot be read.
    pub(crate) fn verify(&self) -> Result<()> {
        let Some((path, _)) = &self.file else {
            return Ok(());
        };
        for (name, file, sha256) in &self.pins {
            check(name, file, sha256).map_err(|message| Error::Inputs {
                path: path.clone(),
                message: format!("{message}: it changed during the run"),
            })?;
        }
        Ok(())
    }
}

/// Fails with a message unless the file at `file` has SHA-256 `sha256`.
fn check(name: &str, file: &Path, sha256: &str) -> std::result::Result<(), String> {
    let found = sha256_file(file).map_err(|error| format!("{name} {}: {error}", file.display()))?;
    if found == sha256 {
        Ok(())
    } else {
        Err(format!(
            "{name} {} has SHA-256 {found}, not the pinned {sha256}",
            file.display()
        ))
    }
}

/// Whether `text` is a SHA-256 as pins write it: 64 lowercase hexadecimal
/// digits.
fn is_sha256(text: &str) -> bool {
    text.len() == 64
        && text
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write(path: &Path, text: &str) -> std::result::Result<(), Box<dyn std::error::Error>> {
        Ok(fs::write(path, text)?)
    }

    /// Pinned files resolve against the inputs file's directory and load
    /// when their content is the pinned one; a changed file (at load or
    /// before `run.json`), another profile than the configuration's, a pin
    /// that is not a SHA-256, an unknown key, or an invalid month is
    /// refused.
    #[test]
    fn loads_pinned_inputs_and_refuses_a_changed_file()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        let root = tempfile::tempdir()?;
        let directory = root.path();
        fs::create_dir(directory.join("profiles"))?;
        for name in [
            "panel.tsv",
            "truth.sqlite",
            "reference.fa",
            "recheck.tsv",
            "profiles/target.toml",
            "other.toml",
        ] {
            write(&directory.join(name), name)?;
        }
        write(
            &directory.join("dna.toml"),
            "schema_version = 1\nprofile = \"profiles/target.toml\"\n",
        )?;
        let pin = |name: &str| -> std::result::Result<String, Box<dyn std::error::Error>> {
            Ok(format!(
                "{{ path = \"{name}\", sha256 = \"{}\" }}",
                sha256_file(&directory.join(name))?
            ))
        };
        let text = format!(
            "current_since = \"2026-01\"\nevidence_cache = \"cache\"\nmanifest = {}\ntruth = {}\nconfig = {}\nprofile = {}\nreference = {}\nrecheck = {}\n",
            pin("panel.tsv")?,
            pin("truth.sqlite")?,
            pin("dna.toml")?,
            pin("profiles/target.toml")?,
            pin("reference.fa")?,
            pin("recheck.tsv")?,
        );
        let file = directory.join("inputs.toml");
        write(&file, &text)?;
        let inputs = Inputs::load(&file)?;
        assert_eq!(inputs.manifest, directory.join("panel.tsv"));
        assert_eq!(inputs.recheck, Some(directory.join("recheck.tsv")));
        assert_eq!(inputs.evidence_cache, Some(directory.join("cache")));
        assert_eq!(inputs.current_since, Some("2026-01".parse()?));
        assert_eq!(inputs.forced_windows, None);
        assert_eq!(inputs.pins.len(), 6);
        assert_eq!(
            inputs.file.as_ref().map(|(_, sha256)| sha256.clone()),
            Some(sha256_file(&file)?)
        );
        assert!(inputs.verify().is_ok());

        write(&directory.join("recheck.tsv"), "edited")?;
        let during = inputs.verify().map_err(|error| error.to_string());
        assert!(
            during.as_ref().is_err_and(
                |message| message.contains("recheck") && message.contains("during the run")
            ),
            "{during:?}"
        );
        let changed = Inputs::load(&file).map_err(|error| error.to_string());
        assert!(
            changed
                .as_ref()
                .is_err_and(|message| message.contains("recheck") && message.contains("pinned")),
            "{changed:?}"
        );
        write(&directory.join("recheck.tsv"), "recheck.tsv")?;
        let sha256 = sha256_file(&directory.join("panel.tsv"))?;
        let other = text.replace(&pin("profiles/target.toml")?, &pin("other.toml")?);
        for (edited, message) in [
            (text.replace("current_since", "since"), "unknown field"),
            (text.replace("2026-01", "2026-13"), "month"),
            (text.replace("manifest = ", "panel = "), "unknown field"),
            (text.replace(&sha256, &sha256.to_uppercase()), "lowercase"),
            (text.replace(&sha256, &sha256[..63]), "lowercase"),
            (other, "not the configuration's profile"),
        ] {
            write(&file, &edited)?;
            let result = Inputs::load(&file).map_err(|error| error.to_string());
            assert!(
                result.as_ref().is_err_and(|error| error.contains(message)),
                "{edited}: {result:?}"
            );
        }
        Ok(())
    }
}
