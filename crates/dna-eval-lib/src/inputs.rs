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
//! reference = { path = "rCRS.fasta", sha256 = "…" }
//! recheck = { path = "recheck.tsv", sha256 = "…" }
//! ```
//!
//! `manifest`, `truth`, `config`, and `reference` are required; `recheck`
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
}

/// The inputs file as written.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawInputs {
    manifest: Pinned,
    truth: Pinned,
    config: Pinned,
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
        let pinned = |name: &str, pinned: Pinned| -> Result<PathBuf> {
            let file = resolve(pinned.path);
            let found = sha256_file(&file)
                .map_err(|error| invalid(format!("{name} {}: {error}", file.display())))?;
            if found != pinned.sha256 {
                return Err(invalid(format!(
                    "{name} {} has SHA-256 {found}, not the pinned {}",
                    file.display(),
                    pinned.sha256
                )));
            }
            Ok(file)
        };
        let current_since = raw
            .current_since
            .map(|month| month.parse::<Month>())
            .transpose()
            .map_err(invalid)?;
        let sha256 =
            sha256_file(path).map_err(|error| invalid(format!("{}: {error}", path.display())))?;
        Ok(Self {
            manifest: pinned("manifest", raw.manifest)?,
            truth: pinned("truth", raw.truth)?,
            config: pinned("config", raw.config)?,
            reference: pinned("reference", raw.reference)?,
            recheck: raw
                .recheck
                .map(|file| pinned("recheck", file))
                .transpose()?,
            forced_windows: raw
                .forced_windows
                .map(|file| pinned("forced_windows", file))
                .transpose()?,
            current_since,
            evidence_cache: raw.evidence_cache.map(resolve),
            run_context: raw.run_context.map(resolve),
            file: Some((path.to_path_buf(), sha256)),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write(path: &Path, text: &str) -> std::result::Result<(), Box<dyn std::error::Error>> {
        Ok(fs::write(path, text)?)
    }

    /// Pinned files resolve against the inputs file's directory and load
    /// when their content is the pinned one; a changed file, an unknown key,
    /// or an invalid month is refused.
    #[test]
    fn loads_pinned_inputs_and_refuses_a_changed_file()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        let root = tempfile::tempdir()?;
        let directory = root.path();
        for name in [
            "panel.tsv",
            "truth.sqlite",
            "dna.toml",
            "reference.fa",
            "recheck.tsv",
        ] {
            write(&directory.join(name), name)?;
        }
        let pin = |name: &str| -> std::result::Result<String, Box<dyn std::error::Error>> {
            Ok(format!(
                "{{ path = \"{name}\", sha256 = \"{}\" }}",
                sha256_file(&directory.join(name))?
            ))
        };
        let text = format!(
            "current_since = \"2026-01\"\nevidence_cache = \"cache\"\nmanifest = {}\ntruth = {}\nconfig = {}\nreference = {}\nrecheck = {}\n",
            pin("panel.tsv")?,
            pin("truth.sqlite")?,
            pin("dna.toml")?,
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
        assert_eq!(
            inputs.file.as_ref().map(|(_, sha256)| sha256.clone()),
            Some(sha256_file(&file)?)
        );

        write(&directory.join("recheck.tsv"), "edited")?;
        let changed = Inputs::load(&file).map_err(|error| error.to_string());
        assert!(
            changed
                .as_ref()
                .is_err_and(|message| message.contains("recheck") && message.contains("pinned")),
            "{changed:?}"
        );
        write(&directory.join("recheck.tsv"), "recheck.tsv")?;
        for edited in [
            text.replace("current_since", "since"),
            text.replace("2026-01", "2026-13"),
            text.replace("manifest = ", "panel = "),
        ] {
            write(&file, &edited)?;
            assert!(
                matches!(Inputs::load(&file), Err(Error::Inputs { .. })),
                "{edited}"
            );
        }
        Ok(())
    }
}
