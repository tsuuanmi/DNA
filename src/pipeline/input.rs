//! CLI/application path naming and publication-target validation.

use std::path::{Path, PathBuf};

use crate::error::{Error, Result};

fn validate_output(output: &Path) -> Result<()> {
    if output.exists() {
        return Err(Error::Path {
            kind: "output",
            path: output.to_path_buf(),
            reason: "target already exists".into(),
        });
    }
    let parent = output
        .parent()
        .filter(|value| !value.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    if parent.exists() && !parent.is_dir() {
        return Err(Error::Path {
            kind: "output directory",
            path: parent.to_path_buf(),
            reason: "path exists but is not a directory".into(),
        });
    }
    Ok(())
}

/// Returns and validates the deterministic CLI publication path for analysis.
pub(crate) fn analysis_output(trace: &Path) -> Result<PathBuf> {
    let output = PathBuf::from("results").join(format!("{}.json", trace_stem(trace)?));
    validate_output(&output)?;
    Ok(output)
}

/// Returns and validates the deterministic CLI publication path for basecalls.
pub(crate) fn basecall_output(trace: &Path) -> Result<PathBuf> {
    let output = PathBuf::from("results").join(format!("{}.basecalls.json", trace_stem(trace)?));
    validate_output(&output)?;
    Ok(output)
}

/// Returns and validates the deterministic CLI publication path for sample evidence.
pub(crate) fn sample_output(sample_id: &str) -> Result<PathBuf> {
    let output = PathBuf::from("results").join(format!("{sample_id}.sample.json"));
    validate_output(&output)?;
    Ok(output)
}

/// Returns the validated UTF-8 trace stem shared by result and log paths.
pub(crate) fn trace_stem(trace: &Path) -> Result<&str> {
    trace
        .file_stem()
        .and_then(|value| value.to_str())
        .filter(|value| !value.is_empty())
        .ok_or_else(|| Error::Path {
            kind: "AB1",
            path: trace.to_path_buf(),
            reason: "file stem must be valid non-empty UTF-8".into(),
        })
}

/// Validates the sample identifier used for deterministic result/log names.
pub(crate) fn validate_sample_id(sample_id: &str) -> Result<()> {
    let mut characters = sample_id.chars();
    let valid_first = characters
        .next()
        .is_some_and(|value| value.is_ascii_alphanumeric());
    let valid_rest =
        characters.all(|value| value.is_ascii_alphanumeric() || matches!(value, '_' | '.' | '-'));
    if sample_id.len() > 128 || !valid_first || !valid_rest {
        return Err(Error::Sample(
            "sample id must be 1..=128 ASCII characters, start with an alphanumeric character, and contain only alphanumeric, '_', '.', or '-'".into(),
        ));
    }
    Ok(())
}
