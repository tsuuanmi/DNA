//! CLI/application filesystem naming and publication-target validation.

use std::path::{Path, PathBuf};

use dna_kernel::error::{Error, Result};

fn validate_output(output: &Path) -> Result<()> {
    if output.exists() {
        return Err(Error::Path {
            kind: "output",
            path: output.to_path_buf(),
            reason: "target already exists",
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
            reason: "path exists but is not a directory",
        });
    }
    Ok(())
}

/// Returns and validates the deterministic CLI publication path for analysis.
pub(super) fn analysis_output(trace: &Path) -> Result<PathBuf> {
    let output = PathBuf::from("results").join(format!("{}.json", trace_stem(trace)?));
    validate_output(&output)?;
    Ok(output)
}

/// Returns and validates the deterministic CLI publication path for basecalls.
pub(super) fn basecall_output(trace: &Path) -> Result<PathBuf> {
    let output = PathBuf::from("results").join(format!("{}.basecalls.json", trace_stem(trace)?));
    validate_output(&output)?;
    Ok(output)
}

/// Returns and validates the deterministic CLI publication path for sample evidence.
pub(super) fn sample_output(sample_id: &str) -> Result<PathBuf> {
    let output = PathBuf::from("results").join(format!("{sample_id}.sample.json"));
    validate_output(&output)?;
    Ok(output)
}

/// Returns and validates the deterministic CLI publication paths for the
/// consensus document and its FASTA.
pub(super) fn consensus_outputs(sample_id: &str) -> Result<(PathBuf, PathBuf)> {
    let document = PathBuf::from("results").join(format!("{sample_id}.consensus.json"));
    let fasta = PathBuf::from("results").join(format!("{sample_id}.consensus.fasta"));
    validate_output(&document)?;
    validate_output(&fasta)?;
    Ok((document, fasta))
}

/// Returns and validates the deterministic CLI publication path for the variants document.
pub(super) fn call_output(sample_id: &str) -> Result<PathBuf> {
    let output = PathBuf::from("results").join(format!("{sample_id}.variants.json"));
    validate_output(&output)?;
    Ok(output)
}

/// Returns and validates the deterministic CLI publication path for notation.
pub(super) fn notation_output(sample_id: &str) -> Result<PathBuf> {
    let output = PathBuf::from("results").join(format!("{sample_id}.notation.json"));
    validate_output(&output)?;
    Ok(output)
}

/// Returns the validated UTF-8 trace stem shared by result and log paths.
pub(super) fn trace_stem(trace: &Path) -> Result<&str> {
    trace
        .file_stem()
        .and_then(|value| value.to_str())
        .filter(|value| !value.is_empty())
        .ok_or_else(|| Error::Path {
            kind: "AB1",
            path: trace.to_path_buf(),
            reason: "file stem must be valid non-empty UTF-8",
        })
}
