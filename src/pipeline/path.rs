//! CLI/application filesystem naming and publication-target validation.

use std::path::{Path, PathBuf};

use crate::error::{KernelError, Result};

fn validate_output(output: &Path) -> Result<()> {
    if output.exists() {
        return Err(KernelError::Path {
            kind: "output",
            path: output.to_path_buf(),
            reason: "target already exists",
        }
        .into());
    }
    let parent = output
        .parent()
        .filter(|value| !value.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    if parent.exists() && !parent.is_dir() {
        return Err(KernelError::Path {
            kind: "output directory",
            path: parent.to_path_buf(),
            reason: "path exists but is not a directory",
        }
        .into());
    }
    Ok(())
}

/// Returns and validates the deterministic CLI publication path for the
/// analysis of the read named `stem`.
pub(super) fn analysis_output(stem: &str) -> Result<PathBuf> {
    let output = PathBuf::from("results").join(format!("{stem}.json"));
    validate_output(&output)?;
    Ok(output)
}

/// Returns and validates the deterministic CLI publication path for basecalls.
pub(super) fn basecall_output(trace: &Path) -> Result<PathBuf> {
    let output = PathBuf::from("results").join(format!("{}.basecalls.json", trace_stem(trace)?));
    validate_output(&output)?;
    Ok(output)
}

/// Returns and validates the deterministic CLI publication path for a read's
/// evidence.
pub(super) fn evidence_output(trace: &Path) -> Result<PathBuf> {
    let output = PathBuf::from("results").join(format!("{}.evidence.json", trace_stem(trace)?));
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

/// Returns and validates the deterministic CLI publication path for the review list.
pub(super) fn review_output(sample_id: &str) -> Result<PathBuf> {
    let output = PathBuf::from("results").join(format!("{sample_id}.review.json"));
    validate_output(&output)?;
    Ok(output)
}

/// Returns and validates the deterministic CLI publication path for a run signal.
pub(super) fn signal_output(run_id: &str) -> Result<PathBuf> {
    let output = PathBuf::from("results").join(format!("{run_id}.signal.json"));
    validate_output(&output)?;
    Ok(output)
}

/// Returns the validated UTF-8 read stem of a read evidence document, shared
/// by result and log paths: its file stem without the `.evidence` suffix
/// `evidence` writes, so it is the stem of the trace it was written from.
pub(super) fn evidence_stem(evidence: &Path) -> Result<&str> {
    evidence
        .file_stem()
        .and_then(|value| value.to_str())
        .map(|value| value.strip_suffix(".evidence").unwrap_or(value))
        .filter(|value| !value.is_empty())
        .ok_or_else(|| {
            KernelError::Path {
                kind: "read evidence",
                path: evidence.to_path_buf(),
                reason: "file stem must be valid non-empty UTF-8",
            }
            .into()
        })
}

/// Returns the validated UTF-8 trace stem shared by result and log paths.
pub(super) fn trace_stem(trace: &Path) -> Result<&str> {
    trace
        .file_stem()
        .and_then(|value| value.to_str())
        .filter(|value| !value.is_empty())
        .ok_or_else(|| {
            KernelError::Path {
                kind: "AB1",
                path: trace.to_path_buf(),
                reason: "file stem must be valid non-empty UTF-8",
            }
            .into()
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    // SRS-IN-010
    #[test]
    fn an_evidence_document_is_named_after_its_trace() {
        assert_eq!(
            evidence_stem(Path::new("results/A01_HV1F.evidence.json")).ok(),
            Some("A01_HV1F")
        );
        assert_eq!(evidence_stem(Path::new("edited.json")).ok(), Some("edited"));
        assert!(evidence_stem(Path::new(".evidence.json")).is_err());
    }
}
