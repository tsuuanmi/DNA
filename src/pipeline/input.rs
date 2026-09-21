//! Path validation and command-specific input loading.

use std::path::{Path, PathBuf};

use crate::cli::{AnalyzeArgs, BasecallArgs, SampleArgs};
use crate::config::{self, Config};
use crate::error::{Error, Result};
use crate::model::reference::Reference;
use crate::model::trace::Chromatogram;
use crate::{reference, trace};

/// Inputs for one reference-guided analysis.
pub(crate) struct AnalysisInputs {
    pub(crate) config: Config,
    pub(crate) trace: Chromatogram,
    pub(crate) reference: Reference,
    pub(crate) output: PathBuf,
}

/// Inputs for one reference-free basecall operation.
pub(crate) struct BasecallInputs {
    pub(crate) config: Config,
    pub(crate) trace: Chromatogram,
    pub(crate) output: PathBuf,
}

/// Inputs for one multi-read sample evidence operation.
pub(crate) struct SampleInputs {
    pub(crate) config: Config,
    pub(crate) traces: Vec<Chromatogram>,
    pub(crate) reference: Reference,
}

/// Validates and loads one trace, one reference, and one configuration.
pub(crate) fn load_analysis(args: &AnalyzeArgs) -> Result<AnalysisInputs> {
    require_regular_file(&args.trace, "AB1")?;
    require_regular_file(&args.reference, "reference")?;
    let config = load_config()?;
    let output = analysis_output_path(&args.trace)?;
    validate_output(&output)?;
    let trace = trace::load(&args.trace)?;
    let reference = reference::load(&args.reference, config.reference.topology)?;
    Ok(AnalysisInputs {
        config,
        trace,
        reference,
        output,
    })
}

/// Validates and loads one trace and one configuration without a reference.
pub(crate) fn load_basecall(args: &BasecallArgs) -> Result<BasecallInputs> {
    require_regular_file(&args.trace, "AB1")?;
    let config = load_config()?;
    let output = basecall_output_path(&args.trace)?;
    validate_output(&output)?;
    let trace = trace::load(&args.trace)?;
    Ok(BasecallInputs {
        config,
        trace,
        output,
    })
}

/// Validates and loads one or more sample traces against one shared reference.
pub(crate) fn load_sample(args: &SampleArgs) -> Result<SampleInputs> {
    validate_sample_id(&args.sample_id)?;
    if args.traces.is_empty() {
        return Err(Error::Sample(
            "sample analysis requires at least one AB1 trace".into(),
        ));
    }
    for trace_path in &args.traces {
        require_regular_file(trace_path, "AB1")?;
    }
    require_regular_file(&args.reference, "reference")?;
    let config = load_config()?;
    let traces = args
        .traces
        .iter()
        .map(|path| trace::load(path))
        .collect::<Result<Vec<_>>>()?;
    let reference = reference::load(&args.reference, config.reference.topology)?;
    Ok(SampleInputs {
        config,
        traces,
        reference,
    })
}

fn load_config() -> Result<Config> {
    let config = config::load()?;
    require_regular_file(&config.source_path, "configuration")?;
    Ok(config)
}

fn require_regular_file(path: &Path, kind: &'static str) -> Result<()> {
    let metadata = path.metadata().map_err(|source| Error::Read {
        kind,
        path: path.to_path_buf(),
        source,
    })?;
    if !metadata.is_file() || metadata.len() == 0 {
        return Err(Error::Path {
            kind,
            path: path.to_path_buf(),
            reason: "path must be a non-empty regular file".into(),
        });
    }
    Ok(())
}

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

fn analysis_output_path(trace: &Path) -> Result<PathBuf> {
    Ok(PathBuf::from("results").join(format!("{}.json", trace_stem(trace)?)))
}

fn basecall_output_path(trace: &Path) -> Result<PathBuf> {
    Ok(PathBuf::from("results").join(format!("{}.basecalls.json", trace_stem(trace)?)))
}

pub(crate) fn sample_output(sample_id: &str) -> Result<PathBuf> {
    let output = PathBuf::from("results").join(format!("{sample_id}.sample.json"));
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
            reason: "file stem must be valid non-empty UTF-8".into(),
        })
}

/// Validates the sample identifier used for deterministic result/log names.
pub(super) fn validate_sample_id(sample_id: &str) -> Result<()> {
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
