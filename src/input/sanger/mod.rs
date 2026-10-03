//! Sanger input adapter for validated scientific inputs.
//!
//! The adapter owns source-path validation plus source-format decoding for
//! Sanger sequencing inputs. ABIF is the current container format; `.ab1` is a
//! common sequencing sample filename/extension using that container.
//!
//! The adapter does not derive result paths, validate publication targets, open
//! operational logs, or serialize outputs.

use std::path::{Path, PathBuf};

use crate::config::{self, Config};
use crate::error::{Error, Result};
use crate::model::reference::Reference;
use crate::model::trace::Chromatogram;
use crate::reference;

pub(crate) mod abif;

/// Validated paths and configuration prepared before decoding one analysis trace.
pub(crate) struct PreparedAnalysisInputs {
    config: Config,
    trace_path: PathBuf,
    reference_path: PathBuf,
}

/// Scientific inputs for one reference-guided Sanger analysis.
pub(crate) struct AnalysisInputs {
    pub(crate) config: Config,
    pub(crate) trace: Chromatogram,
    pub(crate) reference: Reference,
}

/// Validated path and configuration prepared before decoding one basecall trace.
pub(crate) struct PreparedBasecallInputs {
    config: Config,
    trace_path: PathBuf,
}

/// Scientific inputs for one reference-free Sanger basecall.
pub(crate) struct BasecallInputs {
    pub(crate) config: Config,
    pub(crate) trace: Chromatogram,
}

/// Scientific inputs for one multi-read Sanger sample analysis.
pub(crate) struct SampleInputs {
    pub(crate) config: Config,
    pub(crate) traces: Vec<Chromatogram>,
    pub(crate) reference: Reference,
}

/// Validates source paths and loads configuration for one reference-guided trace.
pub(crate) fn prepare_analysis(
    trace_path: &Path,
    reference_path: &Path,
    config_path: &Path,
) -> Result<PreparedAnalysisInputs> {
    require_regular_file(trace_path, "AB1")?;
    require_regular_file(reference_path, "reference")?;
    let config = load_config(config_path)?;
    Ok(PreparedAnalysisInputs {
        config,
        trace_path: trace_path.to_path_buf(),
        reference_path: reference_path.to_path_buf(),
    })
}

/// Decodes one prepared reference-guided Sanger input set.
pub(crate) fn load_analysis(prepared: PreparedAnalysisInputs) -> Result<AnalysisInputs> {
    let PreparedAnalysisInputs {
        config,
        trace_path,
        reference_path,
    } = prepared;
    let trace = abif::load(&trace_path)?;
    let reference = reference::load(&reference_path, config.reference.topology)?;
    Ok(AnalysisInputs {
        config,
        trace,
        reference,
    })
}

/// Validates the source path and loads configuration for one basecall trace.
pub(crate) fn prepare_basecall(
    trace_path: &Path,
    config_path: &Path,
) -> Result<PreparedBasecallInputs> {
    require_regular_file(trace_path, "AB1")?;
    let config = load_config(config_path)?;
    Ok(PreparedBasecallInputs {
        config,
        trace_path: trace_path.to_path_buf(),
    })
}

/// Decodes one prepared reference-free Sanger input set.
pub(crate) fn load_basecall(prepared: PreparedBasecallInputs) -> Result<BasecallInputs> {
    let PreparedBasecallInputs { config, trace_path } = prepared;
    let trace = abif::load(&trace_path)?;
    Ok(BasecallInputs { config, trace })
}

/// Validates and decodes one or more Sanger traces against one shared reference.
pub(crate) fn load_sample(
    trace_paths: &[PathBuf],
    reference_path: &Path,
    config_path: &Path,
) -> Result<SampleInputs> {
    if trace_paths.is_empty() {
        return Err(Error::Sample(
            "sample analysis requires at least one AB1 trace".into(),
        ));
    }
    for trace_path in trace_paths {
        require_regular_file(trace_path, "AB1")?;
    }
    require_regular_file(reference_path, "reference")?;
    let config = load_config(config_path)?;
    let traces = trace_paths
        .iter()
        .map(|path| abif::load(path))
        .collect::<Result<Vec<_>>>()?;
    let reference = reference::load(reference_path, config.reference.topology)?;
    Ok(SampleInputs {
        config,
        traces,
        reference,
    })
}

fn load_config(path: &Path) -> Result<Config> {
    let config = config::load_path(path)?;
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
