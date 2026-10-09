//! Sanger input adapter for validated scientific inputs.
//!
//! The adapter owns source-path validation plus source-format decoding for
//! Sanger sequencing inputs. ABIF is the current container format; `.ab1` is a
//! common sequencing sample filename/extension using that container.
//!
//! The adapter does not derive result paths, validate publication targets, open
//! operational logs, or serialize outputs.

use std::path::{Path, PathBuf};

use crate::config::Config;
use crate::error::{Result, SampleError};
use crate::model::reference::Reference;
use crate::model::sanger::Chromatogram;
use crate::profile::Profile;

use super::{load_config, load_profile, load_reference, require_regular_file};

pub(crate) mod abif;

/// Validated paths, configuration, and profile prepared before decoding one analysis trace.
pub(crate) struct PreparedAnalysisInputs {
    config: Config,
    profile: Profile,
    trace_path: PathBuf,
    reference_path: PathBuf,
}

/// Scientific inputs for one reference-guided Sanger analysis.
pub(crate) struct AnalysisInputs {
    pub(crate) config: Config,
    pub(crate) profile: Profile,
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
    pub(crate) profile: Profile,
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
    let profile = load_profile(&config)?;
    Ok(PreparedAnalysisInputs {
        config,
        profile,
        trace_path: trace_path.to_path_buf(),
        reference_path: reference_path.to_path_buf(),
    })
}

/// Decodes one prepared reference-guided Sanger input set.
pub(crate) fn load_analysis(prepared: PreparedAnalysisInputs) -> Result<AnalysisInputs> {
    let PreparedAnalysisInputs {
        config,
        profile,
        trace_path,
        reference_path,
    } = prepared;
    let trace = abif::load(&trace_path)?;
    let reference = load_reference(&reference_path, &profile)?;
    Ok(AnalysisInputs {
        config,
        profile,
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
        return Err(SampleError::NoTraces.into());
    }
    for trace_path in trace_paths {
        require_regular_file(trace_path, "AB1")?;
    }
    require_regular_file(reference_path, "reference")?;
    let config = load_config(config_path)?;
    let profile = load_profile(&config)?;
    let traces = trace_paths
        .iter()
        .map(|path| abif::load(path))
        .collect::<Result<Vec<_>>>()?;
    let reference = load_reference(reference_path, &profile)?;
    Ok(SampleInputs {
        config,
        profile,
        traces,
        reference,
    })
}
