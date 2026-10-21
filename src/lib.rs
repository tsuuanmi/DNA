//! Library boundary for DNA operations.
//!
//! The source-module graph is routed from `src/README.md` and `crates/README.md`. `lib.rs` remains
//! the minimal dispatcher: it exposes stable CLI and error boundaries, translates
//! parsed command arguments into application inputs, and keeps configuration,
//! scientific capabilities, orchestration, and reporting behind explicit module
//! boundaries.

pub mod cli;
mod config;
mod input;
mod model;
mod operation_log;
mod pipeline;
pub mod recorded;
mod report;
pub mod sample_analysis;
pub mod variant_analysis;

pub mod error;

pub use dna_kernel::variant;

/// Target profiles: knowledge about one sequencing target.
pub mod profile {
    pub use dna_kernel::profile::{Profile, ProfileIdentity};
}

/// Optional haplotype-preserving representation normalization of called
/// variants.
pub mod variant_normalization {
    pub use dna_normalize::variant_normalization::{
        NormalizationPolicy, VariantNormalizationResult, normalize,
    };
}

/// Optional profile-driven target nomenclature of normalized variants.
pub mod variant_nomenclature {
    pub use dna_normalize::variant_nomenclature::{
        NomenclatureInput, VariantNomenclatureResult, apply, from_normalization,
    };
}

#[cfg(feature = "fuzzing")]
#[doc(hidden)]
pub mod fuzzing {
    /// Exercises the bounds-checked ABIF directory parser without filesystem I/O.
    pub fn parse_abif(bytes: &[u8]) {
        let _ = dna_sanger::abif::parse_container(bytes.to_vec());
    }
}

use cli::{Cli, Command};
use error::{KernelError, Result};
use std::path::PathBuf;

/// Dispatches a parsed command through the application boundary.
///
/// # Errors
///
/// Returns [`Error`](error::Error) when the operation fails; the failure is also
/// recorded in the operation log unless the log itself cannot be written.
pub fn run(cli: Cli) -> Result<()> {
    let config_path = config::resolve_path();
    match cli.command {
        Command::Analyze(args) => match (&args.trace, &args.evidence) {
            (Some(trace), None) => pipeline::analyze(trace, &args.reference, &config_path),
            (None, Some(evidence)) => {
                pipeline::analyze_evidence(evidence, &args.reference, &config_path)
            }
            // The command line enforces this; a `Cli` built in code may not.
            _ => Err(KernelError::Path {
                kind: "analysis input",
                path: PathBuf::new(),
                reason: "analyze takes exactly one of a trace and --evidence",
            }
            .into()),
        },
        Command::Basecall(args) => pipeline::basecall(&args.trace, &config_path),
        Command::Evidence(args) => pipeline::evidence(&args.trace, &config_path),
        Command::Sample(args) => pipeline::sample(
            &args.sample_id,
            sample_input(&args.traces, &args.evidence)?,
            &args.reference,
            &config_path,
        ),
        Command::Consensus(args) => pipeline::consensus(
            &args.sample_id,
            sample_input(&args.traces, &args.evidence)?,
            &args.reference,
            args.run_context.as_deref(),
            &config_path,
        ),
        Command::Review(args) => pipeline::review(
            &args.sample_id,
            sample_input(&args.traces, &args.evidence)?,
            &args.reference,
            args.signal.as_deref(),
            args.run_context.as_deref(),
            &config_path,
        ),
        Command::RunContext(args) => pipeline::run_context(
            &args.run_id,
            (args.traces.as_deref(), args.evidence.as_deref()),
            &args.reference,
            &config_path,
        ),
        Command::Signal(args) => pipeline::signal(&args.run_id, &args.traces, &config_path),
        Command::Call(args) => pipeline::call(
            &args.sample_id,
            &args.sequences,
            &args.reference,
            &config_path,
        ),
        Command::Notation(args) => pipeline::notation(
            &args.sample_id,
            &args.variants,
            &args.reference,
            &config_path,
        ),
    }
}

/// A sample command's reads: its traces or their read evidence documents,
/// exactly one of the two.
fn sample_input<'a>(
    traces: &'a [PathBuf],
    evidence: &'a [PathBuf],
) -> Result<pipeline::SampleInput<'a>> {
    match (traces.is_empty(), evidence.is_empty()) {
        (false, true) => Ok(pipeline::SampleInput::Traces(traces)),
        (true, false) => Ok(pipeline::SampleInput::Evidence(evidence)),
        // The command line enforces this; a `Cli` built in code may not.
        _ => Err(KernelError::Path {
            kind: "sample input",
            path: PathBuf::new(),
            reason: "a sample command takes traces or --evidence documents, not both",
        }
        .into()),
    }
}
