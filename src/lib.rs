//! Library boundary for DNA operations.
//!
//! The source-module graph is routed from `src/README.md`. `lib.rs` remains
//! the minimal dispatcher: it exposes stable CLI and error boundaries, translates
//! parsed command arguments into application inputs, and keeps configuration,
//! scientific capabilities, orchestration, and reporting behind explicit module
//! boundaries.

mod alignment;
mod basecalling;
mod callability;
mod checksum;
pub mod cli;
mod config;
mod conformance;
pub mod error;
mod input;
mod locus;
mod model;
mod operation_log;
mod pipeline;
mod plugin;
pub mod profile;
mod quality_control;
mod read_evidence;
mod read_processing;
mod reference;
mod report;
mod sample;
mod signal_processing;
pub mod variant;
pub mod variant_analysis;
mod variant_calling;
pub mod variant_nomenclature;
pub mod variant_normalization;
mod variant_representation;

#[cfg(feature = "fuzzing")]
#[doc(hidden)]
pub mod fuzzing {
    /// Exercises the bounds-checked ABIF directory parser without filesystem I/O.
    pub fn parse_abif(bytes: &[u8]) {
        let _ = crate::input::sanger::abif::parse_container(bytes.to_vec());
    }
}

use cli::{Cli, Command};
use error::Result;

/// Dispatches a parsed command through the application boundary.
///
/// # Errors
///
/// Returns [`Error`](error::Error) when the operation fails; the failure is also
/// recorded in the operation log unless the log itself cannot be written.
pub fn run(cli: Cli) -> Result<()> {
    let config_path = config::resolve_path();
    match cli.command {
        Command::Analyze(args) => pipeline::analyze(&args.trace, &args.reference, &config_path),
        Command::Basecall(args) => pipeline::basecall(&args.trace, &config_path),
        Command::Sample(args) => {
            pipeline::sample(&args.sample_id, &args.traces, &args.reference, &config_path)
        }
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
