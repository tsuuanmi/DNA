#![forbid(unsafe_code)]
#![deny(deprecated)]

//! Library boundary for DNA operations.
//!
//! The source-module graph is routed from `src/README.md`. `lib.rs` remains
//! the minimal dispatcher: it exposes stable CLI and error boundaries, translates
//! parsed command arguments into application inputs, and keeps configuration,
//! scientific capabilities, orchestration, and reporting behind explicit module
//! boundaries.

mod alignment;
mod basecalling;
mod checksum;
pub mod cli;
pub mod config;
pub mod error;
mod input;
mod locus;
mod logger;
pub mod model;
mod pipeline;
mod quality_control;
mod read_processing;
mod reference;
mod report;
mod sample;
mod signal_processing;
pub mod variant_analysis;
mod variant_calling;

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
pub fn run(cli: Cli) -> Result<()> {
    let config_path = config::resolve_path();
    match cli.command {
        Command::Analyze(args) => pipeline::analyze(&args.trace, &args.reference, &config_path),
        Command::Basecall(args) => pipeline::basecall(&args.trace, &config_path),
        Command::Sample(args) => {
            pipeline::sample(&args.sample_id, &args.traces, &args.reference, &config_path)
        }
    }
}
