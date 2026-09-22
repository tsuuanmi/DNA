#![forbid(unsafe_code)]
#![deny(deprecated)]

//! Library boundary for DNA operations.
//!
//! The source-module graph is routed from `src/README.md`. `lib.rs` remains
//! the minimal dispatcher: it exposes stable CLI and error boundaries, routes
//! commands, and keeps configuration, decoding, scientific stages, and reporting
//! behind the pipeline boundary.

mod alignment;
mod basecalling;
mod checksum;
pub mod cli;
pub mod config;
pub mod error;
mod locus;
mod logger;
pub mod model;
mod pipeline;
mod quality_control;
mod reference;
mod report;
mod sample;
mod signal_processing;
mod trace;
mod variant_calling;

use cli::{Cli, Command};
use error::Result;

/// Dispatches a parsed command through the application boundary.
pub fn run(cli: Cli) -> Result<()> {
    match cli.command {
        Command::Analyze(args) => pipeline::analyze(&args),
        Command::Basecall(args) => pipeline::basecall(&args),
        Command::Sample(args) => pipeline::sample(&args),
    }
}
