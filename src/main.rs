#![forbid(unsafe_code)]
#![deny(deprecated)]

//! Operating-system boundary for the DNA binary.
//!
//! This file remains single-purpose: parse one command, call the library, print
//! a concise error, and return an exit code. It never scans `data/`, parses
//! `.env`, loads scientific configuration, or runs algorithms itself.

use std::process::ExitCode;

use clap::Parser;
use dna::cli::Cli;

fn main() -> ExitCode {
    match dna::run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}
