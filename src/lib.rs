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
mod report;
pub mod variant_analysis;

pub use dna_kernel::{error, variant};

/// Target profiles: knowledge about one sequencing target.
pub mod profile {
    pub use dna_kernel::profile::{Profile, ProfileIdentity};
}

/// Optional haplotype-preserving representation normalization of called
/// variants.
pub mod variant_normalization {
    pub use dna_post::variant_normalization::{
        NormalizationPolicy, VariantNormalizationResult, normalize,
    };
}

/// Optional profile-driven target nomenclature of normalized variants.
pub mod variant_nomenclature {
    pub use dna_post::variant_nomenclature::{
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
        Command::Consensus(args) => {
            pipeline::consensus(&args.sample_id, &args.traces, &args.reference, &config_path)
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
