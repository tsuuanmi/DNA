//! Command-line argument definitions for focused DNA operations.
//!
//! Single-read commands accept exactly one AB1 path. The sample command accepts
//! one user-supplied sample identifier plus one or more AB1 paths and derives one
//! sample-evidence result. The call command accepts one sample identifier plus
//! one or more consensus-sequence FASTA files and runs the core caller alone.
//! The notation command accepts one sample identifier plus that sample's
//! variants document and derives its notation with the post-calling plugins.
//! Configuration is selected by `DNA_CONFIG`.

use std::path::PathBuf;

use clap::{Args, Parser, Subcommand};

/// DNA command-line interface.
#[derive(Debug, Parser)]
#[command(
    name = "DNA",
    version,
    about = "Call variants from Sanger traces or reviewed consensus sequences"
)]
pub struct Cli {
    /// Operation to run.
    #[command(subcommand)]
    pub command: Command,
}

/// Supported top-level commands.
#[derive(Debug, Subcommand)]
#[non_exhaustive]
pub enum Command {
    /// Analyze one AB1 trace against a reference FASTA.
    Analyze(AnalyzeArgs),
    /// Re-call and quality-trim one AB1 trace without a reference.
    Basecall(BasecallArgs),
    /// Aggregate independently analyzed AB1 traces into sample evidence.
    Sample(SampleArgs),
    /// Call variants from reviewed consensus sequences with the core alone.
    Call(CallArgs),
    /// Derive notation and conformance findings from a variants document.
    Notation(NotationArgs),
}

/// Arguments for the end-to-end reference analysis pipeline.
#[derive(Debug, Args)]
pub struct AnalyzeArgs {
    /// Input ABIF/AB1 trace.
    pub trace: PathBuf,

    /// Single-contig reference FASTA.
    #[arg(long)]
    pub reference: PathBuf,
}

/// Arguments for reference-free base re-calling and trimming.
#[derive(Debug, Args)]
pub struct BasecallArgs {
    /// Input ABIF/AB1 trace.
    pub trace: PathBuf,
}

/// Arguments for reference-coordinate sample evidence aggregation.
#[derive(Debug, Args)]
pub struct SampleArgs {
    /// Stable sample identifier used only for result/log naming and provenance.
    pub sample_id: String,

    /// Independently processed ABIF/AB1 traces belonging to the sample.
    #[arg(required = true, num_args = 1..)]
    pub traces: Vec<PathBuf>,

    /// Single-contig reference FASTA shared by every trace.
    #[arg(long)]
    pub reference: PathBuf,
}

/// Arguments for the core-only call over reviewed consensus sequences.
#[derive(Debug, Args)]
pub struct CallArgs {
    /// Stable sample identifier used only for result/log naming and provenance.
    pub sample_id: String,

    /// FASTA files whose records are reviewed consensus sequences of the sample.
    #[arg(required = true, num_args = 1..)]
    pub sequences: Vec<PathBuf>,

    /// Single-contig reference FASTA shared by every sequence.
    #[arg(long)]
    pub reference: PathBuf,
}

/// Arguments for post-calling notation of a variants document.
#[derive(Debug, Args)]
pub struct NotationArgs {
    /// Sample identifier; it must match the document's.
    pub sample_id: String,

    /// The sample's `dna.variants/v1` document.
    pub variants: PathBuf,

    /// Single-contig reference FASTA the variants were called against.
    #[arg(long)]
    pub reference: PathBuf,
}
