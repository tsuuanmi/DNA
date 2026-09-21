//! Command-line argument definitions for focused DNA operations.
//!
//! Single-read commands accept exactly one AB1 path. The sample command accepts
//! one user-supplied sample identifier plus one or more AB1 paths and derives one
//! sample-evidence result. Configuration is selected by `DNA_CONFIG`.

use std::path::PathBuf;

use clap::{Args, Parser, Subcommand};

/// DNA command-line interface.
#[derive(Debug, Parser)]
#[command(name = "DNA", version, about = "Process Sanger sequencing traces")]
pub struct Cli {
    /// Operation to run.
    #[command(subcommand)]
    pub command: Command,
}

/// Supported top-level commands.
#[derive(Debug, Subcommand)]
pub enum Command {
    /// Analyze one AB1 trace against a reference FASTA.
    Analyze(AnalyzeArgs),
    /// Re-call and quality-trim one AB1 trace without a reference.
    Basecall(BasecallArgs),
    /// Aggregate independently analyzed AB1 traces into sample evidence.
    Sample(SampleArgs),
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
