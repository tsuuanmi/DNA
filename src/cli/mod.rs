//! Command-line argument definitions for focused DNA operations.
//!
//! Single-read commands accept exactly one AB1 path; analyze accepts the
//! read's evidence document (`--evidence`) in its place. The sample command accepts
//! one user-supplied sample identifier plus one or more AB1 paths, or their
//! read evidence documents (`--evidence`), and derives one sample-evidence
//! result; consensus and review take the same reads. The call command accepts one sample identifier plus
//! one or more consensus-sequence FASTA files and runs the core caller alone.
//! The notation command accepts one sample identifier plus that sample's
//! variants document and derives its notation with the post-calling plugins.
//! The review command takes the sample command's arguments and lists the
//! positions a reviewer should look at.
//! Configuration is selected by `DNA_CONFIG`.

use std::path::PathBuf;

use clap::{ArgGroup, Args, Parser, Subcommand};

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
    /// Analyze one AB1 trace, or one read's recorded evidence, against a reference FASTA.
    Analyze(AnalyzeArgs),
    /// Re-call and quality-trim one AB1 trace without a reference.
    Basecall(BasecallArgs),
    /// Write one AB1 trace's per-call evidence as the core receives it, without a reference.
    Evidence(BasecallArgs),
    /// Aggregate independently analyzed AB1 traces, or their recorded evidence, into sample evidence.
    Sample(SampleArgs),
    /// Adjudicate a sample's AB1 traces, or their recorded evidence, into a consensus sequence.
    Consensus(ConsensusArgs),
    /// List the positions of a sample a reviewer should look at, with their evidence.
    Review(ReviewArgs),
    /// Call variants from reviewed consensus sequences with the core alone.
    Call(CallArgs),
    /// Record the second-channel signal a run's AB1 traces share.
    Signal(SignalArgs),
    /// Record what a laboratory batch's samples give each of them: the run context.
    RunContext(RunContextArgs),
    /// Derive notation and conformance findings from a variants document.
    Notation(NotationArgs),
}

/// Arguments for the end-to-end reference analysis pipeline: a trace, or the
/// read's `dna.read_evidence/v1` document, never both.
#[derive(Debug, Args)]
#[command(group(ArgGroup::new("read").required(true).args(["trace", "evidence"])))]
pub struct AnalyzeArgs {
    /// Input ABIF/AB1 trace.
    pub trace: Option<PathBuf>,

    /// The read's evidence (`dna.read_evidence/v1`, written by `evidence`), in
    /// place of its trace (PROP-0011).
    #[arg(long)]
    pub evidence: Option<PathBuf>,

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

/// Arguments for reference-coordinate sample evidence aggregation: the
/// sample's traces, or their `dna.read_evidence/v1` documents, never both.
#[derive(Debug, Args)]
#[command(group(ArgGroup::new("reads").required(true).args(["traces", "evidence"])))]
pub struct SampleArgs {
    /// Stable sample identifier used only for result/log naming and provenance.
    pub sample_id: String,

    /// Independently processed ABIF/AB1 traces belonging to the sample.
    pub traces: Vec<PathBuf>,

    /// The reads' evidence (`dna.read_evidence/v1`, written by `evidence`),
    /// in place of their traces (PROP-0011).
    #[arg(long, num_args = 1..)]
    pub evidence: Vec<PathBuf>,

    /// Single-contig reference FASTA shared by every trace.
    #[arg(long)]
    pub reference: PathBuf,
}

/// Arguments for a sample consensus: the sample's traces, or their
/// `dna.read_evidence/v1` documents, never both.
#[derive(Debug, Args)]
#[command(group(ArgGroup::new("reads").required(true).args(["traces", "evidence"])))]
pub struct ConsensusArgs {
    /// Stable sample identifier used only for result/log naming and provenance.
    pub sample_id: String,

    /// Independently processed ABIF/AB1 traces belonging to the sample.
    pub traces: Vec<PathBuf>,

    /// The reads' evidence (`dna.read_evidence/v1`, written by `evidence`),
    /// in place of their traces (PROP-0011).
    #[arg(long, num_args = 1..)]
    pub evidence: Vec<PathBuf>,

    /// Single-contig reference FASTA shared by every trace.
    #[arg(long)]
    pub reference: PathBuf,

    /// The run context (`dna.run_context/v1`, written by `run-context`) of
    /// the sample's batch, to calibrate the length rule by (PROP-0014).
    #[arg(long)]
    pub run_context: Option<PathBuf>,
}

/// Arguments for a run context over one laboratory batch: a tab-separated
/// list of `sample_id` and read path, one read per line, of traces or of
/// their `dna.read_evidence/v1` documents, never both.
#[derive(Debug, Args)]
#[command(group(ArgGroup::new("reads").required(true).args(["traces", "evidence"])))]
pub struct RunContextArgs {
    /// Stable batch identifier used for result/log naming and provenance.
    pub run_id: String,

    /// The batch's ABIF/AB1 traces, as `sample_id<TAB>path` lines.
    #[arg(long)]
    pub traces: Option<PathBuf>,

    /// The batch's read evidence documents, as `sample_id<TAB>path` lines.
    #[arg(long)]
    pub evidence: Option<PathBuf>,

    /// Single-contig reference FASTA shared by every read.
    #[arg(long)]
    pub reference: PathBuf,
}

/// Arguments for a sample review: the sample's traces, or their
/// `dna.read_evidence/v1` documents, never both.
#[derive(Debug, Args)]
#[command(group(ArgGroup::new("reads").required(true).args(["traces", "evidence"])))]
pub struct ReviewArgs {
    /// Stable sample identifier used only for result/log naming and provenance.
    pub sample_id: String,

    /// Independently processed ABIF/AB1 traces belonging to the sample.
    pub traces: Vec<PathBuf>,

    /// The reads' evidence (`dna.read_evidence/v1`, written by `evidence`),
    /// in place of their traces (PROP-0011).
    #[arg(long, num_args = 1..)]
    pub evidence: Vec<PathBuf>,

    /// Single-contig reference FASTA shared by every trace.
    #[arg(long)]
    pub reference: PathBuf,

    /// The run signal (`dna.signal/v1`) of the sample's run, to judge second
    /// bases against (PROP-0009).
    #[arg(long)]
    pub signal: Option<PathBuf>,

    /// The run context (`dna.run_context/v1`, written by `run-context`) of
    /// the sample's batch, to calibrate the length rule by (PROP-0014).
    #[arg(long)]
    pub run_context: Option<PathBuf>,
}

/// Arguments for a run signal over one run's traces.
#[derive(Debug, Args)]
pub struct SignalArgs {
    /// Stable run identifier used for result/log naming and provenance.
    pub run_id: String,

    /// Every ABIF/AB1 trace of the run.
    #[arg(required = true, num_args = 1..)]
    pub traces: Vec<PathBuf>,
}

/// Arguments for the core-only call over reviewed consensus sequences.
#[derive(Debug, Args)]
pub struct CallArgs {
    /// Stable sample identifier used only for result/log naming and provenance.
    pub sample_id: String,

    /// FASTA files whose records are reviewed consensus sequences of the sample,
    /// or the sample's `dna.consensus/v1` document alone (its decided edits).
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
