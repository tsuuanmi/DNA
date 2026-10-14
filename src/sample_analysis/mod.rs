//! Public Sample Analysis capability.
//!
//! One placement pass over a sample's traces yields both the sample evidence
//! and the consensus, without CLI logging, JSON serialization, or result
//! publication. The `sample` and `consensus` operations run the same pass.

pub(crate) mod reads;

use std::collections::BTreeMap;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

use tracing::{Dispatch, Span};

use crate::config::Config;
use crate::input::sanger;
use dna_core::consensus;
use dna_core::model::sample_evidence::RejectedSampleRead;
use dna_core::sample;
use dna_kernel::error::Result;
use dna_kernel::model::reference::Reference;
use dna_kernel::profile::ProfileIdentity;
use dna_kernel::variant::ReferenceIdentity;

pub use crate::model::read_observation::{RejectedRead, RejectionCause};
pub use dna_core::model::consensus::Consensus;
pub use dna_core::model::sample_evidence::SampleEvidence;
pub use dna_sanger::model::attachment::{SangerAttachment, SangerRejection};

use reads::{PlacedReads, ReadScope};

/// Typed result of one sample analysis: the sample evidence and the consensus
/// of one placement pass.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct SampleAnalysis {
    /// Identity of the reference every read was placed on.
    pub reference: ReferenceIdentity,
    /// Identity of the target profile the configuration references.
    pub profile: ProfileIdentity,
    /// Core sample evidence: admitted and rejected reads, coverage, overlaps,
    /// locus differences, and variants with their support. It carries the
    /// configuration identity.
    pub evidence: SampleEvidence,
    /// Consensus of the admitted reads.
    pub consensus: Consensus,
    /// Sanger evidence of each admitted read, by read content SHA-256.
    pub reads: BTreeMap<String, SangerAttachment>,
    /// Each rejected read with its Sanger evidence and cause, by read content
    /// SHA-256.
    pub rejected: BTreeMap<String, RejectedRead>,
}

/// Runs the Sanger adapter over one sample's traces: each trace is processed
/// and placed once, then the reads are aggregated into the sample evidence
/// and adjudicated into the consensus.
///
/// Reads are processed on up to `workers` scoped threads; with one worker
/// they run one after another on the calling thread. The result does not
/// depend on the number of workers or on scheduling. Each read's
/// instrumentation goes to the caller's current `tracing` dispatcher, inside
/// the caller's current span, so with more than one worker the events of
/// different reads may interleave.
///
/// # Errors
///
/// Returns [`Error`](dna_kernel::error::Error) when an input or configuration
/// file is unreadable or invalid, when a read fails other than by rejection
/// (the first such read in trace order), or when no read is admitted.
///
/// # Panics
///
/// Resumes a panic of a read's processing, a defect, when it is the first
/// failure in trace order.
pub fn analyze_sanger(
    traces: &[PathBuf],
    reference: &Path,
    config: &Path,
    workers: NonZeroUsize,
) -> Result<SampleAnalysis> {
    let inputs = sanger::load_sample(traces, reference, config)?;
    let placed = reads::place(
        &inputs.traces,
        &inputs.reference,
        &inputs.config,
        &inputs.profile,
        &CallerScope::current(),
        workers,
    )?;
    let evidence = aggregate(&placed, &inputs.config)?;
    let consensus = build_consensus(&placed, &evidence, &inputs.reference, &inputs.config);
    Ok(SampleAnalysis {
        reference: inputs.reference.identity(),
        profile: inputs.profile.identity().clone(),
        evidence,
        consensus,
        reads: placed.sanger,
        rejected: by_identity(placed.rejected),
    })
}

/// Aggregates the placed reads into the sample evidence, which validates the
/// read set and fixes the read order.
pub(crate) fn aggregate(placed: &PlacedReads, config: &Config) -> Result<SampleEvidence> {
    let rejected = placed
        .rejected
        .iter()
        .map(|read| RejectedSampleRead {
            input_name: read.sanger.input_name.clone(),
            input_sha256: read.sanger.input_sha256.clone(),
        })
        .collect::<Vec<_>>();
    sample::aggregate(
        &placed.called.iter().collect::<Vec<_>>(),
        &rejected,
        &config.core.sample_reconciliation,
    )
}

/// Builds the consensus of the admitted reads in the evidence's read order.
pub(crate) fn build_consensus(
    placed: &PlacedReads,
    evidence: &SampleEvidence,
    reference: &Reference,
    config: &Config,
) -> Consensus {
    let by_identity = placed
        .called
        .iter()
        .map(|read| (read.input_sha256.as_str(), read))
        .collect::<BTreeMap<_, _>>();
    let ordered = evidence
        .reads
        .iter()
        .filter_map(|read| by_identity.get(read.input_sha256.as_str()).copied())
        .collect::<Vec<_>>();
    consensus::build(&ordered, reference, &config.core)
}

/// Rejected reads keyed by read content identity.
pub(crate) fn by_identity(rejected: Vec<RejectedRead>) -> BTreeMap<String, RejectedRead> {
    rejected
        .into_iter()
        .map(|read| (read.sanger.input_sha256.clone(), read))
        .collect()
}

/// Sends each read's instrumentation to the dispatcher and span that were
/// current on the caller's thread, so a read's spans have the same parent on
/// every worker.
struct CallerScope {
    dispatch: Dispatch,
    span: Span,
}

impl CallerScope {
    fn current() -> Self {
        Self {
            dispatch: tracing::dispatcher::get_default(Dispatch::clone),
            span: Span::current(),
        }
    }
}

impl ReadScope for CallerScope {
    type Records = ();

    fn record<T>(&self, work: impl FnOnce() -> T) -> (T, ()) {
        let value = tracing::dispatcher::with_default(&self.dispatch, || self.span.in_scope(work));
        (value, ())
    }

    fn replay(&self, (): ()) {}
}
