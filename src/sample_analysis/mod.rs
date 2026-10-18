//! Public Sample Analysis capability.
//!
//! One placement pass over a sample's traces yields both the sample evidence
//! and the consensus, without CLI logging, JSON serialization, or result
//! publication. The `sample` and `consensus` operations run the same pass.
//! [`SampleAnalysis::call`] then calls the consensus as `call` would.

pub(crate) mod reads;

use std::collections::BTreeMap;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

use tracing::{Dispatch, Span};

use crate::config::Config;
use crate::error::{Result, SampleError};
use crate::input::sanger::{self, SampleInputs};
use dna_core::consensus;
use dna_core::model::called_read::CalledRead;
use dna_core::model::sample_evidence::{LocusObservations, RejectedSampleRead};
use dna_core::sample;
use dna_core::variant_calling::read_trust;
use dna_kernel::model::reference::Reference;
pub use dna_kernel::profile::{Convention, ConventionKind, ConventionLabel};
use dna_kernel::profile::{Profile, ProfileIdentity};
use dna_kernel::variant::ReferenceIdentity;

pub use crate::model::read_observation::{RejectedRead, RejectionCause};
pub use dna_core::model::consensus::Consensus;
pub use dna_core::model::sample_evidence::SampleEvidence;
pub use dna_core::model::trust::{CallTrust, TrustScope, Withholder};
pub use dna_empop::export::{Cut, CutReason, Export as EmpopExport};
pub use dna_sanger::model::attachment::{SangerAttachment, SangerRejection};

pub use crate::variant_analysis::ReferenceSegment;

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
    /// Each rejected read with its cause, by read content SHA-256; the cause
    /// carries the read's Sanger evidence unless the read was too short to
    /// process.
    pub rejected: BTreeMap<String, RejectedRead>,
    /// The configuration, profile, and reference the analysis ran with; the
    /// identities above derive from them, and [`SampleAnalysis::call`] calls
    /// the consensus with them.
    pub(crate) loaded: Loaded,
    /// The admitted reads as the core placed them, which
    /// [`SampleAnalysis::observations`] reads.
    pub(crate) called: Vec<CalledRead>,
}

impl SampleAnalysis {
    /// The core's trust in every call of the admitted read whose content
    /// SHA-256 is `input_sha256`, as its variant calling and consensus judge
    /// it; `None` for another read. Diagnostic evidence for evaluation and
    /// review, not a calling contract.
    #[must_use]
    pub fn trust(&self, input_sha256: &str) -> Option<Vec<CallTrust>> {
        self.called
            .iter()
            .find(|read| read.input_sha256 == input_sha256)
            .map(|read| read_trust(&read.evidence, &self.loaded.config.core.variant_calling))
    }

    /// The reportable regions of the target profile the analysis ran with:
    /// inclusive 1-based reference intervals, as the profile declares them.
    #[must_use]
    pub fn regions(&self) -> &[[usize; 2]] {
        &self.loaded.profile.regions
    }

    /// Every admitted read's observation of each reference position from
    /// `first_1based` to `last_1based` (inclusive) that a read covers,
    /// whether or not a read differs there: state, reference-strand base,
    /// source call index, and reference-oriented A/C/G/T profile, as
    /// `evidence.locus_differences` records them, with read indices into
    /// `evidence.reads`, but without their support aggregates. Uncovered
    /// positions have no record. The interval does not wrap: on a circular
    /// reference, a span across the origin is two calls.
    ///
    /// # Errors
    ///
    /// Returns [`SampleError::InvalidInterval`]
    /// when `first_1based` is zero or greater than `last_1based`, and
    /// [`Error`](crate::error::Error) when the placed reads are
    /// inconsistent, which the analysis already rules out.
    pub fn observations(
        &self,
        first_1based: usize,
        last_1based: usize,
    ) -> Result<Vec<LocusObservations>> {
        Ok(sample::observations(
            &self.called.iter().collect::<Vec<_>>(),
            first_1based,
            last_1based,
        )?)
    }
}

/// The configuration, profile, and reference an analysis ran with.
#[derive(Debug, Clone)]
pub(crate) struct Loaded {
    pub(crate) config: Config,
    pub(crate) profile: Profile,
    pub(crate) reference: Reference,
}

/// Typed result of calling a sample's consensus (SRS-API-012).
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct ConsensusCall {
    /// Reference intervals the consensus sequences cover with callable bases,
    /// in consensus read order; `callable_reference_segments` of the
    /// `dna.variants/v1` reads.
    pub callable_segments: Vec<ReferenceSegment>,
    /// The sample's calls in the profile notation, distinct and in reference
    /// order, as `notation.calls` of `dna.variants/v1`; `None` when the profile
    /// declares no notation.
    pub notation: Option<Vec<NotationCall>>,
    /// The EMPOP query `dna notation` would export from the `dna.variants/v1`
    /// document (PROP-0007); `None` when the profile declares no `[empop]`.
    pub empop: Option<EmpopExport>,
}

/// One per-base call in the profile notation.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct NotationCall {
    /// 1-based position of the changed base, or of the base an insertion
    /// follows.
    pub position: usize,
    /// `0` for a substitution or deletion, `1..` for successive inserted bases.
    pub ordinal: usize,
    /// The reference base; `None` for an inserted base.
    pub reference: Option<char>,
    /// The called base; `None` for a deleted base.
    pub alternate: Option<char>,
    /// The rendered call, for example `309.1C`, `249DEL`, or `16189C`.
    pub text: String,
    /// The profile's report convention labelling the call (SRS-NOM-017),
    /// when one does; the call is kept either way.
    pub convention: Option<Convention>,
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
/// Returns [`Error`](crate::error::Error) when an input or configuration
/// file is unreadable or invalid, when a read fails other than by rejection
/// (the first such read in trace order), or when no read is admitted
/// (`SampleError::NoAdmittedReads`; [`analyze_sanger_outcome`] returns those
/// reads instead).
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
    match analyze_sanger_outcome(traces, reference, config, workers)? {
        SampleOutcome::Analyzed(analysis) => Ok(*analysis),
        SampleOutcome::Rejected(rejected) => Err(SampleError::NoAdmittedReads {
            rejected: rejected.len(),
        }
        .into()),
    }
}

/// What one sample analysis yields: the analysis, or, when every read was
/// rejected, those reads. The two cases are complete: a read is admitted or
/// rejected.
#[derive(Debug, Clone)]
pub enum SampleOutcome {
    /// At least one read was admitted.
    Analyzed(Box<SampleAnalysis>),
    /// Every read was rejected; each with its cause, by read content SHA-256.
    Rejected(BTreeMap<String, RejectedRead>),
}

/// Runs [`analyze_sanger`], but returns the rejected reads with their causes
/// when no read is admitted, so a caller can name the traces that failed.
///
/// # Errors
///
/// As [`analyze_sanger`], except that a sample without an admitted read is
/// [`SampleOutcome::Rejected`] rather than an error.
///
/// # Panics
///
/// As [`analyze_sanger`].
pub fn analyze_sanger_outcome(
    traces: &[PathBuf],
    reference: &Path,
    config: &Path,
    workers: NonZeroUsize,
) -> Result<SampleOutcome> {
    analyze(
        sanger::load_sample(traces, reference, config)?,
        &CallerScope::current(),
        workers,
    )
}

/// Runs the placement pass, the aggregation, and the consensus over loaded
/// inputs, with each read's instrumentation going to `scope`; the shared body
/// of [`analyze_sanger_outcome`] and the `review` operation.
pub(crate) fn analyze(
    inputs: SampleInputs,
    scope: &impl ReadScope,
    workers: NonZeroUsize,
) -> Result<SampleOutcome> {
    let placed = reads::place(
        &inputs.traces,
        &inputs.reference,
        &inputs.config,
        &inputs.profile,
        scope,
        workers,
    )?;
    if placed.called.is_empty() && !placed.rejected.is_empty() {
        return Ok(SampleOutcome::Rejected(by_identity(placed.rejected)));
    }
    let evidence = aggregate(&placed, &inputs.config)?;
    let consensus = build_consensus(&placed, &evidence, &inputs.reference, &inputs.config);
    let loaded = Loaded {
        config: inputs.config,
        profile: inputs.profile,
        reference: inputs.reference,
    };
    Ok(SampleOutcome::Analyzed(Box::new(SampleAnalysis {
        reference: loaded.reference.identity(),
        profile: loaded.profile.identity().clone(),
        evidence,
        consensus,
        reads: placed.sanger,
        rejected: by_identity(placed.rejected),
        loaded,
        called: placed.called,
    })))
}

/// Aggregates the placed reads into the sample evidence, which validates the
/// read set and fixes the read order.
pub(crate) fn aggregate(placed: &PlacedReads, config: &Config) -> Result<SampleEvidence> {
    let rejected = placed
        .rejected
        .iter()
        .map(|read| RejectedSampleRead {
            input_name: read.input_name.clone(),
            input_sha256: read.input_sha256.clone(),
        })
        .collect::<Vec<_>>();
    Ok(sample::aggregate(
        &placed.called.iter().collect::<Vec<_>>(),
        &rejected,
        &config.sample,
        &config.core.variant_calling,
    )?)
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
        .map(|read| (read.input_sha256.clone(), read))
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
