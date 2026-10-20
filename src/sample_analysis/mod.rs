//! Public Sample Analysis capability.
//!
//! One placement pass over a sample's traces yields both the sample evidence
//! and the consensus, without CLI logging, JSON serialization, or result
//! publication. The `sample` and `consensus` operations run the same pass,
//! from the traces or from their recorded read evidence (PROP-0011).
//! [`SampleAnalysis::call`] then calls the consensus as `call` would.

pub(crate) mod reads;

use std::collections::BTreeMap;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

use tracing::{Dispatch, Span};

use crate::config::Config;
use crate::error::{Result, SampleError};
use crate::input::sanger::{self, SampleInputs};
use crate::model::sanger_read::{Rejection, SangerRead};
use dna_call::model::called_read::CalledRead;
use dna_call::read_trust;
use dna_kernel::model::reference::Reference;
pub use dna_kernel::profile::{Convention, ConventionKind, ConventionLabel};
use dna_kernel::profile::{Profile, ProfileIdentity};
use dna_kernel::read_evidence::{EvidenceReason, Placement};
pub use dna_kernel::read_evidence::{LengthMarginal, Margin, TractLengths};
use dna_kernel::variant::ReferenceIdentity;
use dna_sample::model::sample_evidence::{LocusObservations, RejectedSampleRead};

pub use crate::model::read_observation::{RejectedRead, RejectionCause};
use dna_align::model::alignment::Orientation;
use dna_call::components::component_calls;
pub use dna_call::components::{ComponentCall, MixtureCalls};
pub use dna_call::model::trust::{CallTrust, TrustScope, Withholder};
pub use dna_consensus::model::consensus::Consensus;
pub use dna_empop::export::{Cut, CutReason, Export as EmpopExport};
pub use dna_sample::model::sample_evidence::SampleEvidence;
pub use dna_sanger::model::attachment::{SangerAttachment, SangerRejection};
pub use dna_sanger::run_signal::{Contributions, RunSignal};

pub use crate::variant_analysis::ReferenceSegment;

use reads::{PlacedReads, ReadScope, Reads};

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
    /// Sanger products of each admitted read, by read content SHA-256: of
    /// every admitted read, since this API analyzes traces.
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
    /// Each admitted read's Sanger evidence as the documents report it, by
    /// read content SHA-256: what every report reads. Unlike [`Self::reads`],
    /// it is complete for an analysis the CLI runs from recorded read
    /// evidence (PROP-0011), which keeps no trace products.
    pub(crate) sanger: BTreeMap<String, SangerRead>,
    /// Each rejected read as the documents report it, by read content
    /// SHA-256; complete as [`Self::sanger`] is.
    pub(crate) rejections: BTreeMap<String, Rejection>,
}

impl SampleAnalysis {
    /// An empty run signal of `run_id`, made with the flank and resolution of
    /// this analysis's configuration, for its reads' contributions to join.
    #[must_use]
    pub fn run_signal(&self, run_id: String) -> RunSignal {
        RunSignal::build(
            run_id,
            std::collections::BTreeSet::new(),
            [],
            &self.loaded.config.sanger.signal_processing.run_signal,
        )
    }

    /// What the admitted reads contribute to their run's signal
    /// (SRS-CALL-019), summed over the sample, with the content SHA-256 of
    /// those reads: what a run signal holds of this sample.
    #[must_use]
    pub fn signal_contributions(&self) -> (Vec<String>, Contributions) {
        let mut traces = Vec::new();
        let mut sum = Contributions::default();
        for read in &self.called {
            let Some(sanger) = self.sanger.get(&read.input_sha256) else {
                continue;
            };
            let own = dna_sanger::run_signal::contributions(
                &sanger.heights(),
                &read.evidence,
                &self.loaded.config.sanger.signal_processing.run_signal,
            );
            traces.push(read.input_sha256.clone());
            sum.calls += own.calls;
            for (key, shares) in own.keys {
                let into = sum.keys.entry(key).or_default();
                for (share, count) in shares {
                    *into.entry(share).or_default() += count;
                }
            }
        }
        (traces, sum)
    }

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

    /// Every admitted read's phase mixtures (ADR-0082) with each component
    /// placed on the reference and called on its own, by read: its name, its
    /// content SHA-256, whether it reads the reference strand, and its
    /// mixtures. Evidence for review and evaluation; no call changes, and a
    /// component that cannot be placed is left out, never failing the sample.
    #[must_use]
    pub fn components(&self) -> Vec<ReadMixtures> {
        let config = &self.loaded.config.core;
        self.called
            .iter()
            .map(|read| ReadMixtures {
                read: read.input_name.clone(),
                sha256: read.input_sha256.clone(),
                forward: read.alignment.orientation() == Orientation::Forward,
                mixtures: component_calls(
                    read,
                    &self.loaded.reference,
                    &config.alignment,
                    &config.variant_calling,
                    &self.loaded.profile.regions,
                ),
            })
            .collect()
    }

    /// Every read placed on the reference and called on its own, by read,
    /// with its tract decodes (PROP-0012): its name, its content SHA-256,
    /// whether it reads the reference strand, and each tract with the
    /// reference position of its first call, the read's reason for masking
    /// its onset call, and whether the read stays out of phase to its end. Evidence for evaluation; the
    /// core reads none of it.
    #[must_use]
    pub fn tracts(&self) -> Vec<ReadTracts> {
        self.called
            .iter()
            .map(|read| {
                let positions = read
                    .alignment
                    .columns()
                    .iter()
                    .filter_map(|column| {
                        Some((
                            column.original_call_index_0based()?,
                            column.reference_index_0based()? + 1,
                        ))
                    })
                    .collect::<BTreeMap<_, _>>();
                let calls = read.evidence.calls();
                ReadTracts {
                    read: read.input_name.clone(),
                    sha256: read.input_sha256.clone(),
                    forward: read.alignment.orientation() == Orientation::Forward,
                    tracts: read
                        .evidence
                        .tracts()
                        .iter()
                        .map(|tract| PlacedTract {
                            reference: positions.get(&tract.call_start).copied(),
                            onset_reason: calls
                                .get(tract.onset)
                                .and_then(|call| call.placement.reason())
                                .map(EvidenceReason::label),
                            to_end: calls
                                .iter()
                                .skip(tract.onset)
                                .all(|call| call.placement != Placement::Scored),
                            tract: tract.clone(),
                        })
                        .collect(),
                }
            })
            .collect()
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
        Ok(dna_sample::observations(
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

/// One admitted read's phase mixtures, each component compared with the
/// reference (ADR-0082).
#[non_exhaustive]
#[derive(Debug, Clone)]
pub struct ReadMixtures {
    /// The read's name.
    pub read: String,
    /// Its content SHA-256.
    pub sha256: String,
    /// It reads the reference strand.
    pub forward: bool,
    /// Its mixtures.
    pub mixtures: Vec<MixtureCalls>,
}

/// One read's tract decodes (PROP-0012).
#[derive(Debug, Clone, PartialEq)]
pub struct ReadTracts {
    /// The read's name.
    pub read: String,
    /// Its content SHA-256.
    pub sha256: String,
    /// It reads the reference strand.
    pub forward: bool,
    /// Its tracts, in onset order.
    pub tracts: Vec<PlacedTract>,
}

/// One tract decode, with where its first call lies on the reference.
#[derive(Debug, Clone, PartialEq)]
pub struct PlacedTract {
    /// The reference position (1-based) the tract's first call is aligned
    /// to; `None` when the call is not aligned to a reference base.
    pub reference: Option<usize>,
    /// The modality's reason for masking the onset call, the read's own
    /// assertion of why it left phase there (`post_homopolymer`,
    /// `dephased_signal`, `weak_signal`, ...); `None` when the call is not
    /// masked.
    pub onset_reason: Option<&'static str>,
    /// The read does not regain phase after the onset: no later call is
    /// scored.
    pub to_end: bool,
    /// The modality's assertion.
    pub tract: TractLengths,
}

/// One item of a sample's review list that a reviewer must look at: a
/// `check` item of `dna review` (SRS-API-015).
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct ReviewCheck {
    /// First reference position of the item (1-based).
    pub first: usize,
    /// Last reference position of the item (1-based, inclusive).
    pub last: usize,
    /// The item's label names, sorted.
    pub labels: Vec<&'static str>,
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
    let SampleInputs {
        config,
        profile,
        traces,
        reference,
    } = sanger::load_sample(traces, reference, config)?;
    let loaded = Loaded {
        config,
        profile,
        reference,
    };
    Ok(
        match analyze(
            Reads::Traces(&traces),
            loaded,
            &CallerScope::current(),
            workers,
        )? {
            Outcome::Analyzed(analysis) => SampleOutcome::Analyzed(analysis),
            Outcome::Rejected { products, .. } => SampleOutcome::Rejected(products),
        },
    )
}

/// What one sample analysis yields inside the crate: the analysis, or, when
/// every read was rejected, those reads as the documents report them and,
/// for reads from their traces, with their Sanger products.
pub(crate) enum Outcome {
    /// At least one read was admitted.
    Analyzed(Box<SampleAnalysis>),
    /// Every read was rejected.
    Rejected {
        /// Each read as the documents report it, by read content SHA-256.
        rejections: BTreeMap<String, Rejection>,
        /// Each read with its Sanger products, by read content SHA-256;
        /// empty for reads from their recorded evidence.
        products: BTreeMap<String, RejectedRead>,
    },
}

/// Runs the placement pass, the aggregation, and the consensus over a
/// sample's reads, from their traces or their recorded evidence, with each
/// read's instrumentation going to `scope`; the shared body of
/// [`analyze_sanger_outcome`] and the `review` operation.
pub(crate) fn analyze(
    reads: Reads<'_>,
    loaded: Loaded,
    scope: &impl ReadScope,
    workers: NonZeroUsize,
) -> Result<Outcome> {
    let placed = reads::place(
        reads,
        &loaded.reference,
        &loaded.config,
        &loaded.profile,
        scope,
        workers,
    )?;
    if placed.called.is_empty() && !placed.rejected.is_empty() {
        return Ok(Outcome::Rejected {
            rejections: by_identity(placed.rejected),
            products: rejected_by_identity(placed.products.rejected),
        });
    }
    let evidence = aggregate(&placed, &loaded.config)?;
    let consensus = build_consensus(&placed, &evidence, &loaded.reference, &loaded.config);
    Ok(Outcome::Analyzed(Box::new(SampleAnalysis {
        reference: loaded.reference.identity(),
        profile: loaded.profile.identity().clone(),
        evidence,
        consensus,
        reads: placed.products.reads,
        rejected: rejected_by_identity(placed.products.rejected),
        loaded,
        called: placed.called,
        sanger: placed.sanger,
        rejections: by_identity(placed.rejected),
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
    Ok(dna_sample::aggregate(
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
    dna_consensus::build(&ordered, reference, &config.core, &config.consensus)
}

/// Rejected reads keyed by read content identity.
pub(crate) fn by_identity(rejected: Vec<Rejection>) -> BTreeMap<String, Rejection> {
    rejected
        .into_iter()
        .map(|read| (read.input_sha256.clone(), read))
        .collect()
}

/// Rejected reads with their Sanger products keyed by read content identity.
fn rejected_by_identity(rejected: Vec<RejectedRead>) -> BTreeMap<String, RejectedRead> {
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
