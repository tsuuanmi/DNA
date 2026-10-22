//! The review list of one sample, `dna.review/v1` (PROP-0004).
//!
//! Items and labels are computed from the sample analysis and its call alone,
//! never from a truth: the consensus sites, the reads' observations and the
//! core's trust in each call, and the Sanger scan position and relative
//! quality of each call, joined by read identity and call index. The rules'
//! version is the `review` plugin's version.

use std::collections::{BTreeMap, BTreeSet};

use crate::config::ReviewConfig;
use crate::error::{ReportError, Result, RunSignalError};
use crate::model::consensus_result::StretchResult;
use crate::model::result::{
    CallabilitySegmentResult, ColumnResult, IntervalResult, ReferenceResult, ScopeResult,
};
use crate::model::review_result::{
    CONTEXT_RANK, Decision, EvidencePositionResult, EvidenceReadResult, GrowthResult, ItemKind,
    ItemResult, LabelCounts, LabelResult, MaskedSpanResult, MixtureResult, PositionCounts,
    ReviewCallResult, ReviewProvenanceResult, ReviewReadResult, ReviewResult, ReviewStatus, Role,
    RunResult, SecondResult, ShareResult, SideResult, SignalProvenanceResult, StrandsResult,
    SummaryResult, TractRunResult, UnreliableSpanResult, ZoneCounts, ZoneKind, ZoneResult,
};
use crate::model::sample_result::ConventionResult;
use crate::model::sanger_read::{Rejection, SangerRead};
use crate::report::json::{project_plugins, project_profile};
use crate::report::length;
use crate::report::notation::project_empop;
use crate::report::sample::{reviewer_names, reviewer_read_names};
use crate::sample_analysis::{
    CallTrust, ConsensusCall, NotationCall, ReviewCheck, ReviewList, SampleAnalysis, TractNote,
    TractNoteRun, TrustScope,
};
use dna_align::model::alignment::{AlignmentColumn, Orientation, ReferenceSegment};
use dna_call::model::called_read::CalledRead;
use dna_call::read_call::CallConfig;
use dna_consensus::decided::not_analysed;
use dna_consensus::model::consensus::{ConsensusRun, ConsensusSite, RunEvidence, SiteState};
use dna_consensus::runs::runs_of;
use dna_consensus::{ConsensusConfig, InsertionView};
use dna_kernel::model::nucleotide::{Nucleotide, is_canonical};
use dna_kernel::model::reference::Reference;
use dna_kernel::plugin::{Contract, PluginDescriptor, PluginFamily};
use dna_kernel::profile::{Convention, ConventionKind, ConventionLabel, LengthTract, Profile};
use dna_kernel::read_evidence::{EvidenceReason, Placement, ReadEvidence};
use dna_normalize::length_convention::ShortestDominant;
use dna_sample::model::sample_evidence::{LocusObservations, LocusState, SampleLocusObservation};
use dna_sanger::model::callability::PhaseState;
use dna_sanger::run_signal::{self, Key, RunSignal, Threshold};
use dna_sanger::signal_processing::SignalProcessingConfig;

/// The review list: items and labels over the consensus, the call, and the
/// reads (PROP-0004). Its version is the labelling rules' version.
pub(crate) const PLUGIN: PluginDescriptor = PluginDescriptor {
    id: "review",
    family: PluginFamily::PostCalling,
    version: 33,
    provides: &[Contract::ReviewList],
    requires: &[
        Contract::ReadEvidence,
        Contract::Consensus,
        Contract::CalledVariants,
    ],
    config_sections: &["review"],
};

/// Priority of an unlabelled variant, after every label's rank.
const UNLABELLED: u8 = CONTEXT_RANK + 1;

/// Inputs consumed to build one review document.
pub(crate) struct CompletedReview {
    pub(crate) sample_id: String,
    pub(crate) reference: Reference,
    pub(crate) profile: Profile,
    pub(crate) configuration_sha256: String,
    pub(crate) outcome: Reviewed,
    /// Plugins of the workflow, in execution order.
    pub(crate) plugins: &'static [&'static PluginDescriptor],
    /// The run signal the calls are judged against, with its file's SHA-256.
    pub(crate) signal: Option<(RunSignal, String)>,
    /// Each admitted read's phase mixtures, by content SHA-256 (ADR-0082).
    pub(crate) mixtures: BTreeMap<String, Vec<MixtureResult>>,
}

/// What the sample analysis and its call gave.
pub(crate) enum Reviewed {
    /// Every read was rejected; each with its cause, by content SHA-256.
    Rejected(BTreeMap<String, Rejection>),
    /// At least one read was admitted; the call, or why it failed.
    Analyzed {
        analysis: Box<SampleAnalysis>,
        call: Result<ConsensusCall>,
    },
}

/// Builds `dna.review/v1` without filesystem side effects.
pub(crate) fn build(completed: CompletedReview) -> Result<ReviewResult> {
    let CompletedReview {
        sample_id,
        reference,
        profile,
        configuration_sha256,
        outcome,
        plugins,
        signal,
        mixtures,
    } = completed;
    let regions = clip(&profile.regions, reference.len());
    let provenance = ReviewProvenanceResult {
        reference: ReferenceResult {
            name: reference.name.clone(),
            topology: reference.topology,
            sha256: reference.sequence_sha256.clone(),
        },
        configuration_sha256,
        profile: project_profile(profile.identity().clone()),
        plugins: project_plugins(plugins),
        rules: PLUGIN.version,
        signal: signal
            .as_ref()
            .map(|(signal, sha256)| SignalProvenanceResult {
                run_id: signal.run_id.clone(),
                sha256: sha256.clone(),
                selection: signal.selection,
            }),
        run_context: match &outcome {
            Reviewed::Analyzed { analysis, .. } => analysis
                .run_context()
                .map(crate::model::run_context_result::ContextUseResult::of),
            Reviewed::Rejected(_) => None,
        },
    };
    let (analysis, call) = match outcome {
        Reviewed::Rejected(rejected) => {
            let names = reviewer_names(
                rejected
                    .values()
                    .map(|read| (read.input_name.as_str(), read.input_sha256.as_str())),
            )?;
            let reads = rejected
                .values()
                .zip(names)
                .map(|(read, name)| ReviewReadResult {
                    name,
                    sha256: read.input_sha256.clone(),
                    strand: None,
                    admitted: false,
                    rejection: Some(read.cause),
                    masked: None,
                    unreliable: None,
                    mixtures: None,
                })
                .collect();
            return Ok(ReviewResult {
                schema_version: "dna.review/v1",
                sample_id,
                status: ReviewStatus::Rejected,
                error: None,
                provenance,
                reads,
                calls: None,
                callable: None,
                reference_only: None,
                not_analysed: None,
                empop: None,
                zones: None,
                summary: summarize(&[], ReviewStatus::Rejected, &regions, &[], &[], 0),
                items: Vec::new(),
            });
        }
        Reviewed::Analyzed { analysis, call } => (analysis, call),
    };

    let reference_only = analysis
        .consensus
        .reference_only
        .iter()
        .map(|position| position + 1)
        .filter(|&position| inside(&regions, position))
        .collect::<Vec<_>>();
    let (names, rejected_names) = reviewer_read_names(&analysis.evidence)?;
    let mut reads = analysis
        .evidence
        .reads
        .iter()
        .zip(&names)
        .map(|(read, name)| {
            let missing = || ReportError::Inconsistent("review read lacks its evidence");
            let called = analysis
                .called
                .iter()
                .find(|called| called.input_sha256 == read.input_sha256)
                .ok_or_else(missing)?;
            let sanger = analysis
                .sanger
                .get(&read.input_sha256)
                .ok_or_else(missing)?;
            Ok(ReviewReadResult {
                name: name.clone(),
                sha256: read.input_sha256.clone(),
                strand: Some(read.alignment.orientation),
                admitted: true,
                rejection: None,
                masked: Some(masked_spans(
                    called.alignment.columns(),
                    &sanger.callability.segments,
                    &analysis.trust(&read.input_sha256).ok_or_else(missing)?,
                    &regions,
                )),
                unreliable: Some(unreliable_spans(
                    called.alignment.columns(),
                    &called.evidence,
                    &regions,
                )),
                mixtures: Some(
                    mixtures
                        .get(&read.input_sha256)
                        .cloned()
                        .unwrap_or_default(),
                ),
            })
        })
        .collect::<Result<Vec<_>>>()?;
    for (read, name) in analysis.evidence.rejected_reads.iter().zip(rejected_names) {
        let rejected =
            analysis
                .rejections
                .get(&read.input_sha256)
                .ok_or(ReportError::Inconsistent(
                    "rejected review read lacks its rejection cause",
                ))?;
        reads.push(ReviewReadResult {
            name,
            sha256: read.input_sha256.clone(),
            strand: None,
            admitted: false,
            rejection: Some(rejected.cause),
            masked: None,
            unreliable: None,
            mixtures: None,
        });
    }

    let (status, error, call) = match call {
        Ok(call) => (ReviewStatus::Called, None, Some(call)),
        Err(error) => (ReviewStatus::CallFailed, Some(error.to_string()), None),
    };
    let callable = merge(match &call {
        Some(call) => segments(
            call.callable_segments
                .iter()
                .map(|segment| (segment.start_0based, segment.end_0based_exclusive)),
        ),
        None => segments(
            analysis
                .evidence
                .reads
                .iter()
                .flat_map(|read| &read.alignment.callable_reference_segments)
                .map(span),
        ),
    });
    let notation = call.as_ref().and_then(|call| call.notation.as_deref());
    // The called sample's analysed regions: the reportable regions without
    // the stretches at their edges the call does not analyse.
    let (analysed, not_analysed) = if call.is_some() {
        analysed_regions(&analysis, &regions)
    } else {
        (regions.clone(), Vec::new())
    };
    let items = Sample::new(
        &analysis,
        &names,
        &reference,
        &analysed,
        &callable,
        notation,
        signal.as_ref().map(|(signal, _)| signal),
    )?
    .items()?;
    let calls = call.as_ref().and_then(|call| {
        call.notation.as_ref().map(|calls| {
            calls
                .iter()
                .map(|call| ReviewCallResult {
                    call: call.text.clone(),
                    position: call.position,
                    ordinal: call.ordinal,
                    reference: call.reference,
                    alternate: call.alternate,
                    convention: call.convention.as_ref().map(|convention| ConventionResult {
                        name: convention.name.clone(),
                        label: convention.label.label(),
                    }),
                })
                .collect()
        })
    });
    let windows = profile
        .windows
        .iter()
        .map(|window| (window.start + 1, window.end(), window.name.as_str()))
        .collect::<Vec<_>>();
    let zones = call
        .as_ref()
        .map(|_| zones(&analysed, &windows, &callable, &reads, &items));
    Ok(ReviewResult {
        schema_version: "dna.review/v1",
        sample_id,
        status,
        error,
        provenance,
        reads,
        calls,
        callable: call.as_ref().map(|_| {
            callable
                .iter()
                .map(|&[start, end]| IntervalResult { start, end })
                .collect()
        }),
        summary: summarize(
            &items,
            status,
            &regions,
            zones.as_deref().unwrap_or_default(),
            if call.is_some() { &reference_only } else { &[] },
            not_analysed
                .iter()
                .map(|stretch| stretch.last - stretch.first + 1)
                .sum(),
        ),
        reference_only: call.as_ref().map(|_| reference_only.clone()),
        not_analysed: call.as_ref().map(|_| not_analysed.clone()),
        empop: call.and_then(|call| call.empop).map(project_empop),
        zones,
        items,
    })
}

/// The review list of an analysed sample whose consensus was called, by the
/// rules of `dna review`: each check item's first and last reference
/// positions (1-based) and label names, and each `tract_length` label with
/// its item's span.
pub(crate) fn list(
    analysis: &SampleAnalysis,
    call: &ConsensusCall,
    signal: Option<&RunSignal>,
) -> Result<ReviewList> {
    let reference = &analysis.loaded.reference;
    let regions = clip(&analysis.loaded.profile.regions, reference.len());
    let (names, _) = reviewer_read_names(&analysis.evidence)?;
    let callable =
        merge(segments(call.callable_segments.iter().map(|segment| {
            (segment.start_0based, segment.end_0based_exclusive)
        })));
    let (analysed, _) = analysed_regions(analysis, &regions);
    let items = Sample::new(
        analysis,
        &names,
        reference,
        &analysed,
        &callable,
        call.notation.as_deref(),
        signal,
    )?
    .items()?;
    let mut review = ReviewList::default();
    for item in &items {
        let (first, last) = (item.first.position, item.last.position);
        if item.kind == ItemKind::Check {
            review.checks.push(ReviewCheck {
                first,
                last,
                labels: item
                    .labels
                    .iter()
                    .map(|label| label.kind().name())
                    .collect::<BTreeSet<_>>()
                    .into_iter()
                    .collect(),
            });
        }
        for label in &item.labels {
            if let LabelResult::TractLength { tract, runs, .. } = label {
                review.notes.push(TractNote {
                    first,
                    last,
                    tract: tract.clone(),
                    runs: runs
                        .iter()
                        .map(|run| TractNoteRun {
                            run: run.run,
                            called: run.called,
                            status: run.status,
                            reading: run.reading,
                            convention: run.convention.as_ref().map(|convention| convention.length),
                        })
                        .collect(),
                });
            }
        }
    }
    Ok(review)
}

/// The reportable `regions` (1-based, inclusive) without the stretches at
/// their ends the call does not analyse (SRS-CONS-012), and those
/// stretches. A region left out whole is dropped.
fn analysed_regions(
    analysis: &SampleAnalysis,
    regions: &[[usize; 2]],
) -> (Vec<[usize; 2]>, Vec<StretchResult>) {
    let stretches = not_analysed(&analysis.consensus, regions);
    let left_out = |position: usize| {
        stretches
            .iter()
            .any(|stretch| stretch.first_0based < position && position <= stretch.last_0based + 1)
    };
    let mut analysed: Vec<[usize; 2]> = Vec::new();
    for &[start, end] in regions {
        for position in (start..=end).filter(|&position| !left_out(position)) {
            match analysed.last_mut() {
                Some(last) if last[1] + 1 == position && last[0] >= start => last[1] = position,
                _ => analysed.push([position, position]),
            }
        }
    }
    let stretches = stretches
        .iter()
        .map(|stretch| StretchResult {
            first: stretch.first_0based + 1,
            last: stretch.last_0based + 1,
            reason: stretch.reason.label(),
        })
        .collect();
    (analysed, stretches)
}

/// The run signal without the sample's own admitted reads, so the sample is
/// never judged against itself. A run signal that holds none of the sample's
/// traces is another run's.
fn own(analysis: &SampleAnalysis, signal: &RunSignal) -> Result<RunSignal> {
    if !analysis
        .called
        .iter()
        .any(|read| signal.traces.contains(&read.input_sha256))
    {
        return Err(RunSignalError::OtherRun {
            run_id: signal.run_id.clone(),
        }
        .into());
    }
    let parameters = &analysis.loaded.config.sanger.signal_processing.run_signal;
    signal.check(parameters)?;
    let mut own = signal.clone();
    for read in &analysis.called {
        if !signal.traces.contains(&read.input_sha256) {
            continue;
        }
        let sanger = analysis
            .sanger
            .get(&read.input_sha256)
            .ok_or(ReportError::Inconsistent("review read lacks its evidence"))?;
        own = own.without(&run_signal::contributions(
            &sanger.heights(),
            &read.evidence,
            parameters,
        ));
    }
    Ok(own)
}

/// One sample's evidence, indexed for the review rules.
struct Sample<'a> {
    /// The review plugin's thresholds (`[review]`).
    config: ReviewConfig,
    /// The profile's report conventions: a run length one covers is not
    /// recorded, so its estimate leaves nothing to assure.
    conventions: &'a [Convention],
    sequence: &'a [u8],
    /// The analysed regions: the reportable regions without the stretches
    /// the call does not analyse (SRS-CONS-012).
    regions: &'a [[usize; 2]],
    /// The profile's reportable regions, whose ends are the amplicons'.
    reportable: Vec<[usize; 2]>,
    /// The stretches at the reportable regions' ends the call does not
    /// analyse (SRS-CONS-012), 1-based and inclusive.
    not_analysed: Vec<[usize; 2]>,
    callable: &'a [[usize; 2]],
    /// Reference segments of the admitted reads, merged, 1-based inclusive.
    covered: Vec<[usize; 2]>,
    calls: &'a [NotationCall],
    /// Whether DNA's calls are in the profile notation; without them, a
    /// difference is the consensus's decision.
    notation: bool,
    sites: &'a [ConsensusSite],
    /// The positions (0-based) reads observe only inside their unreliable
    /// spans (SRS-CONS-012).
    unreliable_only: &'a BTreeMap<usize, EvidenceReason>,
    /// Per read, in evidence order: name, content SHA-256, strand, the core's
    /// trust in each call, and the Sanger evidence.
    reads: Vec<Read<'a>>,
    /// Read index by content SHA-256, as consensus sites name reads.
    by_id: BTreeMap<&'a str, usize>,
    /// Every covered position's observations.
    loci: BTreeMap<usize, LocusObservations>,
    /// The run signal the calls are judged against, the sample's own reads
    /// taken out, when one is given.
    signal: Option<RunSignal>,
    /// The second-channel threshold's floor and run-level parameters.
    second: &'a SignalProcessingConfig,
    /// The profile's anchored repeat windows and length tracts, whose run
    /// lengths the reads' tract decodes judge (SRS-API-020).
    length_windows: Vec<length::LengthWindow>,
    /// The profile's length tracts.
    length_tracts: &'a [LengthTract],
    /// The profile's nomenclature windows and length tracts, 1-based and
    /// inclusive: where a run's length is judged by its own rules.
    repeat_windows: Vec<[usize; 2]>,
    /// The length convention, when configured (`[normalization]`).
    length_convention: Option<ShortestDominant>,
    /// The reference and the core's and the consensus's configuration, for
    /// the consensus's reading of each read's runs (ADR-0074).
    reference: &'a Reference,
    call_config: &'a CallConfig,
    consensus_config: &'a ConsensusConfig,
}

struct Read<'a> {
    name: &'a str,
    /// The read's alignment columns, in reference order.
    columns: &'a [AlignmentColumn],
    sha256: &'a str,
    strand: Orientation,
    trust: Vec<CallTrust>,
    sanger: &'a SangerRead,
    /// The read's calls as its evidence gives them, in read order.
    evidence: &'a ReadEvidence,
    /// The read as the core called it.
    called: &'a CalledRead,
}

/// An item before its labels are complete.
struct Draft<'a> {
    first: usize,
    last: usize,
    site: Option<&'a ConsensusSite>,
    /// The decided sequence: a site's call, or the reference base of a
    /// position outside every site.
    sequence: Option<String>,
    /// Whether the consensus decided the item without an `N`.
    decided: bool,
    labels: Vec<LabelResult>,
}

impl<'a> Sample<'a> {
    fn new(
        analysis: &'a SampleAnalysis,
        names: &'a [String],
        reference: &'a Reference,
        regions: &'a [[usize; 2]],
        callable: &'a [[usize; 2]],
        notation: Option<&'a [NotationCall]>,
        signal: Option<&RunSignal>,
    ) -> Result<Self> {
        let signal = signal.map(|signal| own(analysis, signal)).transpose()?;
        let reads = analysis
            .evidence
            .reads
            .iter()
            .zip(names)
            .map(|(read, name)| {
                let missing = || ReportError::Inconsistent("review read lacks its evidence");
                let called = analysis
                    .called
                    .iter()
                    .find(|called| called.input_sha256 == read.input_sha256)
                    .ok_or_else(missing)?;
                Ok(Read {
                    name,
                    columns: called.alignment.columns(),
                    evidence: &called.evidence,
                    called,
                    sha256: &read.input_sha256,
                    strand: read.alignment.orientation,
                    trust: analysis.trust(&read.input_sha256).ok_or_else(missing)?,
                    sanger: analysis
                        .sanger
                        .get(&read.input_sha256)
                        .ok_or_else(missing)?,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let by_id = analysis
            .evidence
            .reads
            .iter()
            .enumerate()
            .map(|(index, read)| (read.input_sha256.as_str(), index))
            .collect();
        let covered = merge(segments(
            analysis
                .evidence
                .reads
                .iter()
                .flat_map(|read| &read.alignment.reference_segments)
                .map(span),
        ));
        // The reportable positions outside the analysed `regions`: the
        // stretches the call does not analyse.
        let reportable = clip(&analysis.loaded.profile.regions, reference.len());
        let mut not_analysed: Vec<[usize; 2]> = Vec::new();
        for position in reportable
            .iter()
            .flat_map(|&[start, end]| start..=end)
            .filter(|&position| !inside(regions, position))
        {
            match not_analysed.last_mut() {
                Some(last) if last[1] + 1 == position => last[1] = position,
                _ => not_analysed.push([position, position]),
            }
        }
        let loci = analysis
            .observations(1, reference.len())?
            .into_iter()
            .map(|locus| (locus.position_1based, locus))
            .collect();
        Ok(Self {
            config: analysis.loaded.config.review,
            conventions: &analysis.loaded.profile.conventions,
            sequence: reference.sequence.as_bytes(),
            regions,
            reportable,
            not_analysed,
            callable,
            covered,
            calls: notation.unwrap_or_default(),
            notation: notation.is_some(),
            sites: &analysis.consensus.sites,
            unreliable_only: &analysis.consensus.unreliable_only,
            reads,
            by_id,
            loci,
            signal,
            second: &analysis.loaded.config.sanger.signal_processing,
            length_windows: length::windows(&analysis.loaded.profile),
            length_tracts: &analysis.loaded.profile.length_tracts,
            repeat_windows: analysis
                .loaded
                .profile
                .windows
                .iter()
                .map(|window| [window.start + 1, window.end()])
                .chain(
                    analysis
                        .loaded
                        .profile
                        .length_tracts
                        .iter()
                        .map(|tract| tract.window),
                )
                .collect(),
            length_convention: analysis.loaded.config.length_convention,
            reference,
            call_config: &analysis.loaded.config.core,
            consensus_config: &analysis.loaded.config.consensus,
        })
    }

    /// Every item, in priority order.
    fn items(&self) -> Result<Vec<ItemResult>> {
        let in_site = |position: usize| {
            self.sites
                .iter()
                .any(|site| site.first_0based < position && position <= site.last_0based + 1)
        };
        let mut drafts = Vec::new();
        // A: every consensus site, clipped to the regions.
        for site in self.sites {
            if let Some([first, last]) =
                within(self.regions, site.first_0based + 1, site.last_0based + 1)
            {
                let mut labels = self.site_labels(site);
                let mut sequence = site.call.clone();
                if let Some(shown) = self.clipped_decision(site, first, last) {
                    self.clip_disagreement(&mut labels, first, last, &shown);
                    sequence = Some(shown);
                }
                drafts.push(Draft {
                    first,
                    last,
                    site: Some(site),
                    sequence,
                    decided: site.call.as_ref().is_some_and(|call| !call.contains('N')),
                    labels,
                });
            }
        }
        // B: callable positions outside every site where a second base shows;
        // C: maximal runs of positions outside every site and not callable.
        let mut run: Option<[usize; 2]> = None;
        for position in self.regions.iter().flat_map(|&[start, end]| start..=end) {
            let open = !in_site(position) && !inside(self.callable, position);
            match (&mut run, open) {
                (Some([_, last]), true) if *last + 1 == position => *last = position,
                (_, true) => {
                    if let Some([first, last]) = run.replace([position, position]) {
                        drafts.push(self.gap(first, last));
                    }
                }
                (_, false) => {
                    if let Some([first, last]) = run.take() {
                        drafts.push(self.gap(first, last));
                    }
                }
            }
            if !in_site(position)
                && inside(self.callable, position)
                && let Some(label) = self.heteroplasmy(position)
            {
                drafts.push(Draft {
                    first: position,
                    last: position,
                    site: None,
                    sequence: Some(self.base(position).to_string()),
                    decided: true,
                    labels: vec![label],
                });
            }
        }
        if let Some([first, last]) = run {
            drafts.push(self.gap(first, last));
        }
        // A call no other item holds, for example one the notation moved.
        let held = drafts
            .iter()
            .map(|draft| [draft.first, draft.last])
            .collect::<Vec<_>>();
        let lone = self
            .calls
            .iter()
            .map(|call| call.position)
            .filter(|&position| {
                inside(self.regions, position)
                    && !held
                        .iter()
                        .any(|&[first, last]| first <= position && position <= last)
            })
            .collect::<BTreeSet<_>>();
        for position in lone {
            drafts.push(Draft {
                first: position,
                last: position,
                site: None,
                sequence: None,
                decided: true,
                labels: Vec::new(),
            });
        }

        self.length_labels(&mut drafts);
        let mut items = drafts
            .into_iter()
            .map(|draft| self.finish(draft))
            .filter_map(Result::transpose)
            .collect::<Result<Vec<_>>>()?;
        items.sort_by_key(|item| (priority(&item.labels), item.first, item.last));
        Ok(items)
    }

    /// Labels the runs of each judged window (SRS-API-020). In a window whose
    /// run lengths the profile records, every run the reads' tract decodes
    /// do not confirm is `length_unconfirmed`, a check, unless it lies in a
    /// profile length tract and its dominant length is determined
    /// (`length::status`): such a run is noted in `tract_length`, context.
    /// In a window whose run lengths the profile does not record, every
    /// decoded run is `length_not_recorded`, context. An `estimated_length`
    /// or `long_run` run of the window's item that is a tract run so
    /// determined moves to the note too; a plain tract (no anchored window)
    /// gives a note only so. The labels join the first item that overlaps
    /// the window, or a new item over the window. A window outside the
    /// analysed regions gets none.
    fn length_labels(&self, drafts: &mut Vec<Draft<'a>>) {
        for window in &self.length_windows {
            let Some([first, last]) = within(self.regions, window.first, window.last) else {
                continue;
            };
            let decodes = self
                .reads
                .iter()
                .flat_map(|read| {
                    read.evidence.tracts().iter().map(|tract| length::Decode {
                        read: read.name,
                        strand: read.strand,
                        reference: dna_consensus::length::start_reference(
                            read.columns,
                            tract.call_start,
                        ),
                        tract,
                    })
                })
                .collect::<Vec<_>>();
            let haplotype = self.decided_over(window.first, window.last);
            let called = haplotype
                .as_deref()
                .and_then(|haplotype| length::called_runs(haplotype, window));
            let called_reads = self
                .reads
                .iter()
                .map(|read| read.called)
                .collect::<Vec<_>>();
            let views = dna_consensus::read_runs_over(
                &called_reads,
                self.reference,
                self.call_config,
                self.consensus_config,
                window.first - 1,
                window.last - 1,
            );
            let in_phase = self
                .reads
                .iter()
                .zip(&views)
                .filter_map(|(read, runs)| {
                    Some(length::InPhase {
                        read: read.name,
                        strand: read.strand,
                        runs: runs.as_deref()?,
                    })
                })
                .collect::<Vec<_>>();
            let two_strand = self
                .sites
                .iter()
                .find(|site| {
                    site.first_0based < window.first && window.last <= site.last_0based + 1
                })
                .map(|site| length::two_strand_runs(&site.runs, window))
                .unwrap_or_default();
            // DNA reports a run's length only where the run and the
            // positions that bound it are callable; elsewhere its decided
            // length is an estimate. Without DNA's calls in the profile's
            // notation, its runs are not known to lie around the anchor.
            let anchored = self.notation && length::anchored(window, self.calls);
            let reported = length::run_spans(window, called.as_deref(), anchored)
                .into_iter()
                .enumerate()
                .filter(|&(_, [start, end])| {
                    within(self.regions, start, end).is_some_and(|[start, end]| {
                        (start..=end).all(|position| inside(self.callable, position))
                    })
                })
                .map(|(index, _)| index + 1)
                .collect::<Vec<_>>();
            let judged = length::judge(
                window,
                called.as_deref(),
                &reported,
                &two_strand,
                &decodes,
                &in_phase,
                self.length_convention,
            );
            let single = called.as_ref().is_some_and(|runs| runs.len() == 1);
            let draft = drafts
                .iter()
                .position(|draft| draft.first <= last && first <= draft.last);
            // A tract run DNA's `estimated_length` or `long_run` names: the one
            // run of the window with that base and length.
            let named = |base: char, length: usize| {
                let matching = called
                    .iter()
                    .flatten()
                    .enumerate()
                    .filter(|&(_, &decided)| base == window.repeat && decided == length)
                    .map(|(index, _)| index + 1)
                    .collect::<Vec<_>>();
                match matching.as_slice() {
                    [run] if window.in_tract(*run, single) => Some(*run),
                    _ => None,
                }
            };
            let mut noted = BTreeSet::new();
            let mut moved = Vec::new();
            if let Some(index) = draft {
                for label in &drafts[index].labels {
                    let runs = match label {
                        LabelResult::EstimatedLength { runs } => runs
                            .iter()
                            .map(|run| (run.base, run.length))
                            .collect::<Vec<_>>(),
                        LabelResult::LongRun { runs } => {
                            runs.iter().map(|run| (run.base, run.length)).collect()
                        }
                        _ => Vec::new(),
                    };
                    for (base, length) in runs {
                        let Some(run) = named(base, length) else {
                            continue;
                        };
                        let judged_run = judged.iter().find(|judged| judged.result.run == run);
                        // An estimate no read reaches is unconfirmed; one the
                        // reads judge moves only when they determine it.
                        let determined = judged_run.map_or(!reported.contains(&run), |judged| {
                            length::status(&judged.result).is_some()
                        });
                        if determined {
                            noted.insert(run);
                            moved.push((base, length));
                        }
                    }
                }
            }
            let mut unconfirmed = Vec::new();
            let mut not_recorded = Vec::new();
            for judged_run in &judged {
                let run = judged_run.result.run;
                let tract = window.in_tract(run, single);
                let status = tract.then(|| length::status(&judged_run.result)).flatten();
                match (window.anchor.is_some(), window.recorded, status) {
                    (true, false, _) => {
                        let mut result = judged_run.result.clone();
                        // The lab records no length here, so no reading of it.
                        result.convention = None;
                        result.floor = None;
                        not_recorded.push(result);
                    }
                    (true, true, Some((status, _))) if status != "confirmed" => {
                        noted.insert(run);
                    }
                    (true, true, _) if judged_run.result.agreement != "confirmed" => {
                        unconfirmed.push(judged_run.result.clone());
                    }
                    // A plain tract labels only what it moves.
                    _ => {}
                }
            }
            let mut labels = Vec::new();
            if !unconfirmed.is_empty() {
                labels.push(LabelResult::LengthUnconfirmed { runs: unconfirmed });
            }
            if !not_recorded.is_empty() {
                labels.push(LabelResult::LengthNotRecorded { runs: not_recorded });
            }
            if let Some(tract) = window
                .tract
                .and_then(|tract| self.length_tracts.get(tract.index))
                .filter(|_| !noted.is_empty())
            {
                let runs = noted
                    .iter()
                    .map(|&run| {
                        let judged_run = judged.iter().find(|judged| judged.result.run == run);
                        tract_run(
                            run,
                            window,
                            called.as_deref(),
                            reported.contains(&run),
                            judged_run,
                        )
                    })
                    .collect();
                labels.push(LabelResult::TractLength {
                    tract: tract.name.clone(),
                    ignored_in_search: tract.search,
                    runs,
                });
            }
            if labels.is_empty() {
                continue;
            }
            match draft {
                Some(index) => {
                    let draft = &mut drafts[index];
                    for label in &mut draft.labels {
                        match label {
                            LabelResult::EstimatedLength { runs } => {
                                runs.retain(|run| !moved.contains(&(run.base, run.length)));
                            }
                            LabelResult::LongRun { runs } => {
                                runs.retain(|run| !moved.contains(&(run.base, run.length)));
                            }
                            _ => {}
                        }
                    }
                    draft.labels.retain(|label| {
                        !matches!(
                            label,
                            LabelResult::EstimatedLength { runs } if runs.is_empty()
                        ) && !matches!(label, LabelResult::LongRun { runs } if runs.is_empty())
                    });
                    draft.labels.extend(labels);
                }
                None => drafts.push(Draft {
                    first,
                    last,
                    site: None,
                    decided: haplotype.as_deref().is_some_and(|text| !text.contains('N')),
                    sequence: haplotype,
                    labels,
                }),
            }
        }
    }

    /// DNA's decided sequence over `first..=last` (1-based), widened to the
    /// consensus sites it cuts, with each site's call in place of its
    /// reference; `None` when a site there has no call.
    fn decided_over(&self, first: usize, last: usize) -> Option<String> {
        let (mut start, mut end) = (first - 1, last);
        for site in self.sites {
            if site.first_0based < end && start <= site.last_0based {
                start = start.min(site.first_0based);
                end = end.max(site.last_0based + 1);
            }
        }
        let mut text = String::new();
        let mut position = start;
        while position < end {
            if let Some(site) = self.sites.iter().find(|site| site.first_0based == position) {
                text.push_str(site.call.as_deref()?);
                position = site.last_0based + 1;
            } else {
                text.push(char::from(*self.sequence.get(position)?));
                position += 1;
            }
        }
        Some(text)
    }

    /// A run of positions outside every site that are not callable.
    fn gap(&self, first: usize, last: usize) -> Draft<'a> {
        let covered = self
            .covered
            .iter()
            .any(|&[start, end]| start <= last && first <= end);
        Draft {
            first,
            last,
            site: None,
            sequence: None,
            decided: false,
            labels: vec![if covered {
                LabelResult::Undecided
            } else {
                LabelResult::NotCovered
            }],
        }
    }

    /// The labels a consensus site carries by itself.
    fn site_labels(&self, site: &ConsensusSite) -> Vec<LabelResult> {
        let mut labels = Vec::new();
        let sides = |names: &[String], role: Role| {
            names
                .iter()
                .filter_map(|input| self.by_id.get(input.as_str()))
                .map(|&index| SideResult {
                    read: self.reads[index].name.to_owned(),
                    strand: self.reads[index].strand,
                    role,
                })
                .collect::<Vec<_>>()
        };
        labels.extend(dissent(
            sides(&site.supporting, Role::Supporting),
            sides(&site.opposing, Role::Opposing),
        ));
        if undecided(site) {
            labels.push(LabelResult::Undecided);
        }
        let runs = site
            .runs
            .iter()
            .filter(|run| run.evidence.is_estimate() && !not_recorded(self.conventions, site, run))
            .map(|run| RunResult {
                base: run.base,
                length: run.length,
                length_evidence: run.evidence.label(),
            })
            .collect::<Vec<_>>();
        if !runs.is_empty() {
            labels.push(LabelResult::EstimatedLength { runs });
        }
        let anchored = site
            .runs
            .iter()
            .filter(|run| run.evidence == RunEvidence::AnchoredEnd)
            .map(|run| RunResult {
                base: run.base,
                length: run.length,
                length_evidence: run.evidence.label(),
            })
            .collect::<Vec<_>>();
        if !anchored.is_empty() {
            labels.push(LabelResult::DephasedEnd { runs: anchored });
        }
        if !undecided(site) {
            let support = self.one_strand_support(site);
            let growth = long_runs(site, support.is_some(), self.config.run_growth);
            if !growth.is_empty() {
                labels.push(LabelResult::LongRun { runs: growth });
            }
            if let Some((strand, reads)) = support {
                if self.run_deletion(site) {
                    labels.push(LabelResult::RunDeletion {
                        strand,
                        reads: reads.clone(),
                    });
                }
                if self.run_indel(site) {
                    labels.push(LabelResult::RunIndel {
                        strand,
                        reads: reads.clone(),
                    });
                }
                if self.edge_deletion(site) {
                    labels.push(LabelResult::EdgeDeletion { strand, reads });
                }
            }
        }
        labels
    }

    /// The strand every supporting read of `site` lies on, however many, and
    /// their names in the site's order (ADR-0090); `None` without a
    /// supporting read, with one the sample does not hold, or with
    /// supporting reads on both strands. Reads of one strand share its
    /// systematic errors, so their number does not assure what they show.
    fn one_strand_support(&self, site: &ConsensusSite) -> Option<(Orientation, Vec<String>)> {
        let reads = site
            .supporting
            .iter()
            .map(|input| {
                self.by_id
                    .get(input.as_str())
                    .map(|&index| &self.reads[index])
            })
            .collect::<Option<Vec<_>>>()?;
        let strand = one_strand(&reads.iter().map(|read| read.strand).collect::<Vec<_>>())?;
        Some((
            strand,
            reads.iter().map(|read| read.name.to_owned()).collect(),
        ))
    }

    /// Whether a site that rests on one strand (its caller checks) has a
    /// call at least `review.run_growth` bases shorter than its reference,
    /// where the reference holds a run of at least `review.phase_run` equal
    /// bases: a read that loses or compresses part of a run aligns it as a
    /// long deletion, and reads of one strand cannot assure that.
    fn run_deletion(&self, site: &ConsensusSite) -> bool {
        site.call.as_deref().is_some_and(|call| {
            run_lost(
                &site.reference,
                call,
                self.config.run_growth,
                self.config.phase_run,
            )
        })
    }

    /// Whether a site that rests on one strand (its caller checks) changes
    /// the length of a homopolymer run of its reference of at least
    /// `review.phase_run` equal bases (inserts or deletes copies of the run's
    /// base, the runs' bases unchanged) and lies outside every repeat window
    /// and length tract of the profile. A strand slips on such a run, and the
    /// other strand often reads it only out of phase, so one strand cannot
    /// assure its length.
    fn run_indel(&self, site: &ConsensusSite) -> bool {
        let Some(call) = site.call.as_deref() else {
            return false;
        };
        let (first, last) = (site.first_0based + 1, site.last_0based + 1);
        !self
            .repeat_windows
            .iter()
            .any(|&[start, end]| start <= last && first <= end)
            && run_length_change(&site.reference, call, self.config.phase_run)
    }

    /// `one_strand_insertion` at a decided insertion (a notation call with
    /// an ordinal over `first..=last`, or a decision longer than its
    /// reference without notation) outside every repeat window and length
    /// tract, whose supporting reads `supporting` all lie on `strand`, when
    /// no read of the other strand shows it by the consensus's reading
    /// (SRS-CONS-015): the other strand does not reach the site, ends in it,
    /// shows no base there, or reads it through without the insertion where
    /// the consensus leaves a long repeat's length to its run rules. One
    /// strand cannot assure an insertion; on development panel v5 such
    /// insertions near a read's start were often its own artefacts.
    fn one_strand_insertion(
        &self,
        site: &ConsensusSite,
        first: usize,
        last: usize,
        strand: Orientation,
        supporting: &[usize],
    ) -> Option<LabelResult> {
        let call = site.call.as_deref()?;
        let inserts = if self.notation {
            self.calls.iter().any(|notation| {
                notation.ordinal > 0 && first <= notation.position && notation.position <= last
            })
        } else {
            call.len() > site.reference.len()
        };
        if !inserts
            || self
                .repeat_windows
                .iter()
                .any(|&[start, end]| start <= last && first <= end)
        {
            return None;
        }
        let called = self
            .reads
            .iter()
            .map(|read| read.called)
            .collect::<Vec<_>>();
        let views = dna_consensus::insertion_views(
            &called,
            self.reference,
            self.call_config,
            self.consensus_config,
            site.first_0based,
            site.last_0based,
            call,
        );
        let other = self
            .reads
            .iter()
            .zip(&views)
            .filter(|(read, _)| read.strand != strand)
            .map(|(_, view)| *view)
            .max()
            .unwrap_or(InsertionView::NotCovered);
        if other == InsertionView::Shows {
            return None;
        }
        // The fewest calls from a supporting read's first call to its first
        // call over the site.
        let from_start = supporting
            .iter()
            .filter_map(|&index| {
                // A call aligned over the site, or inserted after a
                // position of it.
                let mut within = false;
                self.reads[index]
                    .columns
                    .iter()
                    .filter(|column| {
                        if let Some(reference) = column.reference_index_0based() {
                            within =
                                site.first_0based <= reference && reference <= site.last_0based;
                        }
                        within
                    })
                    .filter_map(AlignmentColumn::original_call_index_0based)
                    .min()
            })
            .min()?;
        Some(LabelResult::OneStrandInsertion {
            strand,
            reads: supporting
                .iter()
                .map(|&index| self.reads[index].name.to_owned())
                .collect(),
            other: other.label(),
            from_start,
        })
    }

    /// Whether a site that rests on one strand (its caller checks) is a
    /// deletion within `deletion_zone` positions of a reportable region's
    /// first or last position, where a strand's own read edge loses or
    /// doubles a base.
    fn edge_deletion(&self, site: &ConsensusSite) -> bool {
        site.call
            .as_deref()
            .is_some_and(|call| call.len() < site.reference.len())
            && near_edge(
                &self.reportable,
                site.first_0based + 1,
                site.last_0based + 1,
                self.config.deletion_zone,
            )
    }

    /// Completes an item's calls, labels, and evidence; `None` for a site
    /// without a call or a label.
    fn finish(&self, draft: Draft<'_>) -> Result<Option<ItemResult>> {
        let Draft {
            first,
            last,
            site,
            sequence,
            decided,
            mut labels,
        } = draft;
        let calls = self
            .calls
            .iter()
            .filter(|call| first <= call.position && call.position <= last)
            .map(|call| call.text.clone())
            .collect::<Vec<_>>();
        // A difference: DNA's calls when it has them in the profile notation,
        // otherwise the consensus's own decision. A decided difference that
        // DNA does not call within a position of the item is `not_called`.
        let decided_difference = site
            .is_some_and(|site| decided && site.call.as_deref() != Some(site.reference.as_str()));
        let differs = if self.notation {
            !calls.is_empty()
        } else {
            decided_difference
        };
        if self.notation && decided_difference && !called_near(self.calls, first, last) {
            labels.push(LabelResult::NotCalled);
        }
        if let (Some(site), Some(sequence), true) = (site, sequence.as_deref(), decided) {
            labels.extend(self.disagreeing(site, first, last, sequence));
            let ends = self.end_disagreeing(first, last, site, sequence, &labels)?;
            labels.extend(ends);
        }
        if let (true, Some(site)) = (differs, site) {
            let indices = site
                .supporting
                .iter()
                .filter_map(|input| self.by_id.get(input.as_str()).copied())
                .collect::<Vec<_>>();
            let supporting = indices
                .iter()
                .map(|&index| &self.reads[index])
                .collect::<Vec<_>>();
            if let [read] = supporting.as_slice() {
                labels.push(LabelResult::OneRead {
                    read: read.name.to_owned(),
                });
            }
            if let Some(strand) = one_strand(
                &supporting
                    .iter()
                    .map(|read| read.strand)
                    .collect::<Vec<_>>(),
            ) {
                labels.push(LabelResult::OneStrand {
                    strand,
                    reads: supporting.iter().map(|read| read.name.to_owned()).collect(),
                });
                for &index in &indices {
                    labels.extend(self.reference_second(index, first, last)?);
                }
                if decided {
                    labels.extend(self.one_strand_insertion(site, first, last, strand, &indices));
                }
            }
        }
        if differs
            && beside_not_analysed(
                &self.reportable,
                &self.not_analysed,
                first,
                last,
                self.config.edge_zone,
            )
        {
            labels.push(LabelResult::BesideNotAnalysed);
        }
        labels.extend((first..=last).filter_map(|position| self.heteroplasmy(position)));
        let (read_labels, evidence) = self.reads_over(first, last)?;
        labels.extend(read_labels);
        // A label that leaves the item unassured makes it a check; the others
        // are context for the reviewer. An item with neither a call nor such a
        // label is DNA's own decision and is not listed.
        // A difference inside a read's poor signal is a check too: there the
        // reads' differences are often wrong (SRS-CALL-024).
        let poor = labels
            .iter()
            .any(|label| matches!(label, LabelResult::PoorSignal { .. }));
        let assured = !(labels.iter().any(LabelResult::unassured) || (differs && poor));
        let noted = labels.iter().any(|label| {
            matches!(
                label,
                LabelResult::PossibleHeteroplasmy {
                    explained: true,
                    ..
                } | LabelResult::TractLength { .. }
            )
        });
        let kind = match (assured, differs, noted) {
            (false, _, _) => ItemKind::Check,
            (true, true, _) => ItemKind::Variant,
            (true, false, true) => ItemKind::Note,
            (true, false, false) => return Ok(None),
        };
        labels.sort_by(|left, right| key(left).cmp(&key(right)));
        labels.dedup();
        let callable = (first..=last)
            .filter(|&position| inside(self.callable, position))
            .count();
        let decision = if callable == 0 {
            Decision::NotCallable
        } else if decided && callable == last - first + 1 {
            Decision::Called
        } else {
            Decision::PartlyCallable
        };
        Ok(Some(ItemResult {
            first: ColumnResult {
                position: first,
                ordinal: 0,
            },
            last: ColumnResult {
                position: last,
                ordinal: 0,
            },
            kind,
            decision,
            sequence,
            calls,
            labels,
            evidence,
        }))
    }

    /// The labels of the reads' calls over `first..=last` that are not
    /// clean, one per label and read, and every covering read's evidence.
    fn reads_over(
        &self,
        first: usize,
        last: usize,
    ) -> Result<(Vec<LabelResult>, Vec<EvidenceReadResult>)> {
        let mut labels = BTreeMap::new();
        let mut positions = BTreeMap::<usize, Vec<EvidencePositionResult>>::new();
        for position in first..=last {
            let Some(locus) = self.loci.get(&position) else {
                continue;
            };
            for observation in &locus.observations {
                let index = observation.read_index;
                let read = self.read(index)?;
                if let Some(call) = observation.call_index_0based {
                    let trust = read.trust.get(call).ok_or(ReportError::Inconsistent(
                        "review observation lacks the trust of its call",
                    ))?;
                    let label = match trust.placement {
                        Placement::Anchoring(_) => Some(LabelResult::PhaseShift {
                            read: read.name.to_owned(),
                            from: self.shift_start(index, position),
                        }),
                        Placement::Unresolved(reason) => Some(LabelResult::NoisySignal {
                            read: read.name.to_owned(),
                            reason: reason.label(),
                        }),
                        Placement::Scored if !trust.trusted && !trust.exempt => {
                            Some(LabelResult::ReadEdge {
                                read: read.name.to_owned(),
                            })
                        }
                        Placement::Scored => None,
                    };
                    if let Some(label) = label {
                        labels
                            .entry((label.rank(), label.kind().name(), index))
                            .or_insert(label);
                    }
                    // Inside the read's unreliable span, where no read shows
                    // the position soundly outside its own (SRS-CONS-012).
                    if read
                        .evidence
                        .unreliable()
                        .iter()
                        .any(|span| span.calls.contains(&call))
                        && self.unreliable_only.contains_key(&(position - 1))
                    {
                        let label = LabelResult::PoorSignal {
                            read: read.name.to_owned(),
                        };
                        labels
                            .entry((label.rank(), label.kind().name(), index))
                            .or_insert(label);
                    }
                }
                let trust = observation
                    .call_index_0based
                    .and_then(|call| read.trust.get(call));
                let scope = trust.map(|trust| match trust.scope {
                    TrustScope::Full => ScopeResult::Full,
                    TrustScope::Base => ScopeResult::Base,
                });
                let cause = observation
                    .call_index_0based
                    .and_then(|call| self.second_cause(read, call));
                let second = trust
                    .and_then(|trust| trust.second)
                    .filter(|_| observation.state != LocusState::Unresolved)
                    .map(|second| SecondResult {
                        base: read.strand.reference_base(second.base),
                        share: round(second.share),
                        cause,
                    });
                positions.entry(index).or_default().push(evidence(
                    position,
                    observation,
                    read.sanger,
                    scope,
                    second,
                )?);
            }
        }
        let evidence = positions
            .into_iter()
            .map(|(index, positions)| {
                let read = self.read(index)?;
                Ok(EvidenceReadResult {
                    read: read.name.to_owned(),
                    sha256: read.sha256.to_owned(),
                    strand: read.strand,
                    positions,
                })
            })
            .collect::<Result<_>>()?;
        Ok((labels.into_values().collect(), evidence))
    }

    /// The decision of a site the regions clip, over its reportable columns
    /// `first..=last`: what most of its supporting reads show there (a tie
    /// goes to the first in sequence order); `None` for an unclipped site, an
    /// undecided one, or one no supporting read bounds over those columns.
    fn clipped_decision(&self, site: &ConsensusSite, first: usize, last: usize) -> Option<String> {
        if site.first_0based + 1 == first && site.last_0based + 1 == last {
            return None;
        }
        site.call.as_ref().filter(|call| !call.contains('N'))?;
        let mut shown = BTreeMap::<String, usize>::new();
        for index in site
            .supporting
            .iter()
            .filter_map(|input| self.by_id.get(input.as_str()).copied())
        {
            if let Some(sequence) = self.shown(index, first, last) {
                *shown.entry(sequence).or_default() += 1;
            }
        }
        let most = shown.values().copied().max()?;
        shown
            .into_iter()
            .find(|(_, count)| *count == most)
            .map(|(sequence, _)| sequence)
    }

    /// Keeps a clipped site's dissent label (`strands_disagree` or
    /// `same_strand_dissent`) only for the opposing reads that may differ
    /// from `decision` over its reportable columns: one that shows the
    /// decision there, its difference lying outside the regions, no longer
    /// opposes the item (a read that does not bound those columns still
    /// does). The label is then the one the remaining reads give
    /// ([`dissent`]), or goes when no opposing read is left.
    fn clip_disagreement(
        &self,
        labels: &mut Vec<LabelResult>,
        first: usize,
        last: usize,
        decision: &str,
    ) {
        let agrees = |name: &str| {
            self.reads
                .iter()
                .position(|read| read.name == name)
                .and_then(|index| self.shown(index, first, last))
                .is_some_and(|shown| shown == decision)
        };
        let Some(at) = labels.iter().position(|label| {
            matches!(
                label,
                LabelResult::StrandsDisagree { .. } | LabelResult::SameStrandDissent { .. }
            )
        }) else {
            return;
        };
        let (LabelResult::StrandsDisagree { reads } | LabelResult::SameStrandDissent { reads }) =
            labels.remove(at)
        else {
            return;
        };
        let (supporting, opposing) = reads
            .into_iter()
            .filter(|side| side.role != Role::Opposing || !agrees(&side.read))
            .partition(|side| side.role == Role::Supporting);
        if let Some(label) = dissent(supporting, opposing) {
            labels.insert(at, label);
        }
    }

    /// `read_disagrees` for each read that supports a decided site, or that
    /// the consensus set aside there (`uninformative`), and still shows the
    /// whole item, every call clean of the modality's reasons, as
    /// another sequence than the decision: the consensus explained the
    /// difference away by the read's bounds. Opposing reads are [`dissent`]'s.
    /// Only an unclipped site is compared.
    fn disagreeing(
        &self,
        site: &ConsensusSite,
        first: usize,
        last: usize,
        sequence: &str,
    ) -> Vec<LabelResult> {
        if site.first_0based + 1 != first || site.last_0based + 1 != last {
            return Vec::new();
        }
        site.supporting
            .iter()
            .chain(&site.uninformative)
            .filter_map(|input| self.by_id.get(input.as_str()).copied())
            .filter_map(|index| {
                let shown = self.shown(index, first, last)?;
                (shown != sequence).then(|| LabelResult::ReadDisagrees {
                    read: self.reads[index].name.to_owned(),
                    strand: self.reads[index].strand,
                    sequence: shown,
                })
            })
            .collect()
    }

    /// `read_disagrees` for a read that reads one end of a decided item in
    /// phase as another base than the decision: from the end its strand
    /// reaches first (a forward read's first position, a reverse read's
    /// last), every call a clean base of the read's full trust, up to where
    /// the read has read `review.phase_run` equal bases in a row. The far end
    /// of such a run is a length question that dephasing leaves open; the
    /// bases before it are not. Only an item without a length change is
    /// compared, position by position; an opposing read is
    /// [`dissent`]'s, and a read that [`Self::disagreeing`] already
    /// labels is left to it.
    fn end_disagreeing(
        &self,
        first: usize,
        last: usize,
        site: &ConsensusSite,
        sequence: &str,
        labelled: &[LabelResult],
    ) -> Result<Vec<LabelResult>> {
        let decided = sequence.as_bytes();
        if decided.len() != last - first + 1 {
            return Ok(Vec::new());
        }
        let mut labels = Vec::new();
        let opposing = site
            .opposing
            .iter()
            .filter_map(|input| self.by_id.get(input.as_str()).copied())
            .collect::<BTreeSet<_>>();
        for (index, read) in self.reads.iter().enumerate() {
            if opposing.contains(&index)
                || labelled.iter().any(|label| {
                    matches!(label, LabelResult::ReadDisagrees { read: name, .. } if name == read.name)
                })
            {
                continue;
            }
            let order: Vec<usize> = match read.strand {
                Orientation::Forward => (first..=last).collect(),
                Orientation::Reverse => (first..=last).rev().collect(),
            };
            let bases = order
                .iter()
                .map(|&position| self.clean_base(index, position))
                .collect::<Result<Vec<_>>>()?;
            let mut shown = in_phase_end(bases, self.config.phase_run);
            let differs = order.iter().zip(&shown).any(|(&position, &base)| {
                u8::try_from(base).ok() != decided.get(position - first).copied()
            });
            if differs && shown.len() < decided.len() {
                if read.strand == Orientation::Reverse {
                    shown.reverse();
                }
                labels.push(LabelResult::ReadDisagrees {
                    read: read.name.to_owned(),
                    strand: read.strand,
                    sequence: shown.into_iter().collect(),
                });
            }
        }
        Ok(labels)
    }

    /// The reference-strand base read `index` shows at `position` when its
    /// call there is a clean observation of full trust: neither masked,
    /// unresolved, nor a deletion.
    fn clean_base(&self, index: usize, position: usize) -> Result<Option<char>> {
        let Some(observation) = self
            .loci
            .get(&position)
            .and_then(|locus| locus.observations.iter().find(|o| o.read_index == index))
        else {
            return Ok(None);
        };
        if !matches!(
            observation.state,
            LocusState::Reference | LocusState::Alternate
        ) {
            return Ok(None);
        }
        let (Some(base), Some(call)) = (observation.base, observation.call_index_0based) else {
            return Ok(None);
        };
        let trust = self
            .read(index)?
            .trust
            .get(call)
            .ok_or(ReportError::Inconsistent(
                "review observation lacks the trust of its call",
            ))?;
        Ok((trust.scope == TrustScope::Full).then_some(base))
    }

    /// The sequence one read shows over `first..=last`; see [`bounded`].
    fn shown(&self, index: usize, first: usize, last: usize) -> Option<String> {
        let read = self.reads.get(index)?;
        bounded(read.columns, &read.trust, first, last)
    }

    /// Why call `call` of `read` shows a second peak: the modality's cause
    /// (`SangerCallSignal::second`), refined with the review's run signal,
    /// the sample's own reads taken out, so that a share its key's run level
    /// holds reads `context`. `None` when the call's second share does not
    /// reach `signal_processing.second_floor` and it is not masked as mixed.
    fn second_cause(&self, read: &Read<'_>, call: usize) -> Option<&'static str> {
        let peak = read.sanger.calls.get(call)?.second.as_ref()?;
        let key = Key::of(read.evidence, call, peak.base, self.second.run_signal.flank);
        Some(run_signal::refine(peak, key.as_ref(), self.signal.as_ref(), self.second).label())
    }

    /// The second base of a call at `position`, on the reference strand, its
    /// peak share (SRS-CALL-015), its threshold (SRS-CALL-020), and the bases
    /// the same read shows at the neighbouring positions.
    fn judged(
        &self,
        position: usize,
        read: usize,
        call: usize,
        trust: &CallTrust,
    ) -> Option<Clean> {
        let read_index = read;
        let read = self.reads.get(read)?;
        let second = trust.second?;
        let key = Key::of(
            read.evidence,
            call,
            second.base,
            self.second.run_signal.flank,
        );
        let neighbour = |other: Option<usize>| {
            other
                .and_then(|other| self.observation(read_index, other))
                .and_then(|other| other.base)
        };
        Some(Clean {
            read: read_index,
            strand: read.strand,
            second: read.strand.reference_base(second.base),
            share: second.share,
            threshold: run_signal::threshold(key.as_ref(), self.signal.as_ref(), self.second),
            neighbours: [
                neighbour(position.checked_sub(1)),
                neighbour(position.checked_add(1)),
            ],
            explained: read
                .evidence
                .calls()
                .get(call)
                .is_some_and(|evidence| !evidence.second_unexplained),
        })
    }

    /// `reference_second` at each position of `first..=last` where read
    /// `index`, one of the reads of the one strand a difference rests on,
    /// reads another base than the reference while its call shows the
    /// reference base as a second peak above the call's threshold
    /// (SRS-CALL-020). One strand cannot tell that base from a broad peak of
    /// another dye over the reference's, nor from a mixture, and its reruns
    /// repeat a dye blob at its primer's start; on development panel v5 most
    /// such bases at reverse read starts were dye blobs over the reference
    /// (SRS-CALL-021).
    fn reference_second(
        &self,
        index: usize,
        first: usize,
        last: usize,
    ) -> Result<Vec<LabelResult>> {
        let read = self.read(index)?;
        let mut labels = Vec::new();
        for position in first..=last {
            let Some(observation) = self.observation(index, position) else {
                continue;
            };
            let reference = self.base(position);
            if observation.state != LocusState::Alternate
                || observation.base.is_none_or(|base| base == reference)
            {
                continue;
            }
            let Some(call) = observation.call_index_0based else {
                continue;
            };
            let trust = read.trust.get(call).ok_or(ReportError::Inconsistent(
                "review observation lacks the trust of its call",
            ))?;
            let Some(judged) = self.judged(position, index, call, trust) else {
                continue;
            };
            if reference_second(&judged, reference) {
                labels.push(LabelResult::ReferenceSecond {
                    position,
                    read: read.name.to_owned(),
                    strand: read.strand,
                    share: round(judged.share),
                    threshold: round(judged.threshold.share),
                    calls: judged.threshold.calls,
                });
            }
        }
        Ok(labels)
    }

    /// `possible_heteroplasmy` at one position, when it fires.
    fn heteroplasmy(&self, position: usize) -> Option<LabelResult> {
        let locus = self.loci.get(&position)?;
        let clean = locus
            .observations
            .iter()
            .filter_map(|observation| {
                if !matches!(
                    observation.state,
                    LocusState::Reference | LocusState::Alternate
                ) {
                    return None;
                }
                let read = self.reads.get(observation.read_index)?;
                let trust = read.trust.get(observation.call_index_0based?)?;
                if trust.placement != Placement::Scored || !(trust.trusted || trust.exempt) {
                    return None;
                }
                self.judged(
                    position,
                    observation.read_index,
                    observation.call_index_0based?,
                    trust,
                )
            })
            .collect::<Vec<_>>();
        let (base, reads) = heteroplasmy(&clean)?;
        Some(LabelResult::PossibleHeteroplasmy {
            position,
            base,
            explained: explained(&reads),
            reads: reads
                .into_iter()
                .map(|clean| ShareResult {
                    read: self.reads[clean.read].name.to_owned(),
                    strand: clean.strand,
                    share: round(clean.share),
                    threshold: round(clean.threshold.share),
                    calls: clean.threshold.calls,
                })
                .collect(),
        })
    }

    /// The first reference position of the masked stretch of read `index`
    /// that reaches `position`, walking back in the read's own direction.
    fn shift_start(&self, index: usize, position: usize) -> usize {
        let back = |position: usize| match self.reads[index].strand {
            Orientation::Forward => position.checked_sub(1),
            Orientation::Reverse => position.checked_add(1),
        };
        let mut from = position;
        while let Some(previous) = back(from) {
            if self
                .observation(index, previous)
                .is_some_and(|observation| observation.state == LocusState::Masked)
            {
                from = previous;
            } else {
                break;
            }
        }
        from
    }

    fn observation(&self, index: usize, position: usize) -> Option<&SampleLocusObservation> {
        self.loci
            .get(&position)?
            .observations
            .iter()
            .find(|observation| observation.read_index == index)
    }

    fn read(&self, index: usize) -> Result<&Read<'a>> {
        Ok(self
            .reads
            .get(index)
            .ok_or(ReportError::MissingRead { index })?)
    }

    fn base(&self, position: usize) -> char {
        char::from(self.sequence[position - 1])
    }
}

/// One read's evidence at one position, with its Sanger scan position and
/// relative quality.
/// One noted run of a length tract: its judgement when a read reaches it,
/// else an unreported run no read reaches. Where the profile records no
/// length (`recorded` false), no reading is given.
fn tract_run(
    run: usize,
    window: &length::LengthWindow,
    called: Option<&[usize]>,
    reported: bool,
    judged: Option<&length::Judged>,
) -> TractRunResult {
    let Some(judged) = judged else {
        return TractRunResult {
            run,
            base: window.repeat,
            called: called
                .filter(|_| reported)
                .and_then(|runs| runs.get(run - 1))
                .copied(),
            status: "unconfirmed",
            reading: None,
            agreement: None,
            strands: StrandsResult::default(),
            in_phase: Vec::new(),
            floor: None,
            convention: None,
            reads: Vec::new(),
        };
    };
    let result = &judged.result;
    let (status, reading) = length::status(result).unwrap_or(("unconfirmed", result.floor));
    let recorded = window.recorded;
    TractRunResult {
        run,
        base: result.base,
        called: result.called,
        status,
        reading: reading.filter(|_| recorded),
        agreement: Some(result.agreement),
        strands: judged.strands,
        in_phase: result
            .reads
            .iter()
            .filter(|read| read.source == "in_phase")
            .map(|read| read.read.clone())
            .collect(),
        floor: result.floor.filter(|_| recorded),
        convention: result.convention.clone().filter(|_| recorded),
        reads: judged.reads.clone(),
    }
}

fn evidence(
    position: usize,
    observation: &SampleLocusObservation,
    sanger: &SangerRead,
    scope: Option<ScopeResult>,
    second: Option<SecondResult>,
) -> Result<EvidencePositionResult> {
    let base = match observation.state {
        LocusState::Deletion => '-',
        LocusState::Unresolved => 'N',
        _ => observation.base.unwrap_or('N'),
    };
    let (scan, quality) = match observation.call_index_0based {
        Some(index) => {
            let Some(call) = sanger.calls.get(index) else {
                return Err(ReportError::Inconsistent(
                    "review evidence call lacks its Sanger record",
                )
                .into());
            };
            (Some(call.scan), Some(call.quality))
        }
        None => (None, None),
    };
    Ok(EvidencePositionResult {
        position,
        ordinal: 0,
        base,
        state: observation.state,
        call: observation.call_index_0based,
        scope,
        scan,
        quality,
        second,
    })
}

/// One read that reads a position cleanly: trusted or exempt, unmasked, with
/// a second base.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Clean {
    /// Read index.
    read: usize,
    strand: Orientation,
    /// The call's second base on the reference strand, and its peak share.
    second: char,
    share: f64,
    /// The share the second base must exceed.
    threshold: Threshold,
    /// The bases the same read shows at the previous and the next reference
    /// position, when it shows one.
    neighbours: [Option<char>; 2],
    /// Whether the modality explains the call's second peak by the signal's
    /// own shape (SRS-CALL-022: not `second_unexplained`).
    explained: bool,
}

/// Whether a call reading another base than `reference` shows the
/// reference base as a second peak above the call's threshold
/// (SRS-CALL-020): `reference_second` when the difference rests on that
/// read's strand.
fn reference_second(call: &Clean, reference: char) -> bool {
    call.second == reference && call.share > call.threshold.share
}

/// Whether a `possible_heteroplasmy`'s qualifying reads leave one strand
/// without a second peak the modality leaves unexplained: every qualifying
/// read of that strand shows a neighbour's tail or a dye blob, not a possible
/// second sequence. The label is then context, not a check (on the
/// development panel, none of 23 such flags was a heteroplasmy).
fn explained(qualifying: &[Clean]) -> bool {
    [Orientation::Forward, Orientation::Reverse]
        .into_iter()
        .any(|strand| {
            qualifying
                .iter()
                .filter(|read| read.strand == strand)
                .all(|read| read.explained)
        })
}

/// `possible_heteroplasmy` over the clean reads of one position: a read
/// qualifies when its second peak share exceeds its threshold (SRS-CALL-020)
/// and its second base is not a base it shows next to the position (a
/// neighbour shadow). When both strands
/// read the position cleanly, a qualifying read on each strand must show the
/// same second base; when one strand does, one qualifying read suffices. The
/// base most qualifying reads show wins, then the first in A/C/G/T order. It
/// returns that base and its qualifying reads, in input order.
fn heteroplasmy(clean: &[Clean]) -> Option<(char, Vec<Clean>)> {
    let qualifying = clean
        .iter()
        .filter(|read| {
            read.share > read.threshold.share && !read.neighbours.contains(&Some(read.second))
        })
        .collect::<Vec<_>>();
    let strands = clean.iter().map(|read| read.strand).collect::<Vec<_>>();
    let both = strands.contains(&Orientation::Forward) && strands.contains(&Orientation::Reverse);
    let shows = |base: char, strand: Orientation| {
        qualifying
            .iter()
            .any(|read| read.second == base && read.strand == strand)
    };
    let mut best: Option<(usize, char)> = None;
    for base in Nucleotide::ALL.map(Nucleotide::as_char) {
        let count = qualifying.iter().filter(|read| read.second == base).count();
        // Both strands must confirm the base: a site only one strand reads
        // cleanly gives no label (the one-strand flags were 3 % precise).
        let fires = both && shows(base, Orientation::Forward) && shows(base, Orientation::Reverse);
        if fires && best.is_none_or(|(most, _)| count > most) {
            best = Some((count, base));
        }
    }
    let (_, base) = best?;
    Some((
        base,
        qualifying
            .into_iter()
            .filter(|read| read.second == base)
            .copied()
            .collect(),
    ))
}

/// The label of a site's clean dissent, over every supporting and opposing
/// read, however many each strand has: `strands_disagree` when a supporting
/// read lies on another strand than an opposing one, otherwise
/// `same_strand_dissent` (all of them lie on one strand); `None` without a
/// supporting or an opposing read. Lists the supporting, then the opposing
/// reads.
fn dissent(supporting: Vec<SideResult>, opposing: Vec<SideResult>) -> Option<LabelResult> {
    if supporting.is_empty() || opposing.is_empty() {
        return None;
    }
    let across = supporting
        .iter()
        .any(|left| opposing.iter().any(|right| left.strand != right.strand));
    let reads = supporting.into_iter().chain(opposing).collect();
    Some(if across {
        LabelResult::StrandsDisagree { reads }
    } else {
        LabelResult::SameStrandDissent { reads }
    })
}

/// The one strand every supporting read lies on; `None` without a read or
/// with reads on both strands.
fn one_strand(supporting: &[Orientation]) -> Option<Orientation> {
    let (&first, rest) = supporting.split_first()?;
    rest.iter().all(|&strand| strand == first).then_some(first)
}

/// The runs of a decided site longer than their reference run by `growth`
/// or more whose length rests on an anchored end or on one strand
/// (`one_strand`: every supporting read lies on one strand, however many):
/// the reference and the call must hold the same runs in the same order. A
/// length both strands show in phase is assured, however long.
fn long_runs(site: &ConsensusSite, one_strand: bool, growth: usize) -> Vec<GrowthResult> {
    let Some(call) = site.call.as_deref() else {
        return Vec::new();
    };
    let reference = runs_of(site.reference.chars());
    let called = runs_of(call.chars());
    if reference.len() != called.len()
        || reference
            .iter()
            .zip(&called)
            .any(|(left, right)| left.0 != right.0)
    {
        return Vec::new();
    }
    // The consensus lists a run's evidence when it decided the site run by
    // run, in call order.
    let evidence = (site.runs.len() == called.len()).then_some(&site.runs);
    reference
        .iter()
        .zip(&called)
        .enumerate()
        .filter(|(_, (base_run, call_run))| call_run.1 >= base_run.1 + growth)
        .filter_map(|(index, (base_run, call_run))| {
            let run = evidence.and_then(|runs| runs.get(index));
            let anchored = run.is_some_and(|run| run.evidence == RunEvidence::AnchoredEnd);
            (anchored || one_strand).then(|| GrowthResult {
                base: call_run.0,
                length: call_run.1,
                reference_length: base_run.1,
                length_evidence: run.map(|run| run.evidence.label()),
            })
        })
        .collect()
}

/// Whether an item over `first..=last` lies within `zone` positions of a
/// reportable region's end that the call leaves out as not analysed
/// (SRS-CONS-012): the reads there sit right after their poor start, where
/// indels are often wrong.
fn beside_not_analysed(
    regions: &[[usize; 2]],
    not_analysed: &[[usize; 2]],
    first: usize,
    last: usize,
    zone: usize,
) -> bool {
    regions.iter().any(|&[start, end]| {
        start <= first
            && last <= end
            && not_analysed.iter().any(|&[cut_first, cut_last]| {
                (cut_first == start && first < start + zone)
                    || (cut_last == end && last + zone > end)
            })
    })
}

/// Whether `first..=last` (1-based) lies in a region and within `zone`
/// positions of its first or last position.
fn near_edge(regions: &[[usize; 2]], first: usize, last: usize, zone: usize) -> bool {
    regions.iter().any(|&[start, end]| {
        start <= first && last <= end && (first - start < zone || end - last < zone)
    })
}

/// Whether `call` changes only the length of runs of `reference` (the same
/// bases in the same run order, another length), one of them a reference
/// run of at least `run` equal bases.
fn run_length_change(reference: &str, call: &str, run: usize) -> bool {
    let (from, to) = (runs_of(reference.chars()), runs_of(call.chars()));
    call.len() != reference.len()
        && from.len() == to.len()
        && from.iter().zip(&to).all(|(left, right)| left.0 == right.0)
        && from
            .iter()
            .zip(&to)
            .any(|(left, right)| left.1 != right.1 && left.1 >= run)
}

/// Whether `call` is at least `growth` bases shorter than `reference` and the
/// reference holds a run of at least `run` equal bases.
fn run_lost(reference: &str, call: &str, growth: usize, run: usize) -> bool {
    call.len() + growth <= reference.len()
        && runs_of(reference.chars())
            .iter()
            .any(|(_, length)| *length >= run)
}

/// Whether the consensus left a site undecided: no clean observation, a
/// tie (even one the consensus settled on the reference, whichever strands
/// the tied reads are on), or an `N` in its call.
fn undecided(site: &ConsensusSite) -> bool {
    matches!(site.state, SiteState::Unresolved | SiteState::Contested)
        || site.call.as_ref().is_some_and(|call| call.contains('N'))
}

/// Whether DNA calls anything within one position of `first..=last`: a
/// call the notation moved next to the item still answers its difference.
fn called_near(calls: &[NotationCall], first: usize, last: usize) -> bool {
    calls
        .iter()
        .any(|call| first <= call.position + 1 && call.position <= last + 1)
}

/// Whether a report convention leaves `run`'s length unrecorded at `site`: a
/// run-length convention whose window the site overlaps, of the run's base,
/// whose `minimum_run` the run reaches without the bases the site adds to the
/// reference. An estimate can only lengthen a run past the convention's
/// bound by adding bases, so a run that reaches it only by them may be one
/// whose length the reviewers do record. Otherwise the estimate leaves the
/// reviewer nothing to record, and `estimated_length` does not list it.
fn not_recorded(conventions: &[Convention], site: &ConsensusSite, run: &ConsensusRun) -> bool {
    let added = site
        .call
        .as_ref()
        .map_or(0, |call| call.len().saturating_sub(site.reference.len()));
    conventions.iter().any(|convention| {
        let [first, last] = convention.window;
        let ConventionKind::RunLength { base, minimum_run } = convention.kind;
        convention.label == ConventionLabel::NotRecorded
            && site.first_0based < last
            && first <= site.last_0based + 1
            && u8::try_from(run.base).ok() == Some(base)
            && run.length.saturating_sub(added) >= minimum_run
    })
}

/// The bases a read shows from the end of an item it reaches first, in its
/// reading order: up to its first unclean position (`None`), and short of the
/// `phase_run`-th equal base in a row, where the run's far end begins.
fn in_phase_end(bases: impl IntoIterator<Item = Option<char>>, phase_run: usize) -> Vec<char> {
    let mut shown = Vec::new();
    let mut run = 0;
    for base in bases {
        let Some(base) = base else { break };
        run = if shown.last() == Some(&base) {
            run + 1
        } else {
            1
        };
        if run >= phase_run {
            break;
        }
        shown.push(base);
    }
    shown
}

/// Whether a call shows its base for `read_disagrees`: no modality reason
/// against its base, and anchoring, or scored clear of the read's edges.
fn shows(trust: &CallTrust) -> bool {
    trust.vetoes.is_empty()
        && match trust.placement {
            Placement::Scored => trust.trusted || trust.exempt,
            Placement::Anchoring(_) => true,
            Placement::Unresolved(_) => false,
        }
}

/// The sequence a read's alignment `columns` show over the 1-based
/// `first..=last` within one linear span of the reference, with every
/// inserted call after a position of the item, masked ones included; `None`
/// unless the read aligns every position of the item and reads the positions
/// on both sides of it, each call a base that [`shows`] it by its `trust`,
/// so that the read bounds what it shows.
fn bounded(
    columns: &[AlignmentColumn],
    trust: &[CallTrust],
    first: usize,
    last: usize,
) -> Option<String> {
    let shows_call = |column: &AlignmentColumn| {
        is_canonical(column.query_base())
            && column
                .original_call_index_0based()
                .and_then(|call| trust.get(call))
                .is_some_and(shows)
    };
    // 0-based reference indices of the flanks.
    let (before, after) = (first.checked_sub(2)?, last);
    let mut shown = String::new();
    let mut covered = 0;
    let mut bounds = [false, false];
    let mut within = false;
    for column in columns {
        match column.reference_index_0based() {
            Some(reference) if reference == before || reference == after => {
                within = false;
                bounds[usize::from(reference == after)] = shows_call(column);
                continue;
            }
            Some(reference) => {
                within = first <= reference + 1 && reference < last;
                if !within {
                    continue;
                }
                covered += 1;
            }
            None if !within => continue,
            None => {}
        }
        if column.query_base() == '-' {
            continue;
        }
        if !shows_call(column) {
            return None;
        }
        shown.push(column.query_base());
    }
    (covered == last - first + 1 && bounds == [true, true]).then_some(shown)
}

/// The reference positions, within `regions`, that a read's aligned calls
/// cover inside one of its unreliable spans (SRS-CALL-024), as maximal
/// stretches of one reason.
fn unreliable_spans(
    columns: &[AlignmentColumn],
    evidence: &ReadEvidence,
    regions: &[[usize; 2]],
) -> Vec<UnreliableSpanResult> {
    let mut spans: Vec<UnreliableSpanResult> = Vec::new();
    for column in columns {
        let (Some(reference), Some(call)) = (
            column.reference_index_0based(),
            column.original_call_index_0based(),
        ) else {
            continue;
        };
        let position = reference + 1;
        let Some(span) = evidence
            .unreliable()
            .iter()
            .find(|span| span.calls.contains(&call))
        else {
            continue;
        };
        if !inside(regions, position) {
            continue;
        }
        let reason = span.reason.label();
        match spans.last_mut() {
            Some(last) if last.last + 1 == position && last.reason == reason => {
                last.last = position;
            }
            _ => spans.push(UnreliableSpanResult {
                first: position,
                last: position,
                reason,
            }),
        }
    }
    spans
}

/// The reference positions, within `regions`, that a read's aligned calls
/// cover while the modality masks its signal and the call stays masked
/// (anchoring or unresolved, not vouched for), as maximal stretches of one
/// phase state: where the read does not show its bases (it reads past a
/// repeat out of phase, its signal is mixed, weak, or irregularly spaced).
/// Recorded so that what the read did not show can be located and read later.
fn masked_spans(
    columns: &[AlignmentColumn],
    segments: &[CallabilitySegmentResult],
    trust: &[CallTrust],
    regions: &[[usize; 2]],
) -> Vec<MaskedSpanResult> {
    let mut spans: Vec<MaskedSpanResult> = Vec::new();
    for column in columns {
        let (Some(reference), Some(call)) = (
            column.reference_index_0based(),
            column.original_call_index_0based(),
        ) else {
            continue;
        };
        let position = reference + 1;
        if !inside(regions, position) {
            continue;
        }
        let Some(segment) = segments
            .iter()
            .find(|segment| segment.calls.start <= call && call < segment.calls.end)
        else {
            continue;
        };
        // A call the modality vouched for despite its segment (a leading
        // call, SRS-CALL-014) is scored: the read shows that base.
        if segment.state == PhaseState::InPhase
            || trust
                .get(call)
                .is_none_or(|trust| trust.placement == Placement::Scored)
        {
            continue;
        }
        match spans.last_mut() {
            Some(span)
                if span.last + 1 == position
                    && span.state == segment.state
                    && span.after_repeat == segment.after_repeat =>
            {
                span.last = position;
            }
            _ => spans.push(MaskedSpanResult {
                first: position,
                last: position,
                state: segment.state,
                after_repeat: segment.after_repeat,
            }),
        }
    }
    spans
}

/// Label order within an item: rank, name, then the position or read.
fn key(label: &LabelResult) -> (u8, &'static str, usize, &str) {
    let (position, read) = match label {
        LabelResult::PossibleHeteroplasmy { position, .. } => (*position, ""),
        LabelResult::ReferenceSecond { position, read, .. } => (*position, read.as_str()),
        LabelResult::EdgeDeletion { reads, .. }
        | LabelResult::RunDeletion { reads, .. }
        | LabelResult::RunIndel { reads, .. } => (0, reads.first().map_or("", String::as_str)),
        LabelResult::OneRead { read }
        | LabelResult::ReadDisagrees { read, .. }
        | LabelResult::NoisySignal { read, .. }
        | LabelResult::PhaseShift { read, .. }
        | LabelResult::ReadEdge { read }
        | LabelResult::PoorSignal { read } => (0, read.as_str()),
        LabelResult::StrandsDisagree { .. }
        | LabelResult::SameStrandDissent { .. }
        | LabelResult::NotCalled
        | LabelResult::LongRun { .. }
        | LabelResult::OneStrandInsertion { .. }
        | LabelResult::EstimatedLength { .. }
        | LabelResult::LengthUnconfirmed { .. }
        | LabelResult::LengthNotRecorded { .. }
        | LabelResult::TractLength { .. }
        | LabelResult::Undecided
        | LabelResult::OneStrand { .. }
        | LabelResult::NotCovered
        | LabelResult::DephasedEnd { .. }
        | LabelResult::BesideNotAnalysed => (0, ""),
    };
    (label.rank(), label.kind().name(), position, read)
}

/// An item's priority: its best label rank, or [`UNLABELLED`].
fn priority(labels: &[LabelResult]) -> u8 {
    labels
        .iter()
        .map(LabelResult::rank)
        .min()
        .unwrap_or(UNLABELLED)
}

/// The reportable positions within `regions` as zones a reviewer opens
/// whole. Each position takes the first kind that holds, in [`ZoneKind`]
/// order: inside a repeat window of the profile (1-based inclusive, with its
/// name), outside the `callable` ranges, inside a read's masked stretch past a
/// repeat, inside another masked stretch, else clean. A zone is a maximal
/// stretch of one kind (and window) within one region; it counts the `check`
/// items touching it. Every boundary comes from the profile and the sample's
/// reads, none from a width.
fn zones(
    regions: &[[usize; 2]],
    windows: &[(usize, usize, &str)],
    callable: &[[usize; 2]],
    reads: &[ReviewReadResult],
    items: &[ItemResult],
) -> Vec<ZoneResult> {
    let masked = reads
        .iter()
        .flat_map(|read| read.masked.iter().flatten())
        .collect::<Vec<_>>();
    let masks = |position: usize, after_repeat: bool| {
        masked.iter().any(|span| {
            (!after_repeat || span.after_repeat) && span.first <= position && position <= span.last
        })
    };
    let mut zones: Vec<ZoneResult> = Vec::new();
    for &[start, end] in regions {
        for position in start..=end {
            let window = windows
                .iter()
                .find(|(first, last, _)| *first <= position && position <= *last);
            let (kind, name) = match window {
                Some(&(_, _, name)) => (ZoneKind::Repeat, Some(name)),
                None if !inside(callable, position) => (ZoneKind::NotCallable, None),
                None if masks(position, true) => (ZoneKind::Dephased, None),
                None if masks(position, false) => (ZoneKind::Masked, None),
                None => (ZoneKind::Clean, None),
            };
            match zones.last_mut() {
                Some(zone)
                    if position != start && zone.kind == kind && zone.name.as_deref() == name =>
                {
                    zone.last = position;
                }
                _ => zones.push(ZoneResult {
                    kind,
                    first: position,
                    last: position,
                    name: name.map(str::to_owned),
                    checks: 0,
                }),
            }
        }
    }
    for zone in &mut zones {
        zone.checks = items
            .iter()
            .filter(|item| {
                item.kind == ItemKind::Check
                    && item.first.position <= zone.last
                    && zone.first <= item.last.position
            })
            .count();
    }
    zones
}

/// Counts over the items, the reportable positions inside `regions` that
/// `check` items hold, the `reference_only` positions, and the `zones` the
/// checks touch; `clear` when the call succeeded and no item is a check.
fn summarize(
    items: &[ItemResult],
    status: ReviewStatus,
    regions: &[[usize; 2]],
    zones: &[ZoneResult],
    reference_only: &[usize],
    not_analysed: usize,
) -> SummaryResult {
    let mut labels = LabelCounts::default();
    for item in items {
        let kinds = item
            .labels
            .iter()
            .map(LabelResult::kind)
            .collect::<BTreeSet<_>>();
        for kind in kinds {
            labels.add(kind);
        }
    }
    let checks = items
        .iter()
        .filter(|item| item.kind == ItemKind::Check)
        .count();
    let check = items
        .iter()
        .filter(|item| item.kind == ItemKind::Check)
        .flat_map(|item| item.first.position..=item.last.position)
        .filter(|&position| inside(regions, position))
        .collect::<BTreeSet<_>>()
        .len();
    SummaryResult {
        items: items.len(),
        checks,
        variants: items
            .iter()
            .filter(|item| item.kind == ItemKind::Variant)
            .count(),
        labels,
        positions: PositionCounts {
            reportable: regions.iter().map(|[start, end]| end - start + 1).sum(),
            check,
            reference_only: reference_only.len(),
            not_analysed,
        },
        zones: ZoneCounts {
            total: zones.len(),
            check: zones.iter().filter(|zone| zone.checks > 0).count(),
        },
        clear: status == ReviewStatus::Called && checks == 0,
    }
}

/// A 0-based half-open reference segment as a pair.
const fn span(segment: &ReferenceSegment) -> (usize, usize) {
    (segment.start_0based, segment.end_0based_exclusive)
}

/// 0-based half-open segments as 1-based inclusive ranges, empty ones left
/// out.
fn segments(segments: impl IntoIterator<Item = (usize, usize)>) -> Vec<[usize; 2]> {
    segments
        .into_iter()
        .filter(|(start, end)| start < end)
        .map(|(start, end)| [start + 1, end])
        .collect()
}

/// 1-based inclusive ranges sorted and merged where they overlap or touch.
fn merge(mut ranges: Vec<[usize; 2]>) -> Vec<[usize; 2]> {
    ranges.sort_unstable();
    let mut merged: Vec<[usize; 2]> = Vec::with_capacity(ranges.len());
    for [start, end] in ranges {
        match merged.last_mut() {
            Some(last) if start <= last[1] + 1 => last[1] = last[1].max(end),
            _ => merged.push([start, end]),
        }
    }
    merged
}

/// The profile's regions within a reference of `length` bases, merged.
fn clip(regions: &[[usize; 2]], length: usize) -> Vec<[usize; 2]> {
    merge(
        regions
            .iter()
            .map(|&[start, end]| [start.max(1), end.min(length)])
            .filter(|[start, end]| start <= end)
            .collect(),
    )
}

/// Whether merged `ranges` hold `position`.
fn inside(ranges: &[[usize; 2]], position: usize) -> bool {
    ranges
        .iter()
        .any(|&[start, end]| start <= position && position <= end)
}

/// The first and last positions of `first..=last` inside merged `regions`.
fn within(regions: &[[usize; 2]], first: usize, last: usize) -> Option<[usize; 2]> {
    let clipped = regions
        .iter()
        .map(|&[start, end]| [start.max(first), end.min(last)])
        .filter(|[start, end]| start <= end)
        .collect::<Vec<_>>();
    Some([clipped.first()?[0], clipped.last()?[1]])
}

/// A share rounded to three decimals.
pub(super) fn round(share: f64) -> f64 {
    (share * 1000.0).round() / 1000.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::review_result::LabelKind;
    use dna_kernel::read_evidence::EvidenceReason;

    const F: Orientation = Orientation::Forward;
    const R: Orientation = Orientation::Reverse;

    fn clean(read: usize, strand: Orientation, second: char, share: f64) -> Clean {
        Clean {
            read,
            strand,
            second,
            share,
            threshold: Threshold {
                share: 0.12,
                calls: None,
            },
            neighbours: [None, None],
            explained: false,
        }
    }

    #[test]
    fn heteroplasmy_on_both_strands_needs_the_same_base_on_each() {
        let reads = [clean(0, F, 'T', 0.2), clean(1, R, 'T', 0.15)];
        let (base, qualifying) = heteroplasmy(&reads).unwrap_or(('-', Vec::new()));
        assert_eq!(base, 'T');
        assert_eq!(qualifying, reads);

        let other = [clean(0, F, 'T', 0.2), clean(1, R, 'G', 0.3)];
        assert_eq!(heteroplasmy(&other), None);
        let one_side = [clean(0, F, 'T', 0.2), clean(1, R, 'T', 0.05)];
        assert_eq!(heteroplasmy(&one_side), None);
    }

    #[test]
    fn heteroplasmy_needs_both_strands_to_confirm() {
        // One strand alone never fires, however strong its share.
        assert_eq!(heteroplasmy(&[clean(0, R, 'A', 0.4)]), None);
        let one_strand = [clean(0, F, 'T', 0.3), clean(1, F, 'T', 0.3)];
        assert_eq!(heteroplasmy(&one_strand), None);
        // Both strands with the same second base fire, and the reads that
        // show it are reported.
        let both = [
            clean(0, F, 'T', 0.3),
            clean(1, R, 'T', 0.25),
            clean(2, R, 'C', 0.01),
        ];
        let (base, qualifying) = heteroplasmy(&both).unwrap_or(('-', Vec::new()));
        assert_eq!(base, 'T');
        assert_eq!(qualifying, [both[0], both[1]]);
        // The base most reads show wins, then A/C/G/T order.
        let tied = [
            clean(0, F, 'T', 0.2),
            clean(1, R, 'T', 0.2),
            clean(2, F, 'C', 0.2),
            clean(3, R, 'C', 0.2),
            clean(4, R, 'T', 0.2),
        ];
        assert_eq!(heteroplasmy(&tied).map(|(base, _)| base), Some('T'));
        let even = [
            clean(0, F, 'T', 0.2),
            clean(1, R, 'T', 0.2),
            clean(2, F, 'C', 0.2),
            clean(3, R, 'C', 0.2),
        ];
        assert_eq!(heteroplasmy(&even).map(|(base, _)| base), Some('C'));
    }

    #[test]
    fn a_reference_second_needs_the_reference_base_above_the_threshold() {
        // One reverse read reads A where the reference is C, a C under it.
        assert!(reference_second(&clean(0, R, 'C', 0.39), 'C'));
        // Another second base, or one at the threshold, is no such label;
        // neighbours do not matter: a blob covers them too.
        assert!(!reference_second(&clean(0, R, 'G', 0.39), 'C'));
        assert!(!reference_second(&clean(0, R, 'C', 0.12), 'C'));
        let mut beside = clean(0, F, 'C', 0.4);
        beside.neighbours = [Some('C'), None];
        assert!(reference_second(&beside, 'C'));
        // A run level above the share holds it back.
        let mut raised = clean(0, R, 'C', 0.3);
        raised.threshold = Threshold {
            share: 0.35,
            calls: Some(40),
        };
        assert!(!reference_second(&raised, 'C'));
    }

    #[test]
    fn heteroplasmy_needs_a_share_above_the_threshold() {
        let pair = |share| [clean(0, F, 'G', share), clean(1, R, 'G', share)];
        assert_eq!(heteroplasmy(&pair(0.12)), None);
        assert!(heteroplasmy(&pair(0.121)).is_some());
        // A run level above the share on one strand holds the site back.
        let mut raised = pair(0.3);
        raised[0].threshold = Threshold {
            share: 0.35,
            calls: Some(40),
        };
        assert_eq!(heteroplasmy(&raised), None);
        assert_eq!(heteroplasmy(&[]), None);
    }

    #[test]
    fn heteroplasmy_excludes_a_neighbour_shadow() {
        let mut pair = [clean(0, F, 'G', 0.4), clean(1, R, 'G', 0.4)];
        pair[0].neighbours = [Some('A'), Some('G')];
        assert_eq!(heteroplasmy(&pair), None);
        pair[0].neighbours = [Some('G'), None];
        assert_eq!(heteroplasmy(&pair), None);
        pair[0].neighbours = [Some('A'), Some('C')];
        assert!(heteroplasmy(&pair).is_some());
    }

    #[test]
    fn a_strand_whose_second_peaks_the_modality_explains_makes_heteroplasmy_context() {
        let shaped = |mut read: Clean| {
            read.explained = true;
            read
        };
        // Centred, unexplained second peaks on both strands: a check.
        assert!(!explained(&[clean(0, F, 'G', 0.3), clean(1, R, 'G', 0.3)]));
        // One strand's only qualifying read is a neighbour's tail or a blob.
        assert!(explained(&[
            clean(0, F, 'G', 0.3),
            shaped(clean(1, R, 'G', 0.3))
        ]));
        // Another read on that strand still shows an unexplained one.
        assert!(!explained(&[
            clean(0, F, 'G', 0.3),
            shaped(clean(1, R, 'G', 0.3)),
            clean(2, R, 'G', 0.25),
        ]));
        let label = |explained| LabelResult::PossibleHeteroplasmy {
            position: 1,
            base: 'G',
            explained,
            reads: Vec::new(),
        };
        assert!(label(false).unassured() && !label(true).unassured());
        assert_eq!((label(false).rank(), label(true).rank()), (1, 4));
    }

    fn site(
        reference: &str,
        call: &str,
        supporting: usize,
        runs: Vec<ConsensusRun>,
    ) -> ConsensusSite {
        ConsensusSite {
            first_0based: 302,
            last_0based: 302 + reference.len() - 1,
            reference: reference.to_owned(),
            call: Some(call.to_owned()),
            state: SiteState::Called,
            supporting: (0..supporting).map(|read| format!("r{read}")).collect(),
            opposing: Vec::new(),
            uninformative: Vec::new(),
            runs,
        }
    }

    fn run(base: char, length: usize, evidence: RunEvidence) -> ConsensusRun {
        ConsensusRun {
            base,
            length,
            evidence,
        }
    }

    #[test]
    fn a_long_run_needs_both_strands_in_phase() {
        // C7 T C5 read as C10 T C5: three more C than the reference run.
        let reference = "CCCCCCCTCCCCC";
        let call = "CCCCCCCCCCTCCCCC";
        let anchored = vec![
            run('C', 10, RunEvidence::AnchoredEnd),
            run('T', 1, RunEvidence::InPhase),
            run('C', 5, RunEvidence::InPhase),
        ];
        let flagged = long_runs(&site(reference, call, 2, anchored), false, 3);
        assert_eq!(
            flagged,
            [GrowthResult {
                base: 'C',
                length: 10,
                reference_length: 7,
                length_evidence: Some("anchored_end"),
            }]
        );
        // Shown in phase on both strands: assured.
        let in_phase = vec![
            run('C', 10, RunEvidence::InPhase),
            run('T', 1, RunEvidence::InPhase),
            run('C', 5, RunEvidence::InPhase),
        ];
        assert!(long_runs(&site(reference, call, 2, in_phase.clone()), false, 3).is_empty());
        // Resting on one strand: a check, however it is known and however
        // many reads of the strand show it (ADR-0090).
        for reads in 1..=3 {
            assert_eq!(
                long_runs(&site(reference, call, reads, in_phase.clone()), true, 3).len(),
                1
            );
        }
        // Two more C stays below the growth.
        assert!(long_runs(&site(reference, "CCCCCCCCCTCCCCC", 1, Vec::new()), true, 3).is_empty());
        // A changed run structure is not a length.
        assert!(long_runs(&site(reference, "CCCCCCCCCCCCCCCC", 1, Vec::new()), true, 3).is_empty());
    }

    #[test]
    fn a_convention_leaves_a_long_run_length_unrecorded() {
        let convention = |window: [usize; 2], base: u8| Convention {
            name: "c_length".to_owned(),
            label: ConventionLabel::NotRecorded,
            kind: ConventionKind::RunLength {
                base,
                minimum_run: 10,
            },
            window,
        };
        // The site lies at 303-315; the call adds three bases.
        let at = site("CCCCCCCTCCCCC", "CCCCCCCCCCCCCCCC", 2, Vec::new());
        let long = run('C', 16, RunEvidence::ReferenceFrame);
        let short = run('C', 9, RunEvidence::ReferenceFrame);
        // Ten bases only with the three the site adds: the length may be one
        // the reviewers record.
        assert!(!not_recorded(
            &[convention([303, 315], b'C')],
            &at,
            &run('C', 12, RunEvidence::ReferenceFrame)
        ));
        assert!(not_recorded(
            &[convention([303, 315], b'C')],
            &at,
            &run('C', 13, RunEvidence::ReferenceFrame)
        ));
        assert!(not_recorded(&[convention([303, 315], b'C')], &at, &long));
        assert!(!not_recorded(&[convention([303, 315], b'C')], &at, &short));
        assert!(!not_recorded(&[convention([303, 315], b'A')], &at, &long));
        assert!(!not_recorded(
            &[convention([16180, 16193], b'C')],
            &at,
            &long
        ));
        assert!(!not_recorded(&[], &at, &long));
    }

    #[test]
    fn a_long_deletion_inside_a_run_is_lost_run() {
        // 2550268: CCCCCTCCCC decided as CTC.
        assert!(run_lost("CCCCCTCCCC", "CTC", 3, 4));
        // Two bases shorter stays below the growth.
        assert!(!run_lost("CCCCCTCCCC", "CCCCCCCC", 3, 4));
        // No run of four in the reference: a deletion elsewhere.
        assert!(!run_lost("CCGGAGCACCC", "CCCCC", 3, 4));
    }

    // SRS-API-014
    #[test]
    fn a_one_read_run_length_change_is_a_run_indel() {
        // LN_26_AE3814's 291DEL: A6 read as A5.
        assert!(run_length_change("AAAAAA", "AAAAA", 4));
        assert!(run_length_change("GCCCCA", "GCCCCCA", 4));
        // A run shorter than `phase_run` (249DEL in AA, 82.1A) is not.
        assert!(!run_length_change("AA", "A", 4));
        assert!(!run_length_change("AT", "AAT", 4));
        // A substitution, or a deletion that removes a whole run, is not.
        assert!(!run_length_change("CCTCC", "CCCCC", 4));
        assert!(!run_length_change("CCCCTC", "CCCCC", 4));
        assert!(!run_length_change("AAAAAA", "AAAAAA", 4));
    }

    #[test]
    fn an_in_phase_end_stops_at_an_unclean_call_or_inside_a_run() {
        let read = |text: &str| {
            text.chars()
                .map(|c| (c != '.').then_some(c))
                .collect::<Vec<_>>()
        };
        // A reverse read in reading order: the exit base, then into the run.
        assert_eq!(in_phase_end(read("ACCCCCCCC"), 4), vec!['A', 'C', 'C', 'C']);
        assert_eq!(in_phase_end(read("AC.CCCC"), 4), vec!['A', 'C']);
        assert_eq!(in_phase_end(read("AAAC"), 3), vec!['A', 'A']);
        assert_eq!(in_phase_end(read("ACGT"), 4), vec!['A', 'C', 'G', 'T']);
        assert!(in_phase_end(read(".ACGT"), 4).is_empty());
    }

    #[test]
    fn an_edge_is_the_first_or_last_positions_of_a_region() {
        let regions = [[16024, 16365], [73, 340]];
        assert!(near_edge(&regions, 16031, 16031, 20));
        assert!(near_edge(&regions, 16350, 16351, 20));
        assert!(!near_edge(&regions, 16044, 16044, 20));
        assert!(!near_edge(&regions, 200, 201, 20));
        assert!(!near_edge(&regions, 10, 10, 20));
        assert!(!near_edge(&regions, 16031, 16031, 0));
        // `106-111DEL` lies outside the deletion zone (20) of HVS-II's start,
        // but inside the wider edge zone beside a region end not analysed.
        assert!(!near_edge(&regions, 106, 111, 20));
        assert!(near_edge(&regions, 106, 111, 40));
    }

    #[test]
    fn dissent_is_across_or_within_one_strand_over_every_read() {
        let sides = |strands: &[Orientation], role: Role| {
            strands
                .iter()
                .enumerate()
                .map(|(index, &strand)| SideResult {
                    read: format!("{role:?}{index}"),
                    strand,
                    role,
                })
                .collect::<Vec<_>>()
        };
        let label = |supporting: &[Orientation], opposing: &[Orientation]| {
            dissent(
                sides(supporting, Role::Supporting),
                sides(opposing, Role::Opposing),
            )
            .map(|label| label.kind().name())
        };
        assert_eq!(label(&[F], &[R]), Some("strands_disagree"));
        assert_eq!(label(&[F, R], &[F]), Some("strands_disagree"));
        // Three reads of one strand against one: no strand contradicts, but
        // a clean read of the supporting strand dissents.
        assert_eq!(label(&[F, F, F], &[F]), Some("same_strand_dissent"));
        assert_eq!(label(&[R, R], &[R, R, R]), Some("same_strand_dissent"));
        assert_eq!(label(&[F, F, F], &[F, R]), Some("strands_disagree"));
        assert_eq!(label(&[F, F, F], &[]), None);
        assert_eq!(label(&[], &[R, R, R]), None);
        // The label lists every read, the supporting first.
        let Some(LabelResult::SameStrandDissent { reads }) = dissent(
            sides(&[F, F, F], Role::Supporting),
            sides(&[F], Role::Opposing),
        ) else {
            panic!("three forward reads against one dissent within the strand");
        };
        assert_eq!(
            reads.iter().map(|side| side.role).collect::<Vec<_>>(),
            [
                Role::Supporting,
                Role::Supporting,
                Role::Supporting,
                Role::Opposing
            ]
        );
        assert_eq!(one_strand(&[R, R]), Some(R));
        assert_eq!(one_strand(&[R, F]), None);
        assert_eq!(one_strand(&[]), None);
    }

    #[test]
    fn a_call_next_to_the_item_answers_its_difference() {
        let call = |position| NotationCall {
            position,
            ordinal: 0,
            reference: Some('C'),
            alternate: Some('T'),
            text: format!("{position}T"),
            convention: None,
        };
        assert!(called_near(&[call(9)], 10, 12));
        assert!(called_near(&[call(13)], 10, 12));
        assert!(called_near(&[call(11)], 10, 12));
        assert!(!called_near(&[call(8), call(14)], 10, 12));
        assert!(!called_near(&[], 10, 12));
    }

    #[test]
    fn a_call_shows_its_base_only_without_reasons_and_clear_of_the_edges() {
        const MIXED: EvidenceReason = EvidenceReason::new("mixed_supporting_dna");
        const DEPHASED: EvidenceReason = EvidenceReason::new("dephased");
        let trust = |placement, vetoes: Vec<EvidenceReason>, trusted, exempt| CallTrust {
            placement,
            unresolved: false,
            vetoes,
            trusted,
            exempt,
            withheld_by: Vec::new(),
            scope: TrustScope::Full,
            second: None,
        };
        assert!(shows(&trust(Placement::Scored, Vec::new(), true, false)));
        assert!(shows(&trust(Placement::Scored, Vec::new(), false, true)));
        assert!(!shows(&trust(Placement::Scored, Vec::new(), false, false)));
        assert!(!shows(&trust(Placement::Scored, vec![MIXED], true, false)));
        assert!(shows(&trust(
            Placement::Anchoring(DEPHASED),
            Vec::new(),
            false,
            false
        )));
        assert!(!shows(&trust(
            Placement::Anchoring(DEPHASED),
            vec![MIXED],
            false,
            false
        )));
        assert!(!shows(&trust(
            Placement::Unresolved(DEPHASED),
            Vec::new(),
            true,
            false
        )));
    }

    fn column(base: char, call: Option<usize>, reference: Option<usize>) -> AlignmentColumn {
        AlignmentColumn::new(base, 'C', call, reference)
    }

    fn clean_call(placement: Placement) -> CallTrust {
        CallTrust {
            placement,
            unresolved: false,
            vetoes: Vec::new(),
            trusted: true,
            exempt: false,
            withheld_by: Vec::new(),
            scope: TrustScope::Full,
            second: None,
        }
    }

    #[test]
    fn a_read_shows_an_item_only_between_two_clean_flanks() {
        const DEPHASED: EvidenceReason = EvidenceReason::new("dephased");
        // Reference positions 10-14 (0-based 9-13) read A C C C G; the item
        // is the C run 11-13. The read inserts a masked C after 13, and one
        // after 10, before the item.
        let columns = |last: char| {
            vec![
                column('A', Some(0), Some(9)),
                column('T', Some(1), None),
                column('C', Some(2), Some(10)),
                column('C', Some(3), Some(11)),
                column('C', Some(4), Some(12)),
                column('C', Some(5), None),
                column(last, (last != '-').then_some(6), Some(13)),
            ]
        };
        let mut trust = vec![clean_call(Placement::Scored); 7];
        trust[5] = clean_call(Placement::Anchoring(DEPHASED));
        // The masked inserted call extends the run; the earlier insertion
        // belongs to the position before the item.
        assert_eq!(
            bounded(&columns('G'), &trust, 11, 13).as_deref(),
            Some("CCCC")
        );
        // A flank the read deletes, or does not resolve, bounds nothing.
        assert_eq!(bounded(&columns('-'), &trust, 11, 13), None);
        assert_eq!(bounded(&columns('N'), &trust, 11, 13), None);
        // A flank or a call inside with a modality reason shows nothing.
        let mut vetoed = trust.clone();
        vetoed[6].vetoes = vec![DEPHASED];
        assert_eq!(bounded(&columns('G'), &vetoed, 11, 13), None);
        let mut inside = trust.clone();
        inside[3].trusted = false;
        assert_eq!(bounded(&columns('G'), &inside, 11, 13), None);
        // A read that does not reach a flank, or the reference's first
        // position, shows nothing.
        assert_eq!(bounded(&columns('G')[2..], &trust, 11, 13), None);
        assert_eq!(bounded(&columns('G'), &trust, 1, 3), None);
    }

    fn read(read: &str) -> String {
        read.to_owned()
    }

    // SRS-CONS-012
    #[test]
    fn a_difference_beside_a_region_end_not_analysed_is_flagged() {
        let regions = [[73, 340], [438, 576]];
        // 73-78 is not analysed: the region's start zone runs to 112.
        let cut = [[73, 78]];
        assert!(beside_not_analysed(&regions, &cut, 111, 111, 40));
        assert!(!beside_not_analysed(&regions, &cut, 113, 113, 40));
        // Another region's ends are not beside the cut.
        assert!(!beside_not_analysed(&regions, &cut, 440, 440, 40));
        // A region end cut at 570-576 flags the end zone.
        assert!(beside_not_analysed(&regions, &[[570, 576]], 540, 540, 40));
        assert!(!beside_not_analysed(&regions, &[[570, 576]], 536, 536, 40));
        // Without a cut, nothing is flagged.
        assert!(!beside_not_analysed(&regions, &[], 80, 80, 40));
    }

    // SRS-API-014
    #[test]
    fn labels_and_items_sort_by_priority() {
        let mut labels = vec![
            LabelResult::ReadEdge { read: read("b") },
            LabelResult::NotCovered,
            LabelResult::ReadEdge { read: read("a") },
            LabelResult::OneStrand {
                strand: F,
                reads: vec![read("a")],
            },
            LabelResult::Undecided,
            LabelResult::NotCalled,
            LabelResult::EstimatedLength { runs: Vec::new() },
            LabelResult::ReadDisagrees {
                read: read("a"),
                strand: F,
                sequence: "CC".to_owned(),
            },
            LabelResult::StrandsDisagree { reads: Vec::new() },
        ];
        labels.sort_by(|left, right| key(left).cmp(&key(right)));
        assert_eq!(
            labels
                .iter()
                .map(|label| label.kind().name())
                .collect::<Vec<_>>(),
            [
                "read_disagrees",
                "strands_disagree",
                "estimated_length",
                "not_called",
                "undecided",
                "one_strand",
                "not_covered",
                "read_edge",
                "read_edge",
            ]
        );
        assert_eq!(labels[8], LabelResult::ReadEdge { read: read("b") });
        assert_eq!(priority(&labels), 0);
        assert_eq!(priority(&labels[5..]), 3);
        assert_eq!(priority(&[]), UNLABELLED);
    }

    fn item(kind: ItemKind, labels: Vec<LabelResult>) -> ItemResult {
        ItemResult {
            first: ColumnResult {
                position: 1,
                ordinal: 0,
            },
            last: ColumnResult {
                position: 1,
                ordinal: 0,
            },
            kind,
            decision: Decision::Called,
            sequence: None,
            calls: Vec::new(),
            labels,
            evidence: Vec::new(),
        }
    }

    #[test]
    fn a_tie_is_undecided_even_when_settled_on_the_reference() {
        let site = |state, call: Option<&str>| ConsensusSite {
            first_0based: 10,
            last_0based: 10,
            reference: "C".to_owned(),
            call: call.map(str::to_owned),
            state,
            supporting: Vec::new(),
            opposing: Vec::new(),
            uninformative: Vec::new(),
            runs: Vec::new(),
        };
        // Two reads of one strand tie one against one: the reference is kept.
        assert!(undecided(&site(SiteState::Contested, Some("C"))));
        assert!(undecided(&site(SiteState::Contested, None)));
        assert!(undecided(&site(SiteState::Unresolved, None)));
        assert!(undecided(&site(SiteState::Called, Some("CN"))));
        assert!(!undecided(&site(SiteState::Called, Some("T"))));
    }

    // SRS-API-014
    #[test]
    fn only_labels_that_leave_an_item_unassured_make_a_check() {
        assert!(LabelResult::Undecided.unassured());
        assert!(LabelResult::NotCovered.unassured());
        assert!(
            !LabelResult::OneRead {
                read: "r".to_owned()
            }
            .unassured()
        );
        assert!(
            !LabelResult::ReadEdge {
                read: "r".to_owned()
            }
            .unassured()
        );
        // SRS-CONS-015: one strand cannot assure an insertion.
        let insertion = LabelResult::OneStrandInsertion {
            strand: F,
            reads: vec!["r".to_owned()],
            other: "not_covered",
            from_start: 25,
        };
        assert!(insertion.unassured());
        assert_eq!(
            (insertion.rank(), insertion.kind().name()),
            (2, "one_strand_insertion")
        );
    }

    // SRS-API-014
    #[test]
    fn summary_counts_items_per_label_and_clears_only_called_samples_without_checks() {
        let variant = item(ItemKind::Variant, Vec::new());
        let one_strand = item(
            ItemKind::Variant,
            vec![
                LabelResult::OneStrand {
                    strand: R,
                    reads: vec![read("a")],
                },
                LabelResult::ReadEdge { read: read("a") },
                LabelResult::ReadEdge { read: read("b") },
            ],
        );
        let summary = summarize(
            &[variant, one_strand],
            ReviewStatus::Called,
            &[[1, 10]],
            &[],
            &[],
            0,
        );
        assert_eq!((summary.items, summary.checks, summary.variants), (2, 0, 2));
        assert_eq!(summary.labels.get(LabelKind::OneStrand), 1);
        assert_eq!(summary.labels.get(LabelKind::ReadEdge), 1);
        assert_eq!(
            summary.positions,
            PositionCounts {
                reportable: 10,
                check: 0,
                reference_only: 0,
                not_analysed: 0
            }
        );
        assert!(summary.clear);

        let check = item(ItemKind::Check, vec![LabelResult::NotCovered]);
        let summary = summarize(&[check], ReviewStatus::Called, &[[1, 10]], &[], &[], 0);
        assert_eq!(
            (summary.checks, summary.labels.get(LabelKind::NotCovered)),
            (1, 1)
        );
        // The check item holds position 1 of the ten reportable positions.
        assert_eq!(
            summary.positions,
            PositionCounts {
                reportable: 10,
                check: 1,
                reference_only: 0,
                not_analysed: 0
            }
        );
        assert!(!summary.clear);
        assert!(!summarize(&[], ReviewStatus::CallFailed, &[[1, 10]], &[], &[], 0).clear);
        assert!(!summarize(&[], ReviewStatus::Rejected, &[[1, 10]], &[], &[], 0).clear);
        assert!(summarize(&[], ReviewStatus::Called, &[[1, 10]], &[], &[], 0).clear);
        // Reference-only positions are counted, and never make a check.
        let summary = summarize(&[], ReviewStatus::Called, &[[1, 10]], &[], &[3, 4], 0);
        assert_eq!((summary.positions.reference_only, summary.clear), (2, true));
        // Positions not analysed are counted, and never make a check
        // (SRS-CONS-012).
        let summary = summarize(&[], ReviewStatus::Called, &[[1, 10]], &[], &[], 6);
        assert_eq!((summary.positions.not_analysed, summary.clear), (6, true));
    }

    // SRS-API-014
    #[test]
    fn zones_follow_the_profile_windows_and_what_the_reads_show() {
        let span = |first, last, state, after_repeat| MaskedSpanResult {
            first,
            last,
            state,
            after_repeat,
        };
        let reads = [ReviewReadResult {
            name: "a".to_owned(),
            sha256: String::new(),
            strand: Some(F),
            admitted: true,
            rejection: None,
            masked: Some(vec![
                span(6, 6, PhaseState::Dephased, true),
                span(7, 7, PhaseState::Mixed, false),
            ]),
            unreliable: Some(Vec::new()),
            mixtures: Some(Vec::new()),
        }];
        let check = item(ItemKind::Check, vec![LabelResult::Undecided]);
        let variant = item(ItemKind::Variant, Vec::new());
        let zones = zones(
            &[[1, 10], [11, 13]],
            &[(4, 5, "W")],
            &[[1, 8], [11, 13]],
            &reads,
            &[check, variant],
        );
        let zone = |kind, first, last, name: Option<&str>, checks| ZoneResult {
            kind,
            first,
            last,
            name: name.map(str::to_owned),
            checks,
        };
        assert_eq!(
            zones,
            [
                // Only the check item at position 1 counts.
                zone(ZoneKind::Clean, 1, 3, None, 1),
                zone(ZoneKind::Repeat, 4, 5, Some("W"), 0),
                zone(ZoneKind::Dephased, 6, 6, None, 0),
                zone(ZoneKind::Masked, 7, 7, None, 0),
                zone(ZoneKind::Clean, 8, 8, None, 0),
                zone(ZoneKind::NotCallable, 9, 10, None, 0),
                // A zone never crosses into the next region.
                zone(ZoneKind::Clean, 11, 13, None, 0),
            ]
        );
        let summary = summarize(&[], ReviewStatus::Called, &[[1, 13]], &zones, &[], 0);
        assert_eq!(summary.zones, ZoneCounts { total: 7, check: 1 });
    }

    #[test]
    fn ranges_merge_clip_and_contain() {
        assert_eq!(segments([(0, 10), (5, 5), (12, 20)]), [[1, 10], [13, 20]]);
        assert_eq!(
            merge(vec![[13, 20], [1, 10], [11, 11], [30, 31]]),
            [[1, 11], [13, 20], [30, 31]]
        );
        assert_eq!(
            clip(&[[16024, 16365], [73, 340], [1, 50_000]], 400),
            [[1, 400]]
        );
        assert_eq!(
            clip(&[[16024, 16365], [73, 340], [438, 576]], 16569),
            [[73, 340], [438, 576], [16024, 16365]]
        );
        assert_eq!(clip(&[[500, 600]], 400), Vec::<[usize; 2]>::new());
        let regions = [[73, 340], [438, 576]];
        assert!(inside(&regions, 73) && inside(&regions, 576));
        assert!(!inside(&regions, 341) && !inside(&regions, 437));
        assert_eq!(within(&regions, 300, 500), Some([300, 500]));
        assert_eq!(within(&regions, 60, 80), Some([73, 80]));
        assert_eq!(within(&regions, 341, 437), None);
        assert_eq!(round(0.153_49), 0.153);
        assert_eq!(round(0.153_5), 0.154);
    }
}
