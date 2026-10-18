//! The review list of one sample, `dna.review/v1` (PROP-0004).
//!
//! Items and labels are computed from the sample analysis and its call alone,
//! never from a truth: the consensus sites, the reads' observations and the
//! core's trust in each call, and the Sanger scan position and relative
//! quality of each call, joined by read identity and call index. The rules are
//! version 1, the version of the `review` plugin.

use std::collections::{BTreeMap, BTreeSet};

use crate::error::{ReportError, Result};
use crate::model::read_observation::RejectedRead;
use crate::model::result::{IntervalResult, ReferenceResult};
use crate::model::review_result::{
    ColumnResult, Decision, EvidencePositionResult, EvidenceReadResult, ItemKind, ItemResult,
    LabelCounts, LabelResult, ReviewCallResult, ReviewProvenanceResult, ReviewReadResult,
    ReviewResult, ReviewStatus, Role, RunResult, ScopeResult, SecondResult, ShareResult,
    SideResult, SummaryResult,
};
use crate::model::sample_result::ConventionResult;
use crate::report::json::{project_plugins, project_profile};
use crate::report::notation::project_empop;
use crate::report::sample::{project_rejection, reviewer_names, reviewer_read_names};
use crate::sample_analysis::{CallTrust, ConsensusCall, NotationCall, SampleAnalysis, TrustScope};
use dna_core::model::alignment::{Orientation, ReferenceSegment};
use dna_core::model::consensus::{ConsensusSite, RunEvidence, SiteState};
use dna_core::model::sample_evidence::{LocusObservations, LocusState, SampleLocusObservation};
use dna_kernel::model::nucleotide::Nucleotide;
use dna_kernel::model::reference::Reference;
use dna_kernel::plugin::{Contract, PluginDescriptor, PluginFamily};
use dna_kernel::profile::Profile;
use dna_kernel::read_evidence::{EvidenceProfile, Placement};
use dna_sanger::model::attachment::SangerAttachment;

/// The review list: items and labels over the consensus, the call, and the
/// reads (PROP-0004). Its version is the labelling rules' version.
pub(crate) const PLUGIN: PluginDescriptor = PluginDescriptor {
    id: "review",
    family: PluginFamily::PostCalling,
    version: 1,
    provides: &[Contract::ReviewList],
    requires: &[
        Contract::ReadEvidence,
        Contract::Consensus,
        Contract::CalledVariants,
    ],
    config_sections: &[],
};

/// Smallest share of a read's second channel for `possible_heteroplasmy`
/// (`>=`).
const MINIMUM_SHARE: f64 = 0.15;

/// Priority of an unlabelled variant, after every label's rank.
const UNLABELLED: u8 = 5;

/// Inputs consumed to build one review document.
pub(crate) struct CompletedReview {
    pub(crate) sample_id: String,
    pub(crate) reference: Reference,
    pub(crate) profile: Profile,
    pub(crate) configuration_sha256: String,
    pub(crate) outcome: Reviewed,
    /// Plugins of the workflow, in execution order.
    pub(crate) plugins: &'static [&'static PluginDescriptor],
}

/// What the sample analysis and its call gave.
pub(crate) enum Reviewed {
    /// Every read was rejected; each with its cause, by content SHA-256.
    Rejected(BTreeMap<String, RejectedRead>),
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
    };
    let (analysis, call) = match outcome {
        Reviewed::Rejected(rejected) => {
            let names = reviewer_names(rejected.values().map(|read| read.input_name.as_str()))?;
            let reads = rejected
                .values()
                .zip(names)
                .map(|(read, name)| ReviewReadResult {
                    name,
                    sha256: read.input_sha256.clone(),
                    strand: None,
                    admitted: false,
                    rejection: Some(project_rejection(&read.cause)),
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
                empop: None,
                summary: summarize(&[], ReviewStatus::Rejected),
                items: Vec::new(),
            });
        }
        Reviewed::Analyzed { analysis, call } => (analysis, call),
    };

    let (names, rejected_names) = reviewer_read_names(&analysis.evidence)?;
    let mut reads = analysis
        .evidence
        .reads
        .iter()
        .zip(&names)
        .map(|(read, name)| ReviewReadResult {
            name: name.clone(),
            sha256: read.input_sha256.clone(),
            strand: Some(read.alignment.orientation),
            admitted: true,
            rejection: None,
        })
        .collect::<Vec<_>>();
    for (read, name) in analysis.evidence.rejected_reads.iter().zip(rejected_names) {
        let rejected =
            analysis
                .rejected
                .get(&read.input_sha256)
                .ok_or(ReportError::Inconsistent(
                    "rejected review read lacks its rejection cause",
                ))?;
        reads.push(ReviewReadResult {
            name,
            sha256: read.input_sha256.clone(),
            strand: None,
            admitted: false,
            rejection: Some(project_rejection(&rejected.cause)),
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
    let sample = Sample::new(&analysis, &names, &reference, &regions, &callable, notation)?;
    let items = sample.items()?;
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
        empop: call.and_then(|call| call.empop).map(project_empop),
        summary: summarize(&items, status),
        items,
    })
}

/// One sample's evidence, indexed for the review rules.
struct Sample<'a> {
    sequence: &'a [u8],
    regions: &'a [[usize; 2]],
    callable: &'a [[usize; 2]],
    /// Reference segments of the admitted reads, merged, 1-based inclusive.
    covered: Vec<[usize; 2]>,
    calls: &'a [NotationCall],
    /// Whether DNA's calls are in the profile notation; without them, a
    /// difference is the consensus's decision.
    notation: bool,
    sites: &'a [ConsensusSite],
    /// Per read, in evidence order: name, content SHA-256, strand, the core's
    /// trust in each call, and the Sanger evidence.
    reads: Vec<Read<'a>>,
    /// Read index by input name, as consensus sites name reads.
    by_input: BTreeMap<&'a str, usize>,
    /// Every covered position's observations.
    loci: BTreeMap<usize, LocusObservations>,
}

struct Read<'a> {
    name: &'a str,
    sha256: &'a str,
    strand: Orientation,
    trust: Vec<CallTrust>,
    sanger: &'a SangerAttachment,
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
    ) -> Result<Self> {
        let reads = analysis
            .evidence
            .reads
            .iter()
            .zip(names)
            .map(|(read, name)| {
                let missing = || ReportError::Inconsistent("review read lacks its evidence");
                Ok(Read {
                    name,
                    sha256: &read.input_sha256,
                    strand: read.alignment.orientation,
                    trust: analysis.trust(&read.input_sha256).ok_or_else(missing)?,
                    sanger: analysis.reads.get(&read.input_sha256).ok_or_else(missing)?,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let by_input = analysis
            .evidence
            .reads
            .iter()
            .enumerate()
            .map(|(index, read)| (read.input_name.as_str(), index))
            .collect();
        let covered = merge(segments(
            analysis
                .evidence
                .reads
                .iter()
                .flat_map(|read| &read.alignment.reference_segments)
                .map(span),
        ));
        let loci = analysis
            .observations(1, reference.len())?
            .into_iter()
            .map(|locus| (locus.position_1based, locus))
            .collect();
        Ok(Self {
            sequence: reference.sequence.as_bytes(),
            regions,
            callable,
            covered,
            calls: notation.unwrap_or_default(),
            notation: notation.is_some(),
            sites: &analysis.consensus.sites,
            reads,
            by_input,
            loci,
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
                drafts.push(Draft {
                    first,
                    last,
                    site: Some(site),
                    sequence: site.call.clone(),
                    decided: site.call.as_ref().is_some_and(|call| !call.contains('N')),
                    labels: self.site_labels(site),
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

        let mut items = drafts
            .into_iter()
            .map(|draft| self.finish(draft))
            .filter_map(Result::transpose)
            .collect::<Result<Vec<_>>>()?;
        items.sort_by_key(|item| (priority(&item.labels), item.first, item.last));
        Ok(items)
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
                .filter_map(|input| self.by_input.get(input.as_str()))
                .map(|&index| SideResult {
                    read: self.reads[index].name.to_owned(),
                    strand: self.reads[index].strand,
                    role,
                })
                .collect::<Vec<_>>()
        };
        let supporting = sides(&site.supporting, Role::Supporting);
        let opposing = sides(&site.opposing, Role::Opposing);
        if disagree(
            &supporting
                .iter()
                .map(|side| side.strand)
                .collect::<Vec<_>>(),
            &opposing.iter().map(|side| side.strand).collect::<Vec<_>>(),
        ) {
            labels.push(LabelResult::StrandsDisagree {
                reads: supporting.into_iter().chain(opposing).collect(),
            });
        }
        if undecided(site) {
            labels.push(LabelResult::Undecided);
        }
        let runs = site
            .runs
            .iter()
            .filter(|run| {
                matches!(
                    run.evidence,
                    RunEvidence::PhaseLoss | RunEvidence::ReferenceFrame
                )
            })
            .map(|run| RunResult {
                base: run.base,
                length: run.length,
                length_evidence: run.evidence.label(),
            })
            .collect::<Vec<_>>();
        if !runs.is_empty() {
            labels.push(LabelResult::EstimatedLength { runs });
        }
        labels
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
        // otherwise the consensus's own decision.
        let differs = if self.notation {
            !calls.is_empty()
        } else {
            site.is_some_and(|site| {
                decided && site.call.as_deref() != Some(site.reference.as_str())
            })
        };
        if let (true, Some(site)) = (differs, site) {
            let supporting = site
                .supporting
                .iter()
                .filter_map(|input| self.by_input.get(input.as_str()))
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
            }
        }
        labels.extend((first..=last).filter_map(|position| self.heteroplasmy(position)));
        let (read_labels, evidence) = self.reads_over(first, last)?;
        labels.extend(read_labels);
        // A label that leaves the item unassured makes it a check; the others
        // are context for the reviewer. An item with neither a call nor such a
        // label is DNA's own decision and is not listed.
        let assured = !labels.iter().any(unassured);
        if !differs && assured {
            return Ok(None);
        }
        let kind = if assured {
            ItemKind::Variant
        } else {
            ItemKind::Check
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
                            .entry((rank(&label), name(&label), index))
                            .or_insert(label);
                    }
                }
                let scope = observation
                    .call_index_0based
                    .and_then(|call| read.trust.get(call))
                    .map(|trust| match trust.scope {
                        TrustScope::Full => ScopeResult::Full,
                        TrustScope::Base => ScopeResult::Base,
                    });
                positions.entry(index).or_default().push(evidence(
                    position,
                    observation,
                    read.sanger,
                    scope,
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
                let (second, share) = second(observation.base?, observation.profile?)?;
                let neighbour = |other: Option<usize>| {
                    other
                        .and_then(|other| self.observation(observation.read_index, other))
                        .and_then(|other| other.base)
                };
                Some(Clean {
                    read: observation.read_index,
                    strand: read.strand,
                    second,
                    share,
                    neighbours: [
                        neighbour(position.checked_sub(1)),
                        neighbour(position.checked_add(1)),
                    ],
                })
            })
            .collect::<Vec<_>>();
        let (base, reads) = heteroplasmy(&clean)?;
        Some(LabelResult::PossibleHeteroplasmy {
            position,
            base,
            reads: reads
                .into_iter()
                .map(|clean| ShareResult {
                    read: self.reads[clean.read].name.to_owned(),
                    strand: clean.strand,
                    share: round(clean.share),
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
fn evidence(
    position: usize,
    observation: &SampleLocusObservation,
    sanger: &SangerAttachment,
    scope: Option<ScopeResult>,
) -> Result<EvidencePositionResult> {
    let base = match observation.state {
        LocusState::Deletion => '-',
        LocusState::Unresolved => 'N',
        _ => observation.base.unwrap_or('N'),
    };
    let (scan, quality) = match observation.call_index_0based {
        Some(index) => {
            let call = sanger
                .calls
                .calls
                .get(index)
                .filter(|call| call.index_0based == index);
            let quality = sanger
                .quality
                .per_call
                .get(index)
                .filter(|quality| quality.index_0based == index);
            let (Some(call), Some(quality)) = (call, quality) else {
                return Err(ReportError::Inconsistent(
                    "review evidence call lacks its Sanger record",
                )
                .into());
            };
            (
                Some(call.locus_position_0based),
                Some(quality.relative_quality_score),
            )
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
        second: observation
            .profile
            .zip(observation.base)
            .filter(|_| observation.state != LocusState::Unresolved)
            .and_then(|(profile, base)| second(base, profile))
            .map(|(base, share)| SecondResult {
                base,
                share: round(share),
            }),
    })
}

/// One read that reads a position cleanly: trusted or exempt, unmasked, with
/// a profile.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Clean {
    /// Read index.
    read: usize,
    strand: Orientation,
    /// The strongest profile channel other than the call, and its weight.
    second: char,
    share: f64,
    /// The bases the same read shows at the previous and the next reference
    /// position, when it shows one.
    neighbours: [Option<char>; 2],
}

/// `possible_heteroplasmy` over the clean reads of one position: a read
/// qualifies when its second share is at least [`MINIMUM_SHARE`] and its
/// second base is not a base it shows next to the position. When both strands
/// read the position cleanly, a qualifying read on each strand must show the
/// same second base; when one strand does, one qualifying read suffices. The
/// base most qualifying reads show wins, then the first in A/C/G/T order. It
/// returns that base and its qualifying reads, in input order.
fn heteroplasmy(clean: &[Clean]) -> Option<(char, Vec<Clean>)> {
    let qualifying = clean
        .iter()
        .filter(|read| read.share >= MINIMUM_SHARE && !read.neighbours.contains(&Some(read.second)))
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
        let fires = if both {
            shows(base, Orientation::Forward) && shows(base, Orientation::Reverse)
        } else {
            count > 0
        };
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

/// The strongest profile channel other than `base`, the first in A/C/G/T
/// order on a tie; `None` for a base that is not A, C, G, or T.
fn second(base: char, profile: EvidenceProfile) -> Option<(char, f64)> {
    let own = Nucleotide::ALL
        .iter()
        .position(|nucleotide| nucleotide.as_char() == base)?;
    let mut best: Option<(char, f64)> = None;
    for (index, nucleotide) in Nucleotide::ALL.iter().enumerate() {
        let weight = profile.weights[index];
        if index != own && best.is_none_or(|(_, most)| weight > most) {
            best = Some((nucleotide.as_char(), weight));
        }
    }
    best
}

/// Whether a supporting read and an opposing read lie on different strands.
fn disagree(supporting: &[Orientation], opposing: &[Orientation]) -> bool {
    supporting
        .iter()
        .any(|left| opposing.iter().any(|right| left != right))
}

/// The one strand every supporting read lies on; `None` without a read or
/// with reads on both strands.
fn one_strand(supporting: &[Orientation]) -> Option<Orientation> {
    let (&first, rest) = supporting.split_first()?;
    rest.iter().all(|&strand| strand == first).then_some(first)
}

/// Whether the consensus left a site undecided: no clean observation, a
/// tie (even one the consensus settled on the reference, whichever strands
/// the tied reads are on), or an `N` in its call.
fn undecided(site: &ConsensusSite) -> bool {
    matches!(site.state, SiteState::Unresolved | SiteState::Contested)
        || site.call.as_ref().is_some_and(|call| call.contains('N'))
}

/// Whether a label leaves its item unassured, so that the item is a check:
/// strands that disagree, a possible heteroplasmy, an estimated length, an
/// undecided position, or no read. The read labels (one read, one strand,
/// noisy signal, read edge, phase shift) are context for the reviewer.
const fn unassured(label: &LabelResult) -> bool {
    matches!(
        label,
        LabelResult::StrandsDisagree { .. }
            | LabelResult::PossibleHeteroplasmy { .. }
            | LabelResult::EstimatedLength { .. }
            | LabelResult::Undecided
            | LabelResult::NotCovered
    )
}

/// A label's rank: the lower, the earlier a reviewer looks.
const fn rank(label: &LabelResult) -> u8 {
    match label {
        LabelResult::StrandsDisagree { .. } => 0,
        LabelResult::PossibleHeteroplasmy { .. } => 1,
        LabelResult::EstimatedLength { .. } | LabelResult::Undecided => 2,
        LabelResult::OneRead { .. } | LabelResult::OneStrand { .. } => 3,
        LabelResult::NoisySignal { .. }
        | LabelResult::NotCovered
        | LabelResult::PhaseShift { .. }
        | LabelResult::ReadEdge { .. } => 4,
    }
}

/// A label's published name.
const fn name(label: &LabelResult) -> &'static str {
    match label {
        LabelResult::StrandsDisagree { .. } => "strands_disagree",
        LabelResult::PossibleHeteroplasmy { .. } => "possible_heteroplasmy",
        LabelResult::EstimatedLength { .. } => "estimated_length",
        LabelResult::Undecided => "undecided",
        LabelResult::OneRead { .. } => "one_read",
        LabelResult::OneStrand { .. } => "one_strand",
        LabelResult::NoisySignal { .. } => "noisy_signal",
        LabelResult::NotCovered => "not_covered",
        LabelResult::PhaseShift { .. } => "phase_shift",
        LabelResult::ReadEdge { .. } => "read_edge",
    }
}

/// Label order within an item: rank, name, then the position or read.
fn key(label: &LabelResult) -> (u8, &'static str, usize, &str) {
    let (position, read) = match label {
        LabelResult::PossibleHeteroplasmy { position, .. } => (*position, ""),
        LabelResult::OneRead { read }
        | LabelResult::NoisySignal { read, .. }
        | LabelResult::PhaseShift { read, .. }
        | LabelResult::ReadEdge { read } => (0, read.as_str()),
        _ => (0, ""),
    };
    (rank(label), name(label), position, read)
}

/// An item's priority: its best label rank, or [`UNLABELLED`].
fn priority(labels: &[LabelResult]) -> u8 {
    labels.iter().map(rank).min().unwrap_or(UNLABELLED)
}

/// Counts over the items; `clear` when the call succeeded and no item is a
/// check.
fn summarize(items: &[ItemResult], status: ReviewStatus) -> SummaryResult {
    let mut labels = LabelCounts::default();
    for item in items {
        let names = item.labels.iter().map(name).collect::<BTreeSet<_>>();
        for label in names {
            let count = match label {
                "strands_disagree" => &mut labels.strands_disagree,
                "possible_heteroplasmy" => &mut labels.possible_heteroplasmy,
                "estimated_length" => &mut labels.estimated_length,
                "undecided" => &mut labels.undecided,
                "one_read" => &mut labels.one_read,
                "one_strand" => &mut labels.one_strand,
                "noisy_signal" => &mut labels.noisy_signal,
                "not_covered" => &mut labels.not_covered,
                "phase_shift" => &mut labels.phase_shift,
                _ => &mut labels.read_edge,
            };
            *count += 1;
        }
    }
    let checks = items
        .iter()
        .filter(|item| item.kind == ItemKind::Check)
        .count();
    SummaryResult {
        items: items.len(),
        checks,
        variants: items.len() - checks,
        labels,
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
fn round(share: f64) -> f64 {
    (share * 1000.0).round() / 1000.0
}

#[cfg(test)]
mod tests {
    use super::*;

    const F: Orientation = Orientation::Forward;
    const R: Orientation = Orientation::Reverse;

    fn clean(read: usize, strand: Orientation, second: char, share: f64) -> Clean {
        Clean {
            read,
            strand,
            second,
            share,
            neighbours: [None, None],
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
    fn heteroplasmy_on_one_strand_needs_one_read() {
        let reads = [clean(0, R, 'A', 0.3), clean(1, R, 'C', 0.01)];
        let (base, qualifying) = heteroplasmy(&reads).unwrap_or(('-', Vec::new()));
        assert_eq!(base, 'A');
        assert_eq!(qualifying, [reads[0]]);
        // The base most reads show wins, then A/C/G/T order.
        let tied = [
            clean(0, F, 'T', 0.2),
            clean(1, F, 'C', 0.2),
            clean(2, F, 'T', 0.2),
        ];
        assert_eq!(heteroplasmy(&tied).map(|(base, _)| base), Some('T'));
        let even = [clean(0, F, 'T', 0.2), clean(1, F, 'C', 0.2)];
        assert_eq!(heteroplasmy(&even).map(|(base, _)| base), Some('C'));
    }

    #[test]
    fn heteroplasmy_excludes_a_neighbour_shadow_and_a_share_below_the_threshold() {
        let mut shadow = clean(0, F, 'G', 0.4);
        shadow.neighbours = [Some('A'), Some('G')];
        assert_eq!(heteroplasmy(&[shadow]), None);
        shadow.neighbours = [Some('G'), None];
        assert_eq!(heteroplasmy(&[shadow]), None);
        shadow.neighbours = [Some('A'), Some('C')];
        assert!(heteroplasmy(&[shadow]).is_some());
        assert_eq!(heteroplasmy(&[clean(0, F, 'G', 0.149)]), None);
        assert!(heteroplasmy(&[clean(0, F, 'G', 0.15)]).is_some());
        assert_eq!(heteroplasmy(&[]), None);
    }

    #[test]
    fn second_is_the_strongest_other_channel() {
        let profile = EvidenceProfile {
            weights: [0.6, 0.1, 0.25, 0.05],
        };
        assert_eq!(second('A', profile), Some(('G', 0.25)));
        assert_eq!(second('G', profile), Some(('A', 0.6)));
        assert_eq!(second('N', profile), None);
        let tie = EvidenceProfile {
            weights: [0.5, 0.25, 0.25, 0.0],
        };
        assert_eq!(second('A', tie), Some(('C', 0.25)));
    }

    #[test]
    fn strands_disagree_needs_opposite_strands_across_the_sides() {
        assert!(disagree(&[F], &[R]));
        assert!(disagree(&[F, R], &[F]));
        assert!(!disagree(&[F, F], &[F]));
        assert!(!disagree(&[F], &[]));
        assert!(!disagree(&[], &[R]));
        assert_eq!(one_strand(&[R, R]), Some(R));
        assert_eq!(one_strand(&[R, F]), None);
        assert_eq!(one_strand(&[]), None);
    }

    fn read(read: &str) -> String {
        read.to_owned()
    }

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
            LabelResult::EstimatedLength { runs: Vec::new() },
            LabelResult::StrandsDisagree { reads: Vec::new() },
        ];
        labels.sort_by(|left, right| key(left).cmp(&key(right)));
        assert_eq!(
            labels.iter().map(name).collect::<Vec<_>>(),
            [
                "strands_disagree",
                "estimated_length",
                "undecided",
                "one_strand",
                "not_covered",
                "read_edge",
                "read_edge",
            ]
        );
        assert_eq!(labels[6], LabelResult::ReadEdge { read: read("b") });
        assert_eq!(priority(&labels), 0);
        assert_eq!(priority(&labels[3..]), 3);
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

    #[test]
    fn only_labels_that_leave_an_item_unassured_make_a_check() {
        assert!(unassured(&LabelResult::Undecided));
        assert!(unassured(&LabelResult::NotCovered));
        assert!(!unassured(&LabelResult::OneRead {
            read: "r".to_owned()
        }));
        assert!(!unassured(&LabelResult::ReadEdge {
            read: "r".to_owned()
        }));
    }

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
        let summary = summarize(&[variant, one_strand], ReviewStatus::Called);
        assert_eq!((summary.items, summary.checks, summary.variants), (2, 0, 2));
        assert_eq!(summary.labels.one_strand, 1);
        assert_eq!(summary.labels.read_edge, 1);
        assert!(summary.clear);

        let check = item(ItemKind::Check, vec![LabelResult::NotCovered]);
        let summary = summarize(&[check], ReviewStatus::Called);
        assert_eq!((summary.checks, summary.labels.not_covered), (1, 1));
        assert!(!summary.clear);
        assert!(!summarize(&[], ReviewStatus::CallFailed).clear);
        assert!(!summarize(&[], ReviewStatus::Rejected).clear);
        assert!(summarize(&[], ReviewStatus::Called).clear);
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
