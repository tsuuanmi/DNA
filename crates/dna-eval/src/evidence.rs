//! The compact triage evidence of one sample: what the draft stage-1
//! detectors of [`crate::triage`] read, taken once by the runner from DNA's
//! sample analysis (stage E), so that every policy runs over it without DNA
//! (stage R).
//!
//! Per reportable region of the target profile it keeps, for each reference
//! position, the reads that observe it cleanly on each strand, the reads that
//! cover it only through a masked call, and the consensus marks over it. Per
//! input trace it keeps the name, the primer the name carries, whether the
//! read was admitted (with its strand and the clean positions it reads in
//! each region) or why it was rejected, and its callable span. Per region
//! position it also keeps the clean reads' minor channels ([`minor`]).

pub(crate) mod minor;

use std::collections::BTreeMap;

use dna::sample_analysis::{Consensus, RejectedRead, RejectionCause, SampleAnalysis};
use dna_align::model::alignment::{Orientation, ReferenceSegment};
use dna_consensus::model::consensus::{ConsensusSite, RunEvidence, SiteState};
use dna_sample::model::sample_evidence::SampleReadEvidence;
use dna_sanger::model::callability::ReadCallability;

/// The strand a read observes the reference on.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Strand {
    Forward,
    Reverse,
}

impl Strand {
    pub(crate) const fn label(self) -> &'static str {
        match self {
            Self::Forward => "forward",
            Self::Reverse => "reverse",
        }
    }

    /// The strand a read placed in `orientation` observes.
    pub(crate) const fn of(orientation: Orientation) -> Self {
        match orientation {
            Orientation::Forward => Self::Forward,
            Orientation::Reverse => Self::Reverse,
        }
    }
}

/// A cause of `region.unreliable`, read from one position's evidence.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Cause {
    /// The consensus has no definite base: an `N` or an undecided site.
    N,
    /// The consensus decided the site by a tie.
    Contested,
    /// A consensus site that a covering read observes only through a masked
    /// or untrusted call (its `uninformative` reads).
    Masked,
    /// Any covering read observes the position only through a masked call.
    MaskedRead,
    /// A run of the site has the length evidence `phase_loss`.
    PhaseLoss,
    /// A run of the site has the length evidence `reference_frame`.
    ReferenceFrame,
}

impl Cause {
    pub(crate) const ALL: [Self; 6] = [
        Self::N,
        Self::Contested,
        Self::Masked,
        Self::MaskedRead,
        Self::PhaseLoss,
        Self::ReferenceFrame,
    ];

    pub(crate) const fn label(self) -> &'static str {
        match self {
            Self::N => "n",
            Self::Contested => "contested",
            Self::Masked => "masked",
            Self::MaskedRead => "masked_read",
            Self::PhaseLoss => "phase_loss",
            Self::ReferenceFrame => "reference_frame",
        }
    }

    /// The consensus mark bit; `None` for a cause read from the reads.
    const fn bit(self) -> Option<u8> {
        match self {
            Self::N => Some(1),
            Self::Contested => Some(2),
            Self::Masked => Some(4),
            Self::PhaseLoss => Some(8),
            Self::ReferenceFrame => Some(16),
            Self::MaskedRead => None,
        }
    }
}

/// The evidence at one reference position.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) struct Position {
    /// Forward reads that observe the position cleanly.
    pub(crate) forward: u8,
    /// Reverse reads that observe the position cleanly.
    pub(crate) reverse: u8,
    /// Reads that cover the position only through a masked call.
    pub(crate) masked_reads: u8,
    /// Consensus marks, one bit per [`Cause`] the consensus shows.
    pub(crate) marks: u8,
}

impl Position {
    /// Whether the position shows `cause`.
    pub(crate) const fn has(self, cause: Cause) -> bool {
        match cause.bit() {
            Some(bit) => self.marks & bit != 0,
            None => self.masked_reads > 0,
        }
    }

    const fn mark(&mut self, cause: Cause) {
        if let Some(bit) = cause.bit() {
            self.marks |= bit;
        }
    }

    /// The strands that observe the position cleanly.
    pub(crate) const fn strands(self) -> usize {
        (self.forward > 0) as usize + (self.reverse > 0) as usize
    }
}

/// One reportable region and the evidence at each of its positions.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Region {
    /// First position, 1-based.
    pub(crate) start: usize,
    /// Last position, 1-based and inclusive.
    pub(crate) end: usize,
    /// One entry per position from `start` to `end`.
    pub(crate) positions: Vec<Position>,
}

impl Region {
    fn new(start: usize, end: usize) -> Self {
        Self {
            start,
            end,
            positions: vec![Position::default(); end.saturating_sub(start) + 1],
        }
    }

    /// `16024-16365`.
    pub(crate) fn label(&self) -> String {
        format!("{}-{}", self.start, self.end)
    }

    /// The positions of a 0-based half-open reference interval that lie in
    /// the region, as indices into `positions`.
    fn overlap(&self, start_0based: usize, end_0based_exclusive: usize) -> std::ops::Range<usize> {
        let first = (start_0based + 1).max(self.start);
        let last = end_0based_exclusive.min(self.end);
        if first > last {
            return 0..0;
        }
        first - self.start..last - self.start + 1
    }
}

/// How one input trace ended up.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum TraceState {
    /// Admitted: its strand and the clean positions it reads per region.
    Admitted { strand: Strand, clean: Vec<usize> },
    /// Fewer calls than one signal window.
    TooShort,
    /// Too few callable calls.
    Weak,
    /// The core could not place it, with the placement rejection's label.
    Unplaced(&'static str),
    /// Rejected for a cause this draft does not know.
    Rejected,
}

/// One input trace.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Trace {
    pub(crate) name: String,
    /// The primer the name carries (`HV1F`), upper case, if any.
    pub(crate) primer: Option<String>,
    pub(crate) state: TraceState,
    /// Callable span in calls, from the first to the last unmasked call, when
    /// the read reached callability.
    pub(crate) span: Option<usize>,
}

/// The triage evidence of one sample.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Evidence {
    pub(crate) regions: Vec<Region>,
    /// Traces sorted by name, then by state.
    pub(crate) traces: Vec<Trace>,
    /// Per region and position, the clean reads' minor channels; empty when
    /// not taken.
    pub(crate) minor: minor::Signals,
    /// Why the minor channels could not be taken, if they could not. A
    /// policy with a minor-signal rule then sends regions that would pass to
    /// review, because the gate could not run.
    pub(crate) minor_failure: Option<String>,
}

impl Evidence {
    /// The region holding the 1-based position `anchor`, and the anchor's
    /// offset from the region's start.
    pub(crate) fn locate(&self, anchor: i64) -> Option<(usize, usize)> {
        let anchor = usize::try_from(anchor).ok()?;
        self.regions.iter().enumerate().find_map(|(index, region)| {
            (region.start..=region.end)
                .contains(&anchor)
                .then(|| (index, anchor - region.start))
        })
    }

    /// The minor channels at position `offset` of region `region`; none when
    /// not taken.
    pub(crate) fn minors(&self, region: usize, offset: usize) -> &[minor::Minor] {
        self.minor
            .get(region)
            .and_then(|positions| positions.get(offset))
            .map_or(&[], Vec::as_slice)
    }
}

/// Takes the triage evidence of one analysed sample. When DNA cannot give
/// the reads' observations, the stage-1 evidence is kept without minor
/// channels and [`Evidence::minor_failure`] says why.
pub(crate) fn extract(analysis: &SampleAnalysis) -> Evidence {
    let spans = analysis
        .reads
        .iter()
        .map(|(sha256, sanger)| (sha256.as_str(), span(&sanger.callability)))
        .collect::<BTreeMap<_, _>>();
    let (minor, minor_failure) = match minor::extract(analysis) {
        Ok(minor) => (minor, None),
        Err(error) => (
            Vec::new(),
            Some(format!("cannot read the observations: {error}")),
        ),
    };
    Evidence {
        minor,
        minor_failure,
        ..build(
            analysis.regions(),
            &analysis.evidence.reads,
            &analysis.consensus,
            &spans,
            rejected(analysis.rejected.values()),
        )
    }
}

/// The evidence of a sample whose every read was rejected: its traces only,
/// with no region evidence.
pub(crate) fn rejected_only<'a>(reads: impl IntoIterator<Item = &'a RejectedRead>) -> Evidence {
    let mut traces = rejected(reads);
    traces.sort_by(|left, right| order(left).cmp(&order(right)));
    Evidence {
        regions: Vec::new(),
        traces,
        minor: Vec::new(),
        minor_failure: None,
    }
}

/// The traces of rejected reads, with their cause and callable span.
fn rejected<'a>(reads: impl IntoIterator<Item = &'a RejectedRead>) -> Vec<Trace> {
    reads
        .into_iter()
        .map(|read| {
            let state = match &read.cause {
                RejectionCause::TooShort(_) => TraceState::TooShort,
                RejectionCause::Callability { .. } => TraceState::Weak,
                RejectionCause::Placement { rejection, .. } => {
                    TraceState::Unplaced(rejection.label())
                }
                _ => TraceState::Rejected,
            };
            Trace {
                name: read.input_name.clone(),
                primer: primer(&read.input_name),
                state,
                span: read.cause.sanger().map(|sanger| span(&sanger.callability)),
            }
        })
        .collect()
}

/// The callable span of a read in calls: from its first to its last unmasked
/// call.
fn span(callability: &ReadCallability) -> usize {
    callability
        .callable_end_0based_exclusive
        .saturating_sub(callability.callable_start_0based)
}

/// Builds the evidence from the profile's regions, the admitted reads, the
/// consensus, the admitted reads' callable spans by content SHA-256, and the
/// rejected traces.
///
/// Each position has one source for its clean and masked reads: at a
/// consensus site, the site's supporting and opposing reads are its clean
/// observers and its uninformative reads its masked ones; elsewhere, the reads'
/// callable reference segments are clean and the rest of their reference
/// segments masked.
pub(crate) fn build(
    regions: &[[usize; 2]],
    reads: &[SampleReadEvidence],
    consensus: &Consensus,
    spans: &BTreeMap<&str, usize>,
    rejected: Vec<Trace>,
) -> Evidence {
    let mut regions = regions
        .iter()
        .map(|[start, end]| Region::new(*start, *end))
        .collect::<Vec<_>>();
    let mut traces = rejected;
    for read in reads {
        let strand = Strand::of(read.alignment.orientation);
        let mut clean = vec![0; regions.len()];
        for (index, region) in regions.iter_mut().enumerate() {
            // Every covered position counts as masked until a clean call
            // observes it.
            for range in overlaps(region, &read.alignment.reference_segments) {
                for position in &mut region.positions[range] {
                    position.masked_reads = position.masked_reads.saturating_add(1);
                }
            }
            for range in overlaps(region, &read.alignment.callable_reference_segments) {
                clean[index] += range.len();
                for position in &mut region.positions[range] {
                    position.masked_reads = position.masked_reads.saturating_sub(1);
                    match strand {
                        Strand::Forward => position.forward = position.forward.saturating_add(1),
                        Strand::Reverse => position.reverse = position.reverse.saturating_add(1),
                    }
                }
            }
        }
        traces.push(Trace {
            name: read.input_name.clone(),
            primer: primer(&read.input_name),
            state: TraceState::Admitted { strand, clean },
            span: spans.get(read.input_sha256.as_str()).copied(),
        });
    }
    // Consensus sites name their reads by content SHA-256, which is unique.
    let strands = reads
        .iter()
        .map(|read| {
            (
                read.input_sha256.as_str(),
                Strand::of(read.alignment.orientation),
            )
        })
        .collect::<BTreeMap<_, _>>();
    for site in &consensus.sites {
        let marks = site_causes(site);
        let observed = observers(&strands, &site.supporting, &site.opposing);
        let masked = u8::try_from(site.uninformative.len()).unwrap_or(u8::MAX);
        for region in &mut regions {
            let range = region.overlap(site.first_0based, site.last_0based + 1);
            for position in &mut region.positions[range] {
                position.forward = observed.forward;
                position.reverse = observed.reverse;
                position.masked_reads = masked;
                for cause in &marks {
                    position.mark(*cause);
                }
            }
        }
    }
    traces.sort_by(|left, right| order(left).cmp(&order(right)));
    Evidence {
        regions,
        traces,
        minor: Vec::new(),
        minor_failure: None,
    }
}

/// The unreliable causes a consensus site shows.
fn site_causes(site: &ConsensusSite) -> Vec<Cause> {
    let mut causes = Vec::new();
    if site.call.as_ref().is_none_or(|call| call.contains('N')) {
        causes.push(Cause::N);
    }
    if site.state == SiteState::Contested {
        causes.push(Cause::Contested);
    }
    if !site.uninformative.is_empty() {
        causes.push(Cause::Masked);
    }
    for run in &site.runs {
        match run.evidence {
            RunEvidence::PhaseLoss => causes.push(Cause::PhaseLoss),
            RunEvidence::ReferenceFrame => causes.push(Cause::ReferenceFrame),
            RunEvidence::InPhase
            | RunEvidence::AnchoredEnd
            | RunEvidence::ConfirmedEnd
            | RunEvidence::TwoStrand => {}
        }
    }
    causes
}

/// The clean observers of a site by strand.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
struct Observers {
    forward: u8,
    reverse: u8,
}

/// Counts a site's clean observers (its supporting and opposing reads, by
/// content SHA-256) by strand.
fn observers(
    strands: &BTreeMap<&str, Strand>,
    supporting: &[String],
    opposing: &[String],
) -> Observers {
    let mut observed = Observers::default();
    for read in supporting.iter().chain(opposing) {
        match strands.get(read.as_str()) {
            Some(Strand::Forward) => observed.forward = observed.forward.saturating_add(1),
            Some(Strand::Reverse) => observed.reverse = observed.reverse.saturating_add(1),
            None => {}
        }
    }
    observed
}

/// The sort key of a trace, so traces sort the same whatever their order.
fn order(trace: &Trace) -> (&str, u8, &str, &[usize], Option<usize>) {
    let (rank, label, clean): (u8, &str, &[usize]) = match &trace.state {
        TraceState::Admitted { strand, clean } => (0, strand.label(), clean),
        TraceState::TooShort => (1, "", &[]),
        TraceState::Weak => (2, "", &[]),
        TraceState::Unplaced(label) => (3, label, &[]),
        TraceState::Rejected => (4, "", &[]),
    };
    (&trace.name, rank, label, clean, trace.span)
}

/// The ranges of `region.positions` that the segments cover.
fn overlaps(region: &Region, segments: &[ReferenceSegment]) -> Vec<std::ops::Range<usize>> {
    segments
        .iter()
        .map(|segment| region.overlap(segment.start_0based, segment.end_0based_exclusive))
        .filter(|range| !range.is_empty())
        .collect()
}

/// The primer a trace name carries: a token `HV<digits><F|R>` between
/// non-alphanumeric characters, in upper case (`D02_x_242285_HV1F_04.ab1` →
/// `HV1F`). This follows the lab's file naming and is a proxy only.
pub(crate) fn primer(name: &str) -> Option<String> {
    name.split(|symbol: char| !symbol.is_ascii_alphanumeric())
        .map(str::to_ascii_uppercase)
        .find(|token| {
            token
                .strip_prefix("HV")
                .and_then(|rest| rest.strip_suffix(['F', 'R']))
                .is_some_and(|digits| {
                    !digits.is_empty() && digits.bytes().all(|byte| byte.is_ascii_digit())
                })
        })
}

/// The strand a primer reads: `F` forward, `R` reverse.
pub(crate) fn primer_strand(primer: &str) -> Option<Strand> {
    match primer.as_bytes().last() {
        Some(b'F') => Some(Strand::Forward),
        Some(b'R') => Some(Strand::Reverse),
        _ => None,
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use dna_consensus::model::consensus::{
        ConsensusRun, ConsensusSegment, ConsensusSite, ConsensusSummary,
    };
    use dna_sample::model::sample_evidence::SampleReadAlignmentEvidence;

    use super::*;

    /// A placed read covering `covered` and reading `clean` cleanly, both
    /// 0-based half-open.
    pub(crate) fn read(
        name: &str,
        orientation: Orientation,
        covered: (usize, usize),
        clean: &[(usize, usize)],
    ) -> SampleReadEvidence {
        let segment = |(start, end): (usize, usize)| ReferenceSegment {
            start_0based: start,
            end_0based_exclusive: end,
        };
        SampleReadEvidence {
            input_name: name.into(),
            input_sha256: format!("sha-{name}"),
            alignment: SampleReadAlignmentEvidence {
                orientation,
                callable_bases: 0,
                identity: 1.0,
                gap_opens: 0,
                unresolved_bases: 0,
                masked_bases: 0,
                reference_segments: vec![segment(covered)],
                callable_reference_segments: clean.iter().copied().map(segment).collect(),
                wraps_origin: false,
            },
        }
    }

    pub(crate) fn site(
        first: usize,
        last: usize,
        call: Option<&str>,
        state: SiteState,
        supporting: &[&str],
        uninformative: &[&str],
        runs: &[RunEvidence],
    ) -> ConsensusSite {
        let names = |names: &[&str]| names.iter().map(|name| format!("sha-{name}")).collect();
        ConsensusSite {
            first_0based: first,
            last_0based: last,
            reference: "A".repeat(last - first + 1),
            call: call.map(str::to_owned),
            state,
            supporting: names(supporting),
            opposing: Vec::new(),
            uninformative: names(uninformative),
            runs: runs
                .iter()
                .map(|evidence| ConsensusRun {
                    base: 'C',
                    length: 3,
                    evidence: *evidence,
                })
                .collect(),
        }
    }

    pub(crate) fn consensus(sites: Vec<ConsensusSite>) -> Consensus {
        Consensus {
            segments: vec![ConsensusSegment {
                start_0based: 0,
                end_0based: 99,
                sequence: String::new(),
            }],
            sites,
            summary: ConsensusSummary::default(),
            reference_only: std::collections::BTreeSet::new(),
            unreliable_only: BTreeMap::new(),
        }
    }

    #[test]
    fn counts_clean_strands_masked_reads_and_consensus_marks() {
        // Region 11–30 (1-based) is 0-based 10..30.
        let reads = [
            read("s_HV1F_1.ab1", Orientation::Forward, (0, 30), &[(5, 25)]),
            read("s_HV1R_1.ab1", Orientation::Reverse, (15, 40), &[(20, 40)]),
        ];
        let sites = vec![
            site(
                21,
                21,
                None,
                SiteState::Unresolved,
                &[],
                &["s_HV1F_1.ab1"],
                &[],
            ),
            site(
                22,
                23,
                Some("CC"),
                SiteState::Contested,
                &["s_HV1R_1.ab1"],
                &[],
                &[RunEvidence::PhaseLoss],
            ),
        ];
        let spans = BTreeMap::from([("sha-s_HV1F_1.ab1", 400)]);
        let rejected = vec![Trace {
            name: "s_HV2R_1.ab1".into(),
            primer: primer("s_HV2R_1.ab1"),
            state: TraceState::Weak,
            span: Some(12),
        }];
        let evidence = build(&[[11, 30]], &reads, &consensus(sites), &spans, rejected);
        let region = &evidence.regions[0];
        assert_eq!(region.positions.len(), 20);
        // 1-based 11..=15: forward only (0-based 10..15).
        assert_eq!(
            (region.positions[0].forward, region.positions[0].reverse),
            (1, 0)
        );
        // 1-based 16..=20: forward clean, reverse covered through a masked call.
        assert_eq!(region.positions[5].masked_reads, 1);
        assert!(region.positions[5].has(Cause::MaskedRead));
        // 1-based 21: both strands.
        assert_eq!(region.positions[10].strands(), 2);
        // 1-based 22: the unresolved site: no clean observer, N and masked.
        let unresolved = region.positions[11];
        assert_eq!(unresolved.strands(), 0);
        assert!(unresolved.has(Cause::N) && unresolved.has(Cause::Masked));
        assert!(!unresolved.has(Cause::Contested));
        // At a site its uninformative reads are the masked ones.
        assert_eq!(unresolved.masked_reads, 1);
        // 1-based 23–24: contested, one reverse observer, phase loss.
        let contested = region.positions[12];
        assert_eq!(
            (contested.forward, contested.reverse, contested.masked_reads),
            (0, 1, 0)
        );
        assert!(contested.has(Cause::Contested) && contested.has(Cause::PhaseLoss));
        assert!(!contested.has(Cause::N) && !contested.has(Cause::ReferenceFrame));
        // 1-based 26..=30: reverse only.
        assert_eq!(region.positions[19].strands(), 1);

        assert_eq!(
            evidence
                .traces
                .iter()
                .map(|trace| trace.name.as_str())
                .collect::<Vec<_>>(),
            ["s_HV1F_1.ab1", "s_HV1R_1.ab1", "s_HV2R_1.ab1"]
        );
        assert_eq!(
            evidence.traces[0].state,
            TraceState::Admitted {
                strand: Strand::Forward,
                clean: vec![15]
            }
        );
        assert_eq!(evidence.traces[0].span, Some(400));
        assert_eq!(evidence.traces[1].span, None);
        assert_eq!(evidence.locate(11), Some((0, 0)));
        assert_eq!(evidence.locate(30), Some((0, 19)));
        assert_eq!(evidence.locate(31), None);
        assert_eq!(evidence.locate(-1), None);

        // The evidence does not depend on the order of the reads or traces.
        let reversed = [reads[1].clone(), reads[0].clone()];
        let rejected = vec![evidence.traces[2].clone()];
        let sites = consensus(vec![]);
        assert_eq!(
            build(&[[11, 30]], &reversed, &sites, &spans, rejected.clone()),
            build(&[[11, 30]], &reads, &sites, &spans, rejected)
        );
    }

    #[test]
    fn reads_of_one_name_count_by_their_own_strands() {
        let ids = |ids: &[&str]| ids.iter().map(|id| (*id).to_owned()).collect::<Vec<_>>();
        // Three reads that share a trace name, told apart by content.
        let strands = BTreeMap::from([
            ("a1", Strand::Forward),
            ("a2", Strand::Reverse),
            ("a3", Strand::Reverse),
        ]);
        assert_eq!(
            observers(&strands, &ids(&["a1", "a2"]), &ids(&["a3"])),
            Observers {
                forward: 1,
                reverse: 2,
            }
        );
    }

    #[test]
    fn reads_the_primer_from_the_trace_name() {
        assert_eq!(
            primer("D02_20240723_242285_HV1F_04.ab1").as_deref(),
            Some("HV1F")
        );
        assert_eq!(primer("x-hv3r.ab1").as_deref(), Some("HV3R"));
        assert_eq!(primer("HVF_1.ab1"), None);
        assert_eq!(primer("HV1X_1.ab1"), None);
        assert_eq!(primer("sample.ab1"), None);
        assert_eq!(primer_strand("HV2R"), Some(Strand::Reverse));
        assert_eq!(primer_strand("HV2F"), Some(Strand::Forward));
    }
}
