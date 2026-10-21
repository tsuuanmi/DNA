//! The labels of a decided site: reads that disagree across or within a
//! strand, a site resting on one strand (a run deletion, a run indel, an edge
//! deletion, a one-strand insertion), a long run, an undecided site, a call
//! beside a region end not analysed, and a length the profile's conventions
//! leave unrecorded.

use super::{Sample, runs_of};
use crate::model::review_result::{GrowthResult, LabelResult, Role, RunResult, SideResult};
use crate::sample_analysis::NotationCall;
use dna_align::model::alignment::{AlignmentColumn, Orientation};
use dna_consensus::InsertionView;
use dna_consensus::model::consensus::{ConsensusRun, ConsensusSite, RunEvidence, SiteState};
use dna_kernel::profile::{Convention, ConventionKind, ConventionLabel};

/// The label of a site's clean dissent, over every supporting and opposing
/// read, however many each strand has: `strands_disagree` when a supporting
/// read lies on another strand than an opposing one, otherwise
/// `same_strand_dissent` (all of them lie on one strand); `None` without a
/// supporting or an opposing read. Lists the supporting, then the opposing
/// reads.
pub(super) fn dissent(
    supporting: Vec<SideResult>,
    opposing: Vec<SideResult>,
) -> Option<LabelResult> {
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
pub(super) fn one_strand(supporting: &[Orientation]) -> Option<Orientation> {
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
    let reference = runs_of(&site.reference);
    let called = runs_of(call);
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
pub(super) fn beside_not_analysed(
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
    let (from, to) = (runs_of(reference), runs_of(call));
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
        && runs_of(reference).iter().any(|(_, length)| *length >= run)
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
pub(super) fn called_near(calls: &[NotationCall], first: usize, last: usize) -> bool {
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

impl Sample<'_> {
    /// The labels a consensus site carries by itself.
    pub(super) fn site_labels(&self, site: &ConsensusSite) -> Vec<LabelResult> {
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
            .filter(|run| {
                matches!(
                    run.evidence,
                    RunEvidence::PhaseLoss | RunEvidence::ReferenceFrame
                ) && !not_recorded(self.conventions, site, run)
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
    pub(super) fn one_strand_insertion(
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
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::report::review::test_support::{F, R};

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
}
