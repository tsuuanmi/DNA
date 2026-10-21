//! What each read shows over an item: its labels and evidence there (one read,
//! one strand, noisy signal, phase shift, read edge, poor signal), the
//! sequence it shows between two clean flanks, whether it disagrees with the
//! decision, and the stretches where it is masked or unreliable.

use super::Sample;
use super::ranges::inside;
use super::round;
use super::sites::dissent;
use crate::error::{ReportError, Result};
use crate::model::result::{CallabilitySegmentResult, PhaseStateResult, ScopeResult};
use crate::model::review_result::{
    EvidencePositionResult, EvidenceReadResult, LabelResult, MaskedSpanResult, Role, SecondResult,
    UnreliableSpanResult,
};
use crate::model::sanger_read::SangerRead;
use crate::sample_analysis::{CallTrust, TrustScope};
use dna_align::model::alignment::{AlignmentColumn, Orientation};
use dna_consensus::model::consensus::ConsensusSite;
use dna_kernel::read_evidence::{Placement, ReadEvidence};
use dna_sample::model::sample_evidence::{LocusState, SampleLocusObservation};
use std::collections::{BTreeMap, BTreeSet};

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
        "ACGT".contains(column.query_base())
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

/// The reference position (1-based) a read's call `call` is aligned to, or,
/// on an inserted column, the position the insertion follows on the
/// reference, on either strand: an inserted base of a run lies at the run's
/// 3' end (SRS-ALN-012), after the run's last position. The first position
/// after it stands in only where no aligned column precedes it.
pub(super) fn start_reference(columns: &[AlignmentColumn], call: usize) -> Option<usize> {
    let index = columns
        .iter()
        .position(|column| column.original_call_index_0based() == Some(call))?;
    let aligned = |column: &AlignmentColumn| column.reference_index_0based();
    columns[..=index]
        .iter()
        .rev()
        .find_map(aligned)
        .or_else(|| columns[index..].iter().find_map(aligned))
        .map(|position| position + 1)
}

/// The reference positions, within `regions`, that a read's aligned calls
/// cover inside one of its unreliable spans (SRS-CALL-024), as maximal
/// stretches of one reason.
pub(super) fn unreliable_spans(
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
pub(super) fn masked_spans(
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
        if segment.state == PhaseStateResult::InPhase
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

impl Sample<'_> {
    /// The labels of the reads' calls over `first..=last` that are not
    /// clean, one per label and read, and every covering read's evidence.
    pub(super) fn reads_over(
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
    pub(super) fn clipped_decision(
        &self,
        site: &ConsensusSite,
        first: usize,
        last: usize,
    ) -> Option<String> {
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
    pub(super) fn clip_disagreement(
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
    pub(super) fn disagreeing(
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
    pub(super) fn end_disagreeing(
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

    pub(super) fn observation(
        &self,
        index: usize,
        position: usize,
    ) -> Option<&SampleLocusObservation> {
        self.loci
            .get(&position)?
            .observations
            .iter()
            .find(|observation| observation.read_index == index)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use dna_kernel::read_evidence::EvidenceReason;

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

    // SRS-API-020
    #[test]
    fn an_inserted_first_call_is_credited_to_the_position_it_follows() {
        // Reference 315-316 (0-based 314-315) read C G; a reverse read's run
        // of C starts, in its order, on the C inserted after 315.
        let columns = [
            column('C', Some(9), Some(313)),
            column('C', Some(8), Some(314)),
            column('C', Some(7), None),
            column('G', Some(6), Some(315)),
        ];
        assert_eq!(start_reference(&columns, 7), Some(315));
        assert_eq!(start_reference(&columns, 6), Some(316));
        // With no aligned column before it, the next one stands in.
        let leading = [column('C', Some(0), None), column('C', Some(1), Some(4))];
        assert_eq!(start_reference(&leading, 0), Some(5));
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
}
