//! The consensus's decided edits as called records (PROP-0006,
//! SRS-API-012): each consensus segment becomes one called record without
//! being placed on the reference again.
//!
//! The segment's positions outside every site read the reference. A site's
//! call is aligned against its own reference interval only: an `N` pairs
//! with the reference base it stands for and is not callable, a substitution
//! pairs with its base, and the rest are insertions and deletions. Every
//! insertion and deletion then takes its 3'-most sequence-equivalent
//! placement over the segment, the move of ADR-0047 that read alignment
//! uses. The variants are then extracted and judged from those columns as
//! a record's are, under the same rules (ADR-0071, SRS-VAR-018): a length the
//! consensus could not bound, or a stretch it decided position by position,
//! stays excluded where a record of the same sequence would exclude it. What
//! changes is the placement: every base sits where the consensus decided it,
//! and an `N` stands only for the position the consensus left undecided. Two
//! rules are the segment's own, both excluding a change as `run_boundary`:
//! - a length change of the reference run at the segment's first or last
//!   position, since no read shows that run's outer end (ADR-0071);
//! - a change at the 3' end of a run whose length the consensus only
//!   estimated or assumed (`phase_loss`, `reference_frame`), where a change
//!   of its length is written; that position is not callable either. The
//!   run's other positions, and so its composition, stay callable
//!   (PROP-0008 step 3).

use crate::model::alignment::{Alignment, AlignmentColumn, ReferenceSegment};
use crate::model::called_read::CalledRead;
use std::collections::BTreeSet;
use std::ops::Range;

use crate::model::consensus::{Consensus, ConsensusSegment, RunEvidence};
use crate::model::variant::{Variant, VariantExclusionReason, VariantKind};
use crate::read_call::{CoreRun, ReadIdentity};
use crate::variant_calling;
use dna_kernel::error::{Result, VariantError};
use dna_kernel::model::nucleotide::is_canonical;
use dna_kernel::model::placement;
use dna_kernel::model::reference::Reference;
use dna_kernel::plugin::{Contract, PluginDescriptor, PluginFamily};
use dna_kernel::read_evidence::ReadEvidence;

/// Calls the consensus segment at `index` from the consensus's decisions.
/// `evidence` is the segment sequence's evidence, as a record of it reads,
/// and `identity` names the record.
///
/// # Errors
///
/// Returns `VariantError` when the segment sequence does not match its
/// sites and the reference, or when a decided edit cannot be anchored.
pub fn call_segment(
    consensus: &Consensus,
    index: usize,
    identity: ReadIdentity,
    evidence: ReadEvidence,
    reference: &Reference,
    run: &CoreRun<'_>,
) -> Result<CalledRead> {
    let segment = consensus
        .segments
        .get(index)
        .ok_or(VariantError::Inconsistent("no consensus segment at index"))?;
    let (columns, estimated) = columns(consensus, segment, reference.sequence.as_bytes())?;
    let columns = right_aligned(columns);
    let cut = estimated
        .iter()
        .filter_map(|range| run_end(&columns, range))
        .collect::<BTreeSet<_>>();
    let callable_segments = callable(&columns)
        .into_iter()
        .flat_map(|segment| {
            let mut pieces = Vec::new();
            let mut start = segment.start_0based;
            for &position in cut.range(segment.start_0based..segment.end_0based_exclusive) {
                if position > start {
                    pieces.push(ReferenceSegment {
                        start_0based: start,
                        end_0based_exclusive: position,
                    });
                }
                start = position + 1;
            }
            if segment.end_0based_exclusive > start {
                pieces.push(ReferenceSegment {
                    start_0based: start,
                    end_0based_exclusive: segment.end_0based_exclusive,
                });
            }
            pieces
        })
        .collect();
    let alignment = Alignment::placed(
        vec![ReferenceSegment {
            start_0based: segment.start_0based,
            end_0based_exclusive: segment.end_0based + 1,
        }],
        callable_segments,
        columns,
    )?;
    let mut variants = variant_calling::call(
        &alignment,
        reference,
        &evidence,
        &run.config.variant_calling,
        run.regions,
    )?;
    // A length change of the reference run at either end of the segment: no
    // read shows that run's outer end, since the segment ends where the reads
    // do (ADR-0071).
    let edges = edge_runs(segment, reference.sequence.as_bytes());
    for observed in &mut variants.observed {
        let variant = &observed.variant;
        let edge =
            variant.kind != VariantKind::Snv && edges.iter().any(|run| touches(variant, *run));
        // A change of the length of a run the consensus only estimated or
        // assumed (`phase_loss`, `reference_frame`): PROP-0008 step 3.
        let estimate = changed(variant)
            .iter()
            .any(|position| cut.contains(position));
        if (edge || estimate)
            && !observed
                .exclusion_reasons
                .contains(&VariantExclusionReason::RunBoundary)
        {
            observed
                .exclusion_reasons
                .push(VariantExclusionReason::RunBoundary);
        }
    }
    variants.reported = variants
        .observed
        .iter()
        .filter(|observed| observed.eligible())
        .map(|observed| observed.variant.clone())
        .collect();
    Ok(CalledRead {
        input_name: identity.input_name,
        input_sha256: identity.input_sha256,
        reference_sha256: reference.sequence_sha256.clone(),
        configuration_sha256: run.configuration_sha256.to_owned(),
        evidence,
        alignment,
        variants,
    })
}

/// The 3'-most reference position of the run whose calls are `range`, the
/// deleted positions beside it included: where a 3'-normalized change of the
/// run's length is written. A run inserted whole ends at the position before
/// it.
fn run_end(columns: &[AlignmentColumn], range: &Range<usize>) -> Option<usize> {
    let marked = |column: &AlignmentColumn| {
        column
            .original_call_index_0based()
            .is_some_and(|call| range.contains(&call))
    };
    let mut end = columns.iter().rposition(marked)?;
    while end + 1 < columns.len() && columns[end + 1].query_base() == '-' {
        end += 1;
    }
    columns[..=end]
        .iter()
        .rev()
        .find_map(AlignmentColumn::reference_index_0based)
}

/// The reference positions a variant changes: a substitution's position, a
/// deletion's deleted positions, both positions beside an insertion.
fn changed(variant: &Variant) -> Vec<usize> {
    let anchor = variant.position_1based.saturating_sub(1);
    match variant.kind {
        VariantKind::Snv => vec![anchor],
        VariantKind::Del => (anchor + 1..anchor + variant.reference.len()).collect(),
        VariantKind::Ins => vec![anchor, anchor + 1],
    }
}

/// The reference runs holding the segment's first and last positions,
/// 0-based inclusive.
fn edge_runs(segment: &ConsensusSegment, bases: &[u8]) -> [(usize, usize); 2] {
    let run = |position: usize| {
        let base = bases[position];
        let mut first = position;
        while first > segment.start_0based && bases[first - 1] == base {
            first -= 1;
        }
        let mut last = position;
        while last < segment.end_0based && bases[last + 1] == base {
            last += 1;
        }
        (first, last)
    };
    [run(segment.start_0based), run(segment.end_0based)]
}

/// Whether a variant changes the length of the 0-based inclusive `run`: a
/// deletion of one of its positions, or an insertion at a junction beside
/// one.
fn touches(variant: &Variant, (first, last): (usize, usize)) -> bool {
    let anchor = variant.position_1based.saturating_sub(1);
    let (start, end) = match variant.kind {
        VariantKind::Del => (
            anchor + 1,
            anchor + variant.reference.len().saturating_sub(1),
        ),
        _ => (anchor, anchor + 1),
    };
    start <= last && first <= end
}

/// The segment's alignment columns: reference bases outside the sites, each
/// site's call aligned against its interval, and the segment sequence's
/// call indices, which must account for every base of it. Also the
/// call-index ranges of the runs whose length the consensus only estimated
/// or assumed.
fn columns(
    consensus: &Consensus,
    segment: &ConsensusSegment,
    bases: &[u8],
) -> Result<(Vec<AlignmentColumn>, Vec<Range<usize>>)> {
    let sequence = segment.sequence.as_bytes();
    let mismatch =
        || VariantError::Inconsistent("consensus segment sequence does not match its sites");
    let mut sites = consensus
        .sites
        .iter()
        .filter(|site| {
            segment.start_0based <= site.first_0based && site.last_0based <= segment.end_0based
        })
        .collect::<Vec<_>>();
    sites.sort_by_key(|site| site.first_0based);
    // (whether a call is read, the reference position), in order.
    let mut steps: Vec<(bool, Option<usize>)> = Vec::with_capacity(sequence.len());
    let mut estimated = Vec::new();
    let mut position = segment.start_0based;
    for site in sites {
        if site.first_0based < position {
            return Err(mismatch().into());
        }
        steps.extend((position..site.first_0based).map(|reference| (true, Some(reference))));
        let interval = bases
            .get(site.first_0based..=site.last_0based)
            .ok_or_else(mismatch)?;
        let decided = site
            .call
            .clone()
            .unwrap_or_else(|| "N".repeat(interval.len()));
        let mut start = steps.iter().filter(|(reads, _)| *reads).count();
        for run in &site.runs {
            if matches!(
                run.evidence,
                RunEvidence::PhaseLoss | RunEvidence::ReferenceFrame
            ) {
                estimated.push(start..start + run.length);
            }
            start += run.length;
        }
        steps.extend(align(decided.as_bytes(), interval).into_iter().map(
            |operation| match operation {
                Operation::Pair(offset) => (true, Some(site.first_0based + offset)),
                Operation::Inserted => (true, None),
                Operation::Deleted(offset) => (false, Some(site.first_0based + offset)),
            },
        ));
        position = site.last_0based + 1;
    }
    steps.extend((position..=segment.end_0based).map(|reference| (true, Some(reference))));
    let mut call = 0;
    let mut columns = Vec::with_capacity(steps.len());
    for (reads, reference) in steps {
        let query = if reads {
            let base = *sequence.get(call).ok_or_else(mismatch)?;
            call += 1;
            Some((call - 1, char::from(base)))
        } else {
            None
        };
        let reference_base = match reference {
            Some(position) => char::from(*bases.get(position).ok_or_else(mismatch)?),
            None => '-',
        };
        columns.push(AlignmentColumn::new(
            query.map_or('-', |(_, base)| base),
            reference_base,
            query.map(|(index, _)| index),
            reference,
        ));
    }
    if call != sequence.len() {
        return Err(mismatch().into());
    }
    Ok((columns, estimated))
}

/// One column of a site's alignment: a call paired with a reference offset,
/// an inserted call, or a deleted reference offset.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Operation {
    Pair(usize),
    Inserted,
    Deleted(usize),
}

/// Cost of a substitution; an `N` pairs with any base at no cost.
const MISMATCH: u32 = 2;
/// Cost of opening and of extending an insertion or a deletion. A gap of two
/// bases plus a substitution costs more than two one-base gaps, so a length
/// change in each of two runs stays two insertions or deletions rather than
/// one gap across the runs' boundary and a substitution; substitutions stay
/// cheaper than an insertion plus a deletion of as many bases.
const GAP_OPEN: u32 = 3;
const GAP_EXTEND: u32 = 2;

/// The least-cost global alignment of a site's `call` against its reference
/// `interval`, affine in gaps. Traced back from the 3' end, a pair comes
/// before a deletion, then an insertion, of equal cost; the placement among
/// sequence-equivalent gaps is settled afterwards by the 3'-most move.
fn align(call: &[u8], interval: &[u8]) -> Vec<Operation> {
    #[derive(Clone, Copy, PartialEq, Eq)]
    enum State {
        Pair,
        Deleted,
        Inserted,
    }
    const UNREACHABLE: u32 = u32::MAX / 4;
    let (rows, columns) = (call.len() + 1, interval.len() + 1);
    // Per state, the best cost of a path ending in it and its predecessor.
    let mut cost = vec![[UNREACHABLE; 3]; rows * columns];
    let mut from = vec![[State::Pair; 3]; rows * columns];
    let at = |row: usize, column: usize| row * columns + column;
    let best = |costs: [u32; 3], extend: [u32; 3]| {
        [State::Pair, State::Deleted, State::Inserted]
            .into_iter()
            .zip(costs.iter().zip(extend))
            .map(|(state, (cost, step))| (cost.saturating_add(step), state))
            .min_by_key(|(cost, _)| *cost)
            .unwrap_or((UNREACHABLE, State::Pair))
    };
    cost[at(0, 0)][State::Pair as usize] = 0;
    for row in 0..rows {
        for column in 0..columns {
            if row > 0 && column > 0 {
                let step = if call[row - 1] == b'N' || call[row - 1] == interval[column - 1] {
                    0
                } else {
                    MISMATCH
                };
                let (value, state) = best(cost[at(row - 1, column - 1)], [step; 3]);
                cost[at(row, column)][State::Pair as usize] = value;
                from[at(row, column)][State::Pair as usize] = state;
            }
            if column > 0 {
                let (value, state) =
                    best(cost[at(row, column - 1)], [GAP_OPEN, GAP_EXTEND, GAP_OPEN]);
                cost[at(row, column)][State::Deleted as usize] = value;
                from[at(row, column)][State::Deleted as usize] = state;
            }
            if row > 0 {
                let (value, state) =
                    best(cost[at(row - 1, column)], [GAP_OPEN, GAP_OPEN, GAP_EXTEND]);
                cost[at(row, column)][State::Inserted as usize] = value;
                from[at(row, column)][State::Inserted as usize] = state;
            }
        }
    }
    let (mut row, mut column) = (call.len(), interval.len());
    let (_, mut state) = best(cost[at(row, column)], [0; 3]);
    let mut operations = Vec::with_capacity(rows + columns);
    while row > 0 || column > 0 {
        let previous = from[at(row, column)][state as usize];
        match state {
            State::Pair => {
                operations.push(Operation::Pair(column - 1));
                row -= 1;
                column -= 1;
            }
            State::Deleted => {
                operations.push(Operation::Deleted(column - 1));
                column -= 1;
            }
            State::Inserted => {
                operations.push(Operation::Inserted);
                row -= 1;
            }
        }
        state = previous;
    }
    operations.reverse();
    operations
}

/// Moves every insertion and deletion to its 3'-most sequence-equivalent
/// placement (ADR-0047): across a following matched column whose base the
/// gap's first base equals, never across an `N`, a substitution, or into
/// another gap.
fn right_aligned(mut columns: Vec<AlignmentColumn>) -> Vec<AlignmentColumn> {
    loop {
        let mut moved = false;
        let mut start = 0;
        while start < columns.len() {
            let kind = gap_kind(&columns[start]);
            if kind.is_none() {
                start += 1;
                continue;
            }
            let mut end = start;
            while end < columns.len() && gap_kind(&columns[end]) == kind {
                end += 1;
            }
            if shift(&mut columns, start, end) {
                moved = true;
                start += 1;
            } else {
                start = end;
            }
        }
        if !moved {
            return columns;
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Gap {
    Inserted,
    Deleted,
}

fn gap_kind(column: &AlignmentColumn) -> Option<Gap> {
    match (column.query_base(), column.reference_base()) {
        ('-', _) => Some(Gap::Deleted),
        (_, '-') => Some(Gap::Inserted),
        _ => None,
    }
}

/// Moves the gap `start..end` one base to the right when the move keeps the
/// sequence: the following column is a match of canonical bases equal to the
/// gap's first base, and the column after it is no gap of the same kind.
fn shift(columns: &mut Vec<AlignmentColumn>, start: usize, end: usize) -> bool {
    let kind = gap_kind(&columns[start]);
    let Some(following) = columns.get(end).cloned() else {
        return false;
    };
    let matched = following.query_base() == following.reference_base()
        && is_canonical(following.query_base());
    let joins = columns
        .get(end + 1)
        .is_some_and(|after| gap_kind(after) == kind);
    let first = columns[start].clone();
    let leading = match kind {
        Some(Gap::Inserted) => first.query_base(),
        Some(Gap::Deleted) => first.reference_base(),
        None => return false,
    };
    if !matched || joins || !placement::crosses(leading as u8, following.reference_base() as u8) {
        return false;
    }
    let gap = columns[start + 1..end].to_vec();
    let moved = match kind {
        // The first inserted call pairs with the crossed reference base, and
        // the crossed call becomes the last inserted one.
        Some(Gap::Inserted) => std::iter::once(AlignmentColumn::new(
            first.query_base(),
            following.reference_base(),
            first.original_call_index_0based(),
            following.reference_index_0based(),
        ))
        .chain(gap)
        .chain(std::iter::once(AlignmentColumn::new(
            following.query_base(),
            '-',
            following.original_call_index_0based(),
            None,
        )))
        .collect::<Vec<_>>(),
        // The first deleted reference base pairs with the crossed call, and
        // the crossed reference base becomes the last deleted one.
        _ => std::iter::once(AlignmentColumn::new(
            following.query_base(),
            first.reference_base(),
            following.original_call_index_0based(),
            first.reference_index_0based(),
        ))
        .chain(gap)
        .chain(std::iter::once(AlignmentColumn::new(
            '-',
            following.reference_base(),
            None,
            following.reference_index_0based(),
        )))
        .collect::<Vec<_>>(),
    };
    columns.splice(start..=end, moved);
    true
}

/// The reference positions the segment calls: every position paired with a
/// canonical base, and every deleted one; an `N` or another code is not
/// callable.
fn callable(columns: &[AlignmentColumn]) -> Vec<ReferenceSegment> {
    let mut segments: Vec<ReferenceSegment> = Vec::new();
    for column in columns {
        let Some(position) = column.reference_index_0based() else {
            continue;
        };
        if column.query_base() != '-' && !is_canonical(column.query_base()) {
            continue;
        }
        match segments.last_mut() {
            Some(last) if last.end_0based_exclusive == position => {
                last.end_0based_exclusive = position + 1;
            }
            _ => segments.push(ReferenceSegment {
                start_0based: position,
                end_0based_exclusive: position + 1,
            }),
        }
    }
    segments
}

/// The consensus's decided edits as called records.
pub const PLUGIN: PluginDescriptor = PluginDescriptor {
    id: "decided",
    family: PluginFamily::Core,
    version: 2,
    provides: &[Contract::CalledVariants],
    requires: &[Contract::Consensus],
    config_sections: &[],
};

#[cfg(test)]
mod tests {
    use dna_kernel::model::reference::ReferenceTopology;

    use super::*;
    use crate::alignment::RawAlignmentConfig;
    use crate::model::consensus::{
        ConsensusRun, ConsensusSite, ConsensusSummary, RunEvidence, SiteState,
    };
    use crate::read_call::{CallConfig, RawCallConfig};

    /// `ACGA TTTTT CAGGCTAGCATCG` at 1..=22: a T run at 5..=9.
    const REFERENCE: &str = "ACGATTTTTCAGGCTAGCATCG";

    fn reference() -> Reference {
        Reference {
            name: "ref".into(),
            sequence: REFERENCE.into(),
            topology: ReferenceTopology::Linear,
            sequence_sha256: String::new(),
        }
    }

    fn config() -> Result<CallConfig> {
        RawCallConfig {
            alignment: RawAlignmentConfig {
                match_score: 3,
                mismatch_score: -5,
                ambiguous_score: 0,
                gap_open_score: -10,
                gap_extension_score: -4,
                minimum_callable_bases: 1,
                minimum_identity: 0.8,
            },
            variant_calling: variant_calling::RawVariantCallingConfig {
                max_indel_length: 50,
                read_end_margin: 0,
            },
        }
        .validate()
    }

    /// A site over 1-based `first..=last` decided as `call`.
    fn site(
        first: usize,
        last: usize,
        call: Option<&str>,
        state: SiteState,
        runs: Vec<ConsensusRun>,
    ) -> ConsensusSite {
        ConsensusSite {
            first_0based: first - 1,
            last_0based: last - 1,
            reference: REFERENCE[first - 1..last].into(),
            call: call.map(Into::into),
            state,
            supporting: Vec::new(),
            opposing: Vec::new(),
            uninformative: Vec::new(),
            runs,
        }
    }

    /// The segment over 1-based `first..=last` with `sequence`, and its sites.
    fn called(
        first: usize,
        last: usize,
        sequence: &str,
        sites: Vec<ConsensusSite>,
    ) -> Result<CalledRead> {
        let consensus = Consensus {
            segments: vec![ConsensusSegment {
                start_0based: first - 1,
                end_0based: last - 1,
                sequence: sequence.into(),
            }],
            sites,
            summary: ConsensusSummary::default(),
        };
        let config = config()?;
        call_segment(
            &consensus,
            0,
            ReadIdentity {
                input_name: "segment".into(),
                input_sha256: "0".into(),
                declared: None,
            },
            ReadEvidence::clean(sequence),
            &reference(),
            &CoreRun {
                config: &config,
                configuration_sha256: "0",
                regions: &[[1, 22]],
            },
        )
    }

    fn edits(read: &CalledRead) -> Vec<(usize, String, String, Vec<&'static str>)> {
        read.variants
            .observed
            .iter()
            .map(|observed| {
                (
                    observed.variant.position_1based,
                    observed.variant.reference.clone(),
                    observed.variant.alternate.clone(),
                    observed
                        .exclusion_reasons
                        .iter()
                        .map(|reason| reason.label())
                        .collect(),
                )
            })
            .collect()
    }

    fn callable(read: &CalledRead) -> Vec<(usize, usize)> {
        read.alignment
            .callable_segments
            .iter()
            .map(|segment| (segment.start_0based, segment.end_0based_exclusive))
            .collect()
    }

    #[test]
    fn a_run_length_change_takes_the_runs_3_prime_end() -> Result<()> {
        let read = called(
            1,
            22,
            "ACGATTTTTTCAGGCTAGCATCG",
            vec![site(
                5,
                9,
                Some("TTTTTT"),
                SiteState::Called,
                vec![ConsensusRun {
                    base: 'T',
                    length: 6,
                    evidence: RunEvidence::InPhase,
                }],
            )],
        )?;
        assert_eq!(edits(&read), [(9, "T".into(), "TT".into(), Vec::new())]);
        assert_eq!(callable(&read), [(0, 22)]);
        Ok(())
    }

    #[test]
    fn an_undecided_site_is_not_callable_and_a_decided_base_is_reported() -> Result<()> {
        // 12A decided; 15 undecided, written N; a resolved tie at 18 keeps the
        // reference.
        let read = called(
            1,
            22,
            "ACGATTTTTCAAGCNAGCATCG",
            vec![
                site(12, 12, Some("A"), SiteState::Called, Vec::new()),
                site(15, 15, None, SiteState::Unresolved, Vec::new()),
                site(18, 18, Some("C"), SiteState::Contested, Vec::new()),
            ],
        )?;
        assert_eq!(edits(&read), [(12, "G".into(), "A".into(), Vec::new())]);
        assert_eq!(callable(&read), [(0, 14), (15, 22)]);
        Ok(())
    }

    #[test]
    fn a_length_change_of_the_run_at_a_segment_end_is_unbounded() -> Result<()> {
        // The segment starts inside the T run: no read shows its 5' end.
        let read = called(
            5,
            22,
            "TTTTCAGGCTAGCATCG",
            vec![site(
                5,
                9,
                Some("TTTT"),
                SiteState::Called,
                vec![ConsensusRun {
                    base: 'T',
                    length: 4,
                    evidence: RunEvidence::InPhase,
                }],
            )],
        )?;
        assert_eq!(
            edits(&read),
            [(8, "TT".into(), "T".into(), vec!["run_boundary"])]
        );
        Ok(())
    }

    #[test]
    fn an_estimated_run_length_is_not_callable() -> Result<()> {
        // A sixth T the consensus only estimated: the run's 3' end is not
        // callable, and the insertion there is excluded.
        let read = called(
            1,
            22,
            "ACGATTTTTTCAGGCTAGCATCG",
            vec![site(
                5,
                9,
                Some("TTTTTT"),
                SiteState::Called,
                vec![ConsensusRun {
                    base: 'T',
                    length: 6,
                    evidence: RunEvidence::PhaseLoss,
                }],
            )],
        )?;
        assert_eq!(
            edits(&read),
            [(9, "T".into(), "TT".into(), vec!["run_boundary"])]
        );
        assert_eq!(callable(&read), [(0, 8), (9, 22)]);
        assert!(read.variants.reported.is_empty());
        Ok(())
    }

    #[test]
    fn a_decided_deletion_between_undecided_positions_is_not_callable_reference() -> Result<()> {
        // 11 and 14 undecided, 12-13 (GG) decided deleted: the deletion's gap
        // could lie anywhere among the `N`s, so it is excluded, and its
        // positions must not read as the reference (ADR-0078, SRS-VAR-019).
        let read = called(
            1,
            22,
            "ACGATTTTTCNNTAGCATCG",
            vec![
                site(11, 11, None, SiteState::Unresolved, Vec::new()),
                site(12, 13, Some(""), SiteState::Called, Vec::new()),
                site(14, 14, None, SiteState::Unresolved, Vec::new()),
            ],
        )?
        .without_excluded_positions();
        assert!(read.variants.reported.is_empty());
        assert!(
            callable(&read)
                .iter()
                .all(|&(start, end)| !(start..end).contains(&11) && !(start..end).contains(&12))
        );
        Ok(())
    }

    #[test]
    fn a_sequence_that_disagrees_with_its_sites_fails() {
        assert!(called(1, 22, "ACGATTTTTCAGGCTAGCATC", Vec::new()).is_err());
    }
}
