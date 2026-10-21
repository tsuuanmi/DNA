//! Forward/reverse profile alignment selection and circular coordinate projection.

use std::cmp::Ordering;
use std::ops::Range;

use crate::AlignmentConfig;
use crate::canonical;
use crate::exact::{self, Problem, Pruned, UpperBoundPlacement, WorkingReference};
use crate::gotoh;
use crate::model::alignment::{Alignment, AlignmentColumn, Orientation, ReferenceSegment};
use crate::runs;
use crate::traceback::RawAlignment;
use dna_kernel::error::{AlignmentError, Result};
use dna_kernel::model::nucleotide::{is_canonical, reverse_complement};
use dna_kernel::model::reference::{Reference, ReferenceTopology};
use dna_kernel::read_evidence::{CallEvidence, EvidenceProfile, Placement, ReadEnds, ReadEvidence};

struct Candidate {
    orientation: Orientation,
    mapping: Vec<usize>,
    placements: Vec<RawAlignment>,
}

/// The selected orientation's outcome: one placement, or two equally scoring
/// ones and the calls they place differently.
enum Selection {
    Unique(Alignment),
    Tied(Vec<usize>),
}

/// Aligns both evidence-profile orientations and returns one unique selected result.
///
/// The query is the read's informative calls. An unresolved masked call enters
/// it as `N` without an evidence profile, so it neither scores nor counts as
/// callable; an anchoring masked call keeps its base and profile.
///
/// When the selected orientation's two best placements tie, and every call
/// from one end of the query to the farthest call they place differently is
/// unreliable — inside one of the read's unreliable spans, or masked — those
/// calls are left unplaced and the rest of the query is aligned again
/// (SRS-ALN-006).
///
/// # Errors
///
/// Returns `AlignmentError` when no unique placement meets the configured
/// minimums.
pub fn align_best(
    evidence: &ReadEvidence,
    reference: &Reference,
    config: &AlignmentConfig,
) -> Result<Alignment> {
    let informative = evidence.informative();
    let disputed = match align_calls(evidence, informative.clone(), reference, config)? {
        Selection::Unique(alignment) => return Ok(alignment),
        Selection::Tied(disputed) => disputed,
    };
    let kept =
        placed_calls(evidence, informative, &disputed).ok_or(AlignmentError::AmbiguousPlacement)?;
    match align_calls(evidence, kept, reference, config)? {
        Selection::Unique(alignment) => Ok(alignment),
        Selection::Tied(_) => Err(AlignmentError::AmbiguousPlacement.into()),
    }
}

/// The calls a tie leaves placed: the query without the stretch from one of
/// its ends to the farthest of the `disputed` calls, when every call of that
/// stretch is unreliable (inside one of the read's unreliable spans, or
/// masked) and some call remains. `None` otherwise.
fn placed_calls(
    evidence: &ReadEvidence,
    query: Range<usize>,
    disputed: &[usize],
) -> Option<Range<usize>> {
    let unreliable = |call: usize| {
        evidence
            .unreliable()
            .iter()
            .any(|span| span.calls.contains(&call))
            || evidence
                .calls()
                .get(call)
                .is_some_and(|call| call.placement != Placement::Scored)
    };
    let first = *disputed.iter().min()?;
    let last = *disputed.iter().max()?;
    if (query.start..=last).all(unreliable) && last + 1 < query.end {
        return Some(last + 1..query.end);
    }
    if (first..query.end).all(unreliable) && query.start < first {
        return Some(query.start..first);
    }
    None
}

/// Aligns the `query` calls in both orientations and selects one.
fn align_calls(
    evidence: &ReadEvidence,
    informative: Range<usize>,
    reference: &Reference,
    config: &AlignmentConfig,
) -> Result<Selection> {
    let calls = evidence.calls();
    let retained = &calls[informative.clone()];
    let forward_query = retained
        .iter()
        .map(|call| {
            if unresolved(Some(call.placement)) {
                'N'
            } else {
                call.base
            }
        })
        .collect::<String>();
    let reverse_query = reverse_complement(&forward_query);
    let forward_profiles = retained
        .iter()
        .map(|call| {
            if unresolved(Some(call.placement)) {
                None
            } else {
                call.profile
            }
        })
        .collect::<Vec<_>>();
    let reverse_profiles = reverse_profiles(&forward_profiles);
    let forward_mapping = informative.clone().collect::<Vec<_>>();
    let reverse_mapping = informative.rev().collect::<Vec<_>>();
    let masked = calls;
    let vouched = evidence.ends() == ReadEnds::Vouched;
    let (working_reference, modulo_length) = match reference.topology {
        ReferenceTopology::Linear => (reference.sequence.clone(), None),
        ReferenceTopology::Circular => (
            format!("{}{}", reference.sequence, reference.sequence),
            Some(reference.len()),
        ),
    };
    let target = WorkingReference::new(&working_reference, modulo_length);
    let forward = Query::new(
        Orientation::Forward,
        &forward_query,
        &forward_profiles,
        forward_mapping,
        &target,
        config,
    );
    let reverse = Query::new(
        Orientation::Reverse,
        &reverse_query,
        &reverse_profiles,
        reverse_mapping,
        &target,
        config,
    );
    if let Some(selected) = select_proven_candidate(
        forward.classify(),
        reverse.classify(),
        &forward.mapping,
        &reverse.mapping,
    )? {
        return finish_alignment(&selected, reference, config, masked, vouched);
    }

    // Evaluate one orientation exactly, preferring a proven one, then the one
    // with the better seed lower bound; its optimum is a threshold the other
    // orientation must reach.
    let (first, first_placements, second, second_proven) =
        match (forward.pruned()?, reverse.pruned()?) {
            (Pruned::Proven(placements), other) => {
                (forward, placements, reverse, other.into_proven())
            }
            (Pruned::Unproven { .. }, Pruned::Proven(placements)) => {
                (reverse, placements, forward, None)
            }
            (
                Pruned::Unproven {
                    lower_bound: forward_lower,
                },
                Pruned::Unproven {
                    lower_bound: reverse_lower,
                },
            ) => {
                let (first, second) = if reverse_lower > forward_lower {
                    (reverse, forward)
                } else {
                    (forward, reverse)
                };
                let placements = first.gotoh()?;
                (first, placements, second, None)
            }
        };
    let first = first.candidate(first_placements);
    let second_placements = match second_proven {
        Some(placements) => placements,
        None => match second.at_or_above(first.placements[0].score)? {
            Some(placements) if placements.is_empty() => {
                return finish_alignment(&first, reference, config, masked, vouched);
            }
            Some(placements) => placements,
            None => second.gotoh()?,
        },
    };
    let second = second.candidate(second_placements);
    let (forward, reverse) = if first.orientation == Orientation::Forward {
        (first, second)
    } else {
        (second, first)
    };
    let ordering = compare(&forward.placements[0], &reverse.placements[0]);
    let selected = match ordering {
        Ordering::Greater => &forward,
        Ordering::Less => &reverse,
        Ordering::Equal => {
            return Err(AlignmentError::OrientationTie.into());
        }
    };
    finish_alignment(selected, reference, config, masked, vouched)
}

/// One orientation's alignment query, with its proof problem when one applies.
struct Query<'a> {
    orientation: Orientation,
    bases: &'a str,
    profiles: &'a [Option<EvidenceProfile>],
    mapping: Vec<usize>,
    reference: &'a WorkingReference<'a>,
    config: &'a AlignmentConfig,
    problem: Option<Problem<'a>>,
}

impl<'a> Query<'a> {
    fn new(
        orientation: Orientation,
        bases: &'a str,
        profiles: &'a [Option<EvidenceProfile>],
        mapping: Vec<usize>,
        reference: &'a WorkingReference<'a>,
        config: &'a AlignmentConfig,
    ) -> Self {
        Self {
            orientation,
            bases,
            profiles,
            mapping,
            reference,
            config,
            problem: Problem::new(bases, profiles, reference, config),
        }
    }

    fn classify(&self) -> UpperBoundPlacement {
        self.problem
            .as_ref()
            .map_or(UpperBoundPlacement::Unproven, exact::classify)
    }

    fn pruned(&self) -> Result<Pruned> {
        match &self.problem {
            Some(problem) => exact::align_pruned(problem),
            None => Ok(Pruned::Unproven { lower_bound: None }),
        }
    }

    fn at_or_above(&self, threshold: i64) -> Result<Option<Vec<RawAlignment>>> {
        match &self.problem {
            Some(problem) => exact::align_at_or_above(problem, threshold),
            None => Ok(None),
        }
    }

    fn gotoh(&self) -> Result<Vec<RawAlignment>> {
        gotoh::align(
            self.bases,
            self.profiles,
            self.reference.sequence(),
            self.config,
            self.reference.modulo_length(),
        )
    }

    fn candidate(self, placements: Vec<RawAlignment>) -> Candidate {
        Candidate {
            orientation: self.orientation,
            mapping: self.mapping,
            placements,
        }
    }
}

fn select_proven_candidate(
    forward: UpperBoundPlacement,
    reverse: UpperBoundPlacement,
    forward_mapping: &[usize],
    reverse_mapping: &[usize],
) -> Result<Option<Candidate>> {
    use UpperBoundPlacement::{Ambiguous, Unattained, Unique, Unproven};

    match (forward, reverse) {
        (Unproven, _) | (_, Unproven) | (Unattained, Unattained) => Ok(None),
        (Unique(raw), Unattained) => Ok(Some(Candidate {
            orientation: Orientation::Forward,
            mapping: forward_mapping.to_vec(),
            placements: vec![raw],
        })),
        (Unattained, Unique(raw)) => Ok(Some(Candidate {
            orientation: Orientation::Reverse,
            mapping: reverse_mapping.to_vec(),
            placements: vec![raw],
        })),
        (Ambiguous, Unattained) | (Unattained, Ambiguous) => {
            Err(AlignmentError::AmbiguousPlacement.into())
        }
        (Unique(_) | Ambiguous, Unique(_) | Ambiguous) => {
            Err(AlignmentError::OrientationTie.into())
        }
    }
}

/// The calls two placements of one orientation place at different reference
/// positions (an inserted call at none), in query order.
fn disputed_calls(
    selected: &Candidate,
    left: &RawAlignment,
    right: &RawAlignment,
    modulo_length: Option<usize>,
) -> Vec<usize> {
    let positions = |alignment: &RawAlignment| {
        let mut positions = vec![None; selected.mapping.len()];
        for column in &alignment.columns {
            if let Some(query) = column.query_index
                && let Some(slot) = positions.get_mut(query)
            {
                *slot = column
                    .reference_index
                    .map(|index| modulo_length.map_or(index, |length| index % length));
            }
        }
        positions
    };
    positions(left)
        .into_iter()
        .zip(positions(right))
        .zip(&selected.mapping)
        .filter(|((left, right), _)| left != right)
        .map(|(_, &call)| call)
        .collect()
}

/// Checks the selected placement and builds the alignment, or reports the
/// calls two equally scoring placements dispute. On a vouched read, a
/// deletion beside an `N` without evidence moves 3' across it
/// ([`canonical::cross_unresolved`]).
fn finish_alignment(
    selected: &Candidate,
    reference: &Reference,
    config: &AlignmentConfig,
    calls: &[CallEvidence],
    vouched: bool,
) -> Result<Selection> {
    let modulo_length = match reference.topology {
        ReferenceTopology::Linear => None,
        ReferenceTopology::Circular => Some(reference.len()),
    };
    let raw = match selected.placements.as_slice() {
        [raw] => raw,
        [left, right] => {
            return Ok(Selection::Tied(disputed_calls(
                selected,
                left,
                right,
                modulo_length,
            )));
        }
        _ => return Err(AlignmentError::AmbiguousPlacement.into()),
    };
    if raw.metrics.callable_columns < config.minimum_callable_bases {
        return Err(AlignmentError::TooFewCallableColumns {
            found: raw.metrics.callable_columns,
            minimum: config.minimum_callable_bases,
        }
        .into());
    }
    if raw.metrics.callable_identity < config.minimum_identity {
        return Err(AlignmentError::LowIdentity {
            identity: raw.metrics.callable_identity,
            minimum: config.minimum_identity,
        }
        .into());
    }
    let mut raw = raw.clone();
    if vouched {
        let unevidenced = |query: usize| {
            selected
                .mapping
                .get(query)
                .and_then(|&call| calls.get(call))
                .is_some_and(|call| call.base == 'N' && call.profile.is_none())
        };
        canonical::cross_unresolved(&mut raw, unevidenced, modulo_length);
    }
    let masked = |query: usize| {
        selected
            .mapping
            .get(query)
            .and_then(|&call| calls.get(call))
            .is_some_and(|call| call.placement != Placement::Scored)
    };
    runs::reexpress(&mut raw, reference, &masked, config.max_stretch);
    let raw = &raw;
    let (segments, wraps_origin) = segments(raw, reference);
    let columns = raw
        .columns
        .iter()
        .map(|column| AlignmentColumn {
            query_base: column.query_base,
            reference_base: column.reference_base,
            original_call_index_0based: column
                .query_index
                .and_then(|index| selected.mapping.get(index).copied()),
            reference_index_0based: column
                .reference_index
                .map(|index| match reference.topology {
                    ReferenceTopology::Linear => index,
                    ReferenceTopology::Circular => index % reference.len(),
                }),
        })
        .collect::<Vec<_>>();
    let mask_of = |column: &AlignmentColumn| {
        column
            .original_call_index_0based
            .and_then(|index| calls.get(index))
            .map(|call| call.placement)
            .filter(|placement| *placement != Placement::Scored)
    };
    let masked_query_bases = columns
        .iter()
        .filter(|column| mask_of(column).is_some())
        .count();
    let hidden = columns
        .iter()
        .filter(|column| unresolved(mask_of(column)))
        .count();
    let mut metrics = raw.metrics.clone();
    metrics.unresolved_query_bases =
        metrics
            .unresolved_query_bases
            .checked_sub(hidden)
            .ok_or(AlignmentError::Inconsistent(
                "masked query bases exceed unresolved query bases",
            ))?;
    metrics.masked_query_bases = masked_query_bases;
    let callable_segments = callable_segments(&columns, |column| {
        mask_of(column).is_none() && is_canonical(column.query_base)
    });
    Ok(Selection::Unique(Alignment {
        orientation: selected.orientation,
        score: raw.score,
        reference_segments: segments,
        callable_segments,
        wraps_origin,
        metrics,
        columns,
    }))
}

/// Reference segments observed by callable calls — unmasked calls that read
/// a resolved base: a call column counts when its call is callable, a deletion
/// column when the call columns on both sides are; runs of consecutive
/// reference indexes form one segment. An unresolved call (`N`, or an IUPAC
/// code a modality presents as `N`) observes no base, masked or not.
fn callable_segments(
    columns: &[AlignmentColumn],
    callable: impl Fn(&AlignmentColumn) -> bool,
) -> Vec<ReferenceSegment> {
    let calls = columns
        .iter()
        .enumerate()
        .filter(|(_, column)| column.original_call_index_0based.is_some())
        .map(|(index, column)| (index, callable(column)))
        .collect::<Vec<_>>();
    let callable_call_around = |index: usize| {
        let before = calls.iter().rev().find(|(position, _)| *position < index);
        let after = calls.iter().find(|(position, _)| *position > index);
        before.is_some_and(|(_, callable)| *callable)
            && after.is_some_and(|(_, callable)| *callable)
    };
    let mut segments: Vec<ReferenceSegment> = Vec::new();
    for (index, column) in columns.iter().enumerate() {
        let Some(reference) = column.reference_index_0based else {
            continue;
        };
        let observed = if column.original_call_index_0based.is_some() {
            callable(column)
        } else {
            callable_call_around(index)
        };
        if !observed {
            continue;
        }
        match segments.last_mut() {
            Some(last) if last.end_0based_exclusive == reference => {
                last.end_0based_exclusive += 1;
            }
            _ => segments.push(ReferenceSegment {
                start_0based: reference,
                end_0based_exclusive: reference + 1,
            }),
        }
    }
    segments
}

/// Whether a call aligns as unresolved.
fn unresolved(placement: Option<Placement>) -> bool {
    matches!(placement, Some(Placement::Unresolved(_)))
}

fn reverse_profiles(profiles: &[Option<EvidenceProfile>]) -> Vec<Option<EvidenceProfile>> {
    profiles
        .iter()
        .rev()
        .map(|profile| profile.map(EvidenceProfile::complemented))
        .collect()
}

fn compare(left: &RawAlignment, right: &RawAlignment) -> Ordering {
    left.score.cmp(&right.score)
}

fn segments(alignment: &RawAlignment, reference: &Reference) -> (Vec<ReferenceSegment>, bool) {
    match reference.topology {
        ReferenceTopology::Linear => (
            vec![ReferenceSegment {
                start_0based: alignment.start_reference,
                end_0based_exclusive: alignment.end_reference,
            }],
            false,
        ),
        ReferenceTopology::Circular => {
            let length = reference.len();
            let start = alignment.start_reference % length;
            let span = alignment.end_reference - alignment.start_reference;
            let unwrapped_end = start + span;
            if unwrapped_end <= length {
                (
                    vec![ReferenceSegment {
                        start_0based: start,
                        end_0based_exclusive: unwrapped_end,
                    }],
                    false,
                )
            } else {
                (
                    vec![
                        ReferenceSegment {
                            start_0based: start,
                            end_0based_exclusive: length,
                        },
                        ReferenceSegment {
                            start_0based: 0,
                            end_0based_exclusive: unwrapped_end - length,
                        },
                    ],
                    true,
                )
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::model::alignment::AlignmentMetrics;
    use dna_kernel::model::reference::ReferenceTopology;
    use dna_kernel::read_evidence::{EvidenceReason, UnreliableSpan};

    use super::*;

    fn raw(score: i64, exact_matches: usize, mismatches: usize, gap_opens: usize) -> RawAlignment {
        RawAlignment {
            score,
            start_reference: 0,
            end_reference: 1,
            columns: Vec::new(),
            metrics: AlignmentMetrics {
                exact_matches,
                mismatches,
                gap_opens,
                callable_columns: exact_matches + mismatches,
                callable_identity: 0.0,
                unresolved_query_bases: 0,
                masked_query_bases: 0,
            },
        }
    }

    #[test]
    fn reverse_profiles_reverse_order_and_complement_channels() {
        let profiles = vec![
            Some(EvidenceProfile {
                weights: [1.0, 0.0, 0.0, 0.0],
            }),
            None,
            Some(EvidenceProfile {
                weights: [0.0, 1.0, 0.0, 0.0],
            }),
        ];
        let reversed = reverse_profiles(&profiles);
        assert_eq!(
            reversed[0].map(|profile| profile.weights),
            Some([0.0, 0.0, 1.0, 0.0])
        );
        assert!(reversed[1].is_none());
        assert_eq!(
            reversed[2].map(|profile| profile.weights),
            Some([0.0, 0.0, 0.0, 1.0])
        );
    }

    fn config() -> AlignmentConfig {
        AlignmentConfig {
            match_score: 3,
            mismatch_score: -5,
            ambiguous_score: 0,
            gap_open_score: -10,
            gap_extension_score: -4,
            minimum_callable_bases: 1,
            minimum_identity: 0.8,
            max_stretch: 64,
        }
    }

    fn reference() -> Reference {
        Reference {
            name: "ref".into(),
            sequence: "GCCAAAAGTT".into(),
            topology: ReferenceTopology::Linear,
            sequence_sha256: String::new(),
        }
    }

    // SRS-ALN-012
    #[test]
    fn forward_and_reverse_reads_share_canonical_deletion_coordinate() -> Result<()> {
        let forward = align_best(&ReadEvidence::clean("GCCAAAGTT"), &reference(), &config())?;
        let reverse = align_best(&ReadEvidence::clean("AACTTTGGC"), &reference(), &config())?;

        assert_eq!(forward.orientation, Orientation::Forward);
        assert_eq!(reverse.orientation, Orientation::Reverse);

        let forward_deleted = forward
            .columns
            .iter()
            .find(|column| column.query_base == '-')
            .and_then(|column| column.reference_index_0based);
        let reverse_deleted = reverse
            .columns
            .iter()
            .find(|column| column.query_base == '-')
            .and_then(|column| column.reference_index_0based);
        assert_eq!(forward_deleted, Some(6));
        assert_eq!(reverse_deleted, forward_deleted);
        Ok(())
    }

    // SRS-ALN-006, SRS-ALN-013
    #[test]
    fn upper_bound_orientation_tie_preserves_existing_error() {
        let reference = Reference {
            name: "ref".into(),
            sequence: "ACGT".into(),
            topology: ReferenceTopology::Linear,
            sequence_sha256: String::new(),
        };
        let error = align_best(&ReadEvidence::clean("ACGT"), &reference, &config())
            .err()
            .map(|error| error.to_string());
        assert_eq!(
            error.as_deref(),
            Some("alignment failed: forward and reverse evidence-profile scores are tied")
        );
    }

    // SRS-ALN-006, SRS-ALN-013
    #[test]
    fn repeated_upper_bound_placement_preserves_existing_error() {
        let reference = Reference {
            name: "ref".into(),
            sequence: "AAAAA".into(),
            topology: ReferenceTopology::Linear,
            sequence_sha256: String::new(),
        };
        let error = align_best(&ReadEvidence::clean("AAA"), &reference, &config())
            .err()
            .map(|error| error.to_string());
        assert_eq!(
            error.as_deref(),
            Some("alignment failed: selected orientation has multiple equally scoring placements")
        );
    }

    /// Deterministic xorshift bases.
    fn random_bases(mut state: u64, length: usize) -> String {
        (0..length)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                ['A', 'C', 'G', 'T'][usize::from(u8::try_from(state % 4).unwrap_or(0))]
            })
            .collect()
    }

    fn working(reference: &Reference) -> (String, Option<usize>) {
        match reference.topology {
            ReferenceTopology::Linear => (reference.sequence.clone(), None),
            ReferenceTopology::Circular => (
                format!("{}{}", reference.sequence, reference.sequence),
                Some(reference.len()),
            ),
        }
    }

    /// Both orientations' unmasked queries.
    fn queries(evidence: &ReadEvidence) -> [(String, Vec<Option<EvidenceProfile>>); 2] {
        let forward = evidence
            .calls()
            .iter()
            .map(|call| call.base)
            .collect::<String>();
        let profiles = evidence
            .calls()
            .iter()
            .map(|call| call.profile)
            .collect::<Vec<_>>();
        let reverse = (reverse_complement(&forward), reverse_profiles(&profiles));
        [(forward, profiles), reverse]
    }

    /// Seeded pruning of both orientations of unmasked evidence.
    fn pruned_orientations(
        evidence: &ReadEvidence,
        reference: &Reference,
        config: &AlignmentConfig,
    ) -> Result<Vec<Pruned>> {
        let (sequence, modulo_length) = working(reference);
        let target = WorkingReference::new(&sequence, modulo_length);
        queries(evidence)
            .iter()
            .map(|(query, profiles)| {
                Problem::new(query, profiles, &target, config)
                    .map_or(Ok(Pruned::Unproven { lower_bound: None }), |problem| {
                        exact::align_pruned(&problem)
                    })
            })
            .collect()
    }

    /// The selection full-reference Gotoh in both orientations makes for
    /// unmasked evidence.
    fn gotoh_selection(
        evidence: &ReadEvidence,
        reference: &Reference,
        config: &AlignmentConfig,
    ) -> Result<Alignment> {
        let (sequence, modulo_length) = working(reference);
        let [
            (forward_query, forward_profiles),
            (reverse_query, reverse_profiles),
        ] = queries(evidence);
        let forward = gotoh::align(
            &forward_query,
            &forward_profiles,
            &sequence,
            config,
            modulo_length,
        )?;
        let reverse = gotoh::align(
            &reverse_query,
            &reverse_profiles,
            &sequence,
            config,
            modulo_length,
        )?;
        let length = evidence.calls().len();
        let selected = match compare(&forward[0], &reverse[0]) {
            Ordering::Greater => Candidate {
                orientation: Orientation::Forward,
                mapping: (0..length).collect(),
                placements: forward,
            },
            Ordering::Less => Candidate {
                orientation: Orientation::Reverse,
                mapping: (0..length).rev().collect(),
                placements: reverse,
            },
            Ordering::Equal => return Err(AlignmentError::OrientationTie.into()),
        };
        match finish_alignment(
            &selected,
            reference,
            config,
            evidence.calls(),
            evidence.ends() == ReadEnds::Vouched,
        )? {
            Selection::Unique(alignment) => Ok(alignment),
            Selection::Tied(_) => Err(AlignmentError::AmbiguousPlacement.into()),
        }
    }

    fn assert_matches_gotoh_selection(
        evidence: &ReadEvidence,
        reference: &Reference,
        config: &AlignmentConfig,
    ) -> Result<Alignment> {
        let expected = gotoh_selection(evidence, reference, config)?;
        let actual = align_best(evidence, reference, config)?;
        assert_eq!(format!("{actual:?}"), format!("{expected:?}"));
        Ok(actual)
    }

    // SRS-ALN-014
    /// Neither orientation is proven when the read's probe windows cover most
    /// of the reference: three near copies of the read's source. Forward,
    /// first on equal lower bounds, runs full Gotoh, and its optimum proves
    /// that reverse cannot reach it.
    #[test]
    fn both_unproven_orientations_reject_the_other_by_threshold_proof() -> Result<()> {
        let source = random_bases(0x0dd5_eed5, 300);
        let copy = |substitutions: &[usize]| {
            let mut copy = source.clone().into_bytes();
            for &position in substitutions {
                copy[position] = if copy[position] == b'A' { b'C' } else { b'A' };
            }
            String::from_utf8(copy).unwrap_or_default()
        };
        let sequence = format!(
            "{}{}{}{}{}{}",
            random_bases(1, 100),
            source,
            random_bases(2, 50),
            copy(&[60, 200]),
            random_bases(3, 50),
            copy(&[70, 220]),
        );
        let reference = Reference {
            name: "ref".into(),
            sequence,
            topology: ReferenceTopology::Linear,
            sequence_sha256: String::new(),
        };
        let evidence = ReadEvidence::clean(&copy(&[40]));

        let pruned = pruned_orientations(&evidence, &reference, &config())?;
        assert!(matches!(
            pruned.as_slice(),
            [
                Pruned::Unproven { lower_bound: None },
                Pruned::Unproven { lower_bound: None }
            ]
        ));
        let forward = gotoh_selection(&evidence, &reference, &config())?;
        let (sequence, modulo_length) = working(&reference);
        let target = WorkingReference::new(&sequence, modulo_length);
        let [_, (reverse_query, reverse_profiles)] = queries(&evidence);
        let reverse = Problem::new(&reverse_query, &reverse_profiles, &target, &config())
            .map(|problem| exact::align_at_or_above(&problem, forward.score))
            .transpose()?;
        assert!(matches!(reverse, Some(Some(placements)) if placements.is_empty()));

        let selected = assert_matches_gotoh_selection(&evidence, &reference, &config())?;
        assert_eq!(selected.orientation, Orientation::Forward);
        assert_eq!(selected.reference_segments[0].start_0based, 100);
        Ok(())
    }

    // SRS-ALN-014
    /// Neither orientation is proven when most read loci are weak calls that
    /// the reference contradicts: too many seeds would be needed. The reverse
    /// read's seeds give a lower bound, so reverse runs full Gotoh first, and
    /// forward, unable to prove the threshold, runs full Gotoh too.
    #[test]
    fn both_unproven_orientations_start_from_the_higher_lower_bound() -> Result<()> {
        let sequence = random_bases(0xfeed_beef, 600);
        let source = &sequence[150..450];
        let clean = |index: usize| (20..40).contains(&(index % 75));
        let mut evidence = ReadEvidence::clean(&reverse_complement(source));
        for (index, base) in source.bytes().enumerate() {
            if clean(index) {
                continue;
            }
            let channel = "ACGT".find(char::from(base)).unwrap_or(0);
            let wrong = (channel + 1) % 4;
            let mut weights = [0.0; 4];
            weights[wrong] = 0.85;
            weights[channel] = 0.15;
            let call_index = source.len() - 1 - index;
            let mut call = evidence.calls()[call_index];
            call.base = reverse_complement(&"ACGT"[wrong..=wrong])
                .chars()
                .next()
                .unwrap_or('N');
            call.profile = Some(EvidenceProfile { weights }.complemented());
            evidence = evidence.with_call(call_index, call);
        }
        let reference = Reference {
            name: "ref".into(),
            sequence,
            topology: ReferenceTopology::Linear,
            sequence_sha256: String::new(),
        };
        let mut config = config();
        config.minimum_identity = 0.0;

        let pruned = pruned_orientations(&evidence, &reference, &config)?;
        assert!(matches!(
            pruned.as_slice(),
            [
                Pruned::Unproven { lower_bound: None },
                Pruned::Unproven {
                    lower_bound: Some(_)
                }
            ]
        ));

        let selected = assert_matches_gotoh_selection(&evidence, &reference, &config)?;
        assert_eq!(selected.orientation, Orientation::Reverse);
        Ok(())
    }

    // SRS-ALN-015
    #[test]
    fn callable_segments_skip_uncallable_calls_and_deletions_beside_them() {
        let column = |call: Option<usize>, reference: Option<usize>| AlignmentColumn {
            query_base: if call.is_some() { 'A' } else { '-' },
            reference_base: if reference.is_some() { 'A' } else { '-' },
            original_call_index_0based: call,
            reference_index_0based: reference,
        };
        let columns = [
            column(Some(0), Some(10)),
            column(None, Some(11)),
            column(Some(1), Some(12)),
            column(Some(2), Some(13)),
            column(None, Some(14)),
            column(Some(3), Some(15)),
            column(Some(4), None),
            column(Some(5), Some(16)),
        ];
        let callable = |column: &AlignmentColumn| column.original_call_index_0based != Some(3);
        assert_eq!(
            callable_segments(&columns, callable),
            [
                ReferenceSegment {
                    start_0based: 10,
                    end_0based_exclusive: 14,
                },
                ReferenceSegment {
                    start_0based: 16,
                    end_0based_exclusive: 17,
                },
            ]
        );
    }

    // SRS-ALN-015
    #[test]
    fn unresolved_unmasked_calls_are_not_callable() -> Result<()> {
        let reference = Reference {
            name: "ref".into(),
            sequence: "TTGACGTACGGATCCATGCAGTCA".into(),
            topology: ReferenceTopology::Linear,
            sequence_sha256: String::new(),
        };
        // An `N` run and an ambiguity code: unresolved bases, the second with
        // a shared-weight profile, neither masked.
        let mut evidence = ReadEvidence::clean("TTGACGNNNGGATCCATGCAGTCA");
        let mut call = evidence.calls()[14];
        call.base = 'N';
        call.profile = Some(EvidenceProfile {
            weights: [0.0, 0.5, 0.0, 0.5],
        });
        evidence = evidence.with_call(14, call);
        let alignment = align_best(&evidence, &reference, &config())?;
        let segment = |start_0based, end_0based_exclusive| ReferenceSegment {
            start_0based,
            end_0based_exclusive,
        };
        assert_eq!(alignment.reference_segments, [segment(0, 24)]);
        assert_eq!(
            alignment.callable_segments,
            [segment(0, 6), segment(9, 14), segment(15, 24)]
        );
        assert_eq!(alignment.metrics.unresolved_query_bases, 4);
        assert_eq!(alignment.metrics.masked_query_bases, 0);
        Ok(())
    }

    // SRS-ALN-011
    #[test]
    fn metrics_derive_from_the_retained_primary_sequence() -> Result<()> {
        let reference = Reference {
            name: "ref".into(),
            sequence: "TTGACGTACGGATCCATGCAGTCA".into(),
            topology: ReferenceTopology::Linear,
            sequence_sha256: String::new(),
        };
        // One substitution (T>A at 6) and one unresolved base (N at 14).
        let alignment = align_best(
            &ReadEvidence::clean("TTGACGAACGGATCNATGCAGTCA"),
            &reference,
            &config(),
        )?;
        let metrics = alignment.metrics();
        assert_eq!(metrics.callable_columns, 23);
        assert_eq!(metrics.exact_matches(), 22);
        assert_eq!(metrics.mismatches(), 1);
        assert_eq!(metrics.unresolved_query_bases, 1);
        assert!((metrics.callable_identity - 22.0 / 23.0).abs() < 1e-12);
        // The same counts follow from the selected traceback's columns.
        let canonical = |base: char| "ACGT".contains(base);
        let callable = alignment
            .columns
            .iter()
            .filter(|column| canonical(column.query_base) && canonical(column.reference_base))
            .collect::<Vec<_>>();
        assert_eq!(callable.len(), metrics.callable_columns);
        assert_eq!(
            callable
                .iter()
                .filter(|column| column.query_base == column.reference_base)
                .count(),
            metrics.exact_matches()
        );
        Ok(())
    }

    // SRS-ALN-006
    #[test]
    fn orientation_comparison_uses_profile_score_only() {
        let left = raw(100, 1, 9, 5);
        let right = raw(100, 10, 0, 0);
        assert_eq!(compare(&left, &right), Ordering::Equal);
        assert_eq!(compare(&raw(101, 0, 10, 10), &right), Ordering::Greater);
    }

    // SRS-ALN-016
    #[test]
    fn re_expresses_a_masked_repeat_stretch_one_length_edit_per_run() -> Result<()> {
        let reference = Reference {
            name: "ref".into(),
            sequence: "TTGACCCCCCCTCCCCCGTAGG".into(),
            topology: ReferenceTopology::Linear,
            sequence_sha256: String::new(),
        };
        let clean = ReadEvidence::clean("TTGACCCCCCCCTCCCCCCGTAGG");
        let mut call = clean.calls()[12];
        call = call.masked(Placement::Anchoring(EvidenceReason::new("dephased_signal")));
        let masked = clean.clone().with_call(12, call);
        let inserted_after = |alignment: &Alignment| {
            alignment
                .columns
                .windows(2)
                .filter(|pair| pair[1].reference_base == '-')
                .filter_map(|pair| pair[0].reference_index_0based)
                .collect::<Vec<_>>()
        };

        // Unmasked evidence keeps the optimal path: one merged gap and a
        // substitution.
        let optimal = align_best(&clean, &reference, &config())?;
        assert_eq!(optimal.metrics.mismatches, 1);
        // With the T masked, each C run gains one base at its 3' end.
        let realigned = align_best(&masked, &reference, &config())?;
        assert_eq!(realigned.metrics.mismatches, 0);
        assert_eq!(inserted_after(&realigned), [10, 16]);
        assert_eq!(realigned.score, optimal.score);
        // A reverse read lands on the same reference coordinates.
        let reverse = ReadEvidence::clean(&reverse_complement("TTGACCCCCCCCTCCCCCCGTAGG"));
        let mut call = reverse.calls()[11];
        call = call.masked(Placement::Anchoring(EvidenceReason::new("dephased_signal")));
        let reverse = align_best(&reverse.clone().with_call(11, call), &reference, &config())?;
        assert_eq!(reverse.orientation, Orientation::Reverse);
        assert_eq!(inserted_after(&reverse), [10, 16]);
        Ok(())
    }

    /// A read whose last three calls tie: two inserted A's beside a matched C
    /// (the read ends at reference 25) or three substitutions (it ends at
    /// reference 27). Both score 51 and agree on every earlier call.
    fn tied_end() -> (Reference, ReadEvidence) {
        let reference = Reference {
            name: "ref".into(),
            sequence: "TTGATTCAGCTAACGGTCATCAAACCATTTT".into(),
            topology: ReferenceTopology::Linear,
            sequence_sha256: String::new(),
        };
        (reference, ReadEvidence::clean("GATTCAGCTAACGGTCATCAAAAAC"))
    }

    fn unreliable(evidence: ReadEvidence, calls: Range<usize>) -> Result<ReadEvidence> {
        evidence.with_unreliable(vec![UnreliableSpan {
            calls,
            reason: EvidenceReason::new("poor_start_signal"),
        }])
    }

    fn placed(alignment: &Alignment) -> Vec<usize> {
        alignment
            .columns
            .iter()
            .filter_map(|column| column.original_call_index_0based)
            .collect()
    }

    fn ambiguous(result: Result<Alignment>) -> bool {
        result.err().map(|error| error.to_string()).as_deref()
            == Some(
                "alignment failed: selected orientation has multiple equally scoring placements",
            )
    }

    // SRS-ALN-006
    #[test]
    fn a_tie_outside_unreliable_calls_stays_rejected() -> Result<()> {
        let (reference, evidence) = tied_end();
        assert!(ambiguous(align_best(&evidence, &reference, &config())));
        // A span that stops short of the first disputed call leaves a
        // reliable call in dispute.
        let short = unreliable(evidence, 23..25)?;
        assert!(ambiguous(align_best(&short, &reference, &config())));
        Ok(())
    }

    // SRS-ALN-006
    #[test]
    fn a_tie_confined_to_an_unreliable_span_leaves_its_calls_unplaced() -> Result<()> {
        let (reference, evidence) = tied_end();
        let alignment = align_best(&unreliable(evidence, 20..25)?, &reference, &config())?;
        assert_eq!(alignment.orientation, Orientation::Forward);
        assert_eq!(placed(&alignment), (0..22).collect::<Vec<_>>());
        assert_eq!(
            alignment.reference_segments,
            [ReferenceSegment {
                start_0based: 2,
                end_0based_exclusive: 24,
            }]
        );
        Ok(())
    }

    // SRS-ALN-006
    #[test]
    fn a_tie_confined_to_masked_calls_leaves_them_unplaced() -> Result<()> {
        let (reference, mut evidence) = tied_end();
        for index in 22..25 {
            let call = evidence.calls()[index]
                .masked(Placement::Anchoring(EvidenceReason::new("dephased_signal")));
            evidence = evidence.with_call(index, call);
        }
        let alignment = align_best(&evidence, &reference, &config())?;
        assert_eq!(placed(&alignment), (0..22).collect::<Vec<_>>());
        Ok(())
    }

    // SRS-ALN-006
    #[test]
    fn a_reverse_read_leaves_its_unreliable_start_unplaced() -> Result<()> {
        let (reference, evidence) = tied_end();
        let read = evidence
            .calls()
            .iter()
            .map(|call| call.base)
            .collect::<String>();
        let reverse = ReadEvidence::clean(&reverse_complement(&read));
        assert!(ambiguous(align_best(&reverse, &reference, &config())));
        let alignment = align_best(&unreliable(reverse, 0..5)?, &reference, &config())?;
        assert_eq!(alignment.orientation, Orientation::Reverse);
        assert_eq!(placed(&alignment), (3..25).rev().collect::<Vec<_>>());
        Ok(())
    }
}
