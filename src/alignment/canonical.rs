//! Canonical right alignment for repeat-equivalent optimal gap placements.

use crate::alignment::scoring::{scaled, substitution};
use crate::alignment::traceback::{self, RawAlignment, RawColumn};
use crate::config::AlignmentConfig;
use crate::error::{Error, Result};
use crate::model::locus_evidence::EvidenceProfile;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum GapKind {
    Insertion,
    Deletion,
}

#[derive(Debug, Clone, Copy)]
struct GapRun {
    kind: GapKind,
    start: usize,
    end: usize,
}

pub(crate) fn right_align(
    alignment: &mut RawAlignment,
    profiles: &[Option<EvidenceProfile>],
    config: &AlignmentConfig,
    modulo_length: Option<usize>,
) -> Result<()> {
    let expected_score = alignment.score;
    let observed_score = score(&alignment.columns, profiles, config)?;
    if observed_score != expected_score {
        return Err(Error::Alignment(format!(
            "traceback score {observed_score} disagrees with DP score {expected_score}"
        )));
    }

    loop {
        let runs = gap_runs(&alignment.columns);
        let mut shifted = false;
        for run in runs.into_iter().rev() {
            let Some(candidate) = shifted_right(&alignment.columns, run, modulo_length) else {
                continue;
            };
            if score(&candidate, profiles, config)? != expected_score {
                continue;
            }
            alignment.columns = candidate;
            alignment.metrics = traceback::metrics(&alignment.columns);
            shifted = true;
            break;
        }
        if !shifted {
            return Ok(());
        }
    }
}

fn gap_runs(columns: &[RawColumn]) -> Vec<GapRun> {
    let mut runs = Vec::new();
    let mut index = 0;
    while index < columns.len() {
        let Some(kind) = gap_kind(&columns[index]) else {
            index += 1;
            continue;
        };
        let start = index;
        index += 1;
        while index < columns.len() && gap_kind(&columns[index]) == Some(kind) {
            index += 1;
        }
        runs.push(GapRun {
            kind,
            start,
            end: index,
        });
    }
    runs
}

fn shifted_right(
    columns: &[RawColumn],
    run: GapRun,
    modulo_length: Option<usize>,
) -> Option<Vec<RawColumn>> {
    let following = columns.get(run.end)?;
    if gap_kind(following).is_some() {
        return None;
    }

    match run.kind {
        GapKind::Deletion => shift_deletion(columns, run, following, modulo_length),
        GapKind::Insertion => shift_insertion(columns, run, following, modulo_length),
    }
}

fn shift_deletion(
    columns: &[RawColumn],
    run: GapRun,
    following: &RawColumn,
    modulo_length: Option<usize>,
) -> Option<Vec<RawColumn>> {
    let first = columns.get(run.start)?;
    if first.reference_base != following.reference_base {
        return None;
    }
    if crosses_circular_seam(
        first.reference_index?,
        following.reference_index?,
        modulo_length,
    ) {
        return None;
    }

    let mut shifted = columns.to_vec();
    shifted[run.start].query_base = following.query_base;
    shifted[run.start].query_index = following.query_index;
    shifted[run.end].query_base = '-';
    shifted[run.end].query_index = None;
    Some(shifted)
}

fn shift_insertion(
    columns: &[RawColumn],
    run: GapRun,
    following: &RawColumn,
    modulo_length: Option<usize>,
) -> Option<Vec<RawColumn>> {
    let first = columns.get(run.start)?;
    if first.query_base != following.reference_base {
        return None;
    }
    let following_reference = following.reference_index?;
    if modulo_length.is_some_and(|length| following_reference % length == 0) {
        return None;
    }

    let mut shifted = columns.to_vec();
    shifted[run.start].reference_base = following.reference_base;
    shifted[run.start].reference_index = following.reference_index;
    shifted[run.end].reference_base = '-';
    shifted[run.end].reference_index = None;
    Some(shifted)
}

fn crosses_circular_seam(
    first_reference: usize,
    following_reference: usize,
    modulo_length: Option<usize>,
) -> bool {
    modulo_length.is_some_and(|length| first_reference / length != following_reference / length)
}

fn score(
    columns: &[RawColumn],
    profiles: &[Option<EvidenceProfile>],
    config: &AlignmentConfig,
) -> Result<i64> {
    let mut total = 0_i64;
    let mut previous_gap = None;

    for column in columns {
        let gap = gap_kind(column);
        if let Some(kind) = gap {
            if previous_gap != Some(kind) {
                total = total
                    .checked_add(scaled(config.gap_open_score))
                    .ok_or_else(|| Error::Alignment("canonical alignment score overflow".into()))?;
            }
            total = total
                .checked_add(scaled(config.gap_extension_score))
                .ok_or_else(|| Error::Alignment("canonical alignment score overflow".into()))?;
            previous_gap = Some(kind);
            continue;
        }

        previous_gap = None;
        let query_index = column
            .query_index
            .ok_or_else(|| Error::Alignment("aligned canonical column lacks query index".into()))?;
        let profile = profiles.get(query_index).copied().ok_or_else(|| {
            Error::Alignment(format!(
                "canonical alignment query index {query_index} is out of profile bounds"
            ))
        })?;
        total = total
            .checked_add(substitution(profile, column.reference_base as u8, config))
            .ok_or_else(|| Error::Alignment("canonical alignment score overflow".into()))?;
    }

    Ok(total)
}

const fn gap_kind(column: &RawColumn) -> Option<GapKind> {
    if column.query_base == '-' {
        Some(GapKind::Deletion)
    } else if column.reference_base == '-' {
        Some(GapKind::Insertion)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn config() -> AlignmentConfig {
        AlignmentConfig {
            match_score: 3,
            mismatch_score: -5,
            ambiguous_score: 0,
            gap_open_score: -10,
            gap_extension_score: -4,
            minimum_callable_bases: 1,
            minimum_identity: 0.8,
        }
    }

    fn one_hot(sequence: &str) -> Vec<Option<EvidenceProfile>> {
        sequence
            .bytes()
            .map(|base| {
                let weights = match base {
                    b'A' => [1.0, 0.0, 0.0, 0.0],
                    b'C' => [0.0, 1.0, 0.0, 0.0],
                    b'G' => [0.0, 0.0, 1.0, 0.0],
                    b'T' => [0.0, 0.0, 0.0, 1.0],
                    _ => return None,
                };
                Some(EvidenceProfile { weights })
            })
            .collect()
    }

    fn raw(columns: Vec<RawColumn>, score: i64) -> RawAlignment {
        let metrics = traceback::metrics(&columns);
        RawAlignment {
            score,
            start_reference: 0,
            end_reference: columns
                .iter()
                .filter_map(|column| column.reference_index)
                .max()
                .map_or(0, |index| index + 1),
            columns,
            metrics,
        }
    }

    fn column(
        query_base: char,
        reference_base: char,
        query_index: Option<usize>,
        reference_index: Option<usize>,
    ) -> RawColumn {
        RawColumn {
            query_base,
            reference_base,
            query_index,
            reference_index,
        }
    }

    #[test]
    fn right_shifts_homopolymer_deletion() -> Result<()> {
        let profiles = one_hot("CAAAG");
        let columns = vec![
            column('C', 'C', Some(0), Some(0)),
            column('-', 'A', None, Some(1)),
            column('A', 'A', Some(1), Some(2)),
            column('A', 'A', Some(2), Some(3)),
            column('A', 'A', Some(3), Some(4)),
            column('G', 'G', Some(4), Some(5)),
        ];
        let expected = score(&columns, &profiles, &config())?;
        let mut alignment = raw(columns, expected);
        right_align(&mut alignment, &profiles, &config(), None)?;

        assert_eq!(
            alignment
                .columns
                .iter()
                .position(|column| column.query_base == '-'),
            Some(4)
        );
        Ok(())
    }

    #[test]
    fn right_shifts_homopolymer_insertion_when_profile_score_is_equal() -> Result<()> {
        let profiles = one_hot("CAAAAAG");
        let columns = vec![
            column('C', 'C', Some(0), Some(0)),
            column('A', '-', Some(1), None),
            column('A', 'A', Some(2), Some(1)),
            column('A', 'A', Some(3), Some(2)),
            column('A', 'A', Some(4), Some(3)),
            column('A', 'A', Some(5), Some(4)),
            column('G', 'G', Some(6), Some(5)),
        ];
        let expected = score(&columns, &profiles, &config())?;
        let mut alignment = raw(columns, expected);
        right_align(&mut alignment, &profiles, &config(), None)?;

        assert_eq!(
            alignment
                .columns
                .iter()
                .position(|column| column.reference_base == '-'),
            Some(5)
        );
        Ok(())
    }

    #[test]
    fn profile_score_can_prevent_equivalent_primary_insertion_shift() -> Result<()> {
        let mut profiles = one_hot("CAAAG");
        profiles[1] = Some(EvidenceProfile {
            weights: [0.0, 1.0, 0.0, 0.0],
        });
        let columns = vec![
            column('C', 'C', Some(0), Some(0)),
            column('A', '-', Some(1), None),
            column('A', 'A', Some(2), Some(1)),
            column('A', 'A', Some(3), Some(2)),
            column('G', 'G', Some(4), Some(3)),
        ];
        let expected = score(&columns, &profiles, &config())?;
        let mut alignment = raw(columns, expected);
        right_align(&mut alignment, &profiles, &config(), None)?;

        assert_eq!(
            alignment
                .columns
                .iter()
                .position(|column| column.reference_base == '-'),
            Some(1)
        );
        Ok(())
    }

    #[test]
    fn right_shifts_tandem_repeat_deletion_by_whole_motif() -> Result<()> {
        let profiles = one_hot("CATG");
        let columns = vec![
            column('C', 'C', Some(0), Some(0)),
            column('-', 'A', None, Some(1)),
            column('-', 'T', None, Some(2)),
            column('A', 'A', Some(1), Some(3)),
            column('T', 'T', Some(2), Some(4)),
            column('G', 'G', Some(3), Some(5)),
        ];
        let expected = score(&columns, &profiles, &config())?;
        let mut alignment = raw(columns, expected);
        right_align(&mut alignment, &profiles, &config(), None)?;

        let deleted = alignment
            .columns
            .iter()
            .filter(|column| column.query_base == '-')
            .filter_map(|column| column.reference_index)
            .collect::<Vec<_>>();
        assert_eq!(deleted, vec![3, 4]);
        Ok(())
    }

    #[test]
    fn right_shifts_tandem_repeat_insertion_by_whole_motif() -> Result<()> {
        let profiles = one_hot("CATATATG");
        let columns = vec![
            column('C', 'C', Some(0), Some(0)),
            column('A', '-', Some(1), None),
            column('T', '-', Some(2), None),
            column('A', 'A', Some(3), Some(1)),
            column('T', 'T', Some(4), Some(2)),
            column('A', 'A', Some(5), Some(3)),
            column('T', 'T', Some(6), Some(4)),
            column('G', 'G', Some(7), Some(5)),
        ];
        let expected = score(&columns, &profiles, &config())?;
        let mut alignment = raw(columns, expected);
        right_align(&mut alignment, &profiles, &config(), None)?;

        let insertion = alignment
            .columns
            .iter()
            .enumerate()
            .filter(|(_, column)| column.reference_base == '-')
            .map(|(index, column)| (index, column.query_base))
            .collect::<Vec<_>>();
        assert_eq!(insertion, vec![(5, 'A'), (6, 'T')]);
        assert_eq!(alignment.columns[4].reference_index, Some(4));
        assert_eq!(alignment.columns[7].reference_index, Some(5));
        Ok(())
    }

    #[test]
    fn canonicalizes_multiple_gap_runs_from_right_to_left() -> Result<()> {
        let profiles = one_hot("CAAGTTG");
        let columns = vec![
            column('C', 'C', Some(0), Some(0)),
            column('-', 'A', None, Some(1)),
            column('A', 'A', Some(1), Some(2)),
            column('A', 'A', Some(2), Some(3)),
            column('G', 'G', Some(3), Some(4)),
            column('-', 'T', None, Some(5)),
            column('T', 'T', Some(4), Some(6)),
            column('T', 'T', Some(5), Some(7)),
            column('G', 'G', Some(6), Some(8)),
        ];
        let expected = score(&columns, &profiles, &config())?;
        let mut alignment = raw(columns, expected);
        right_align(&mut alignment, &profiles, &config(), None)?;

        let deleted = alignment
            .columns
            .iter()
            .filter(|column| column.query_base == '-')
            .filter_map(|column| column.reference_index)
            .collect::<Vec<_>>();
        assert_eq!(deleted, vec![3, 7]);
        Ok(())
    }

    #[test]
    fn does_not_shift_non_repeat_deletion() -> Result<()> {
        let profiles = one_hot("CAG");
        let columns = vec![
            column('C', 'C', Some(0), Some(0)),
            column('-', 'A', None, Some(1)),
            column('A', 'T', Some(1), Some(2)),
            column('G', 'G', Some(2), Some(3)),
        ];
        let expected = score(&columns, &profiles, &config())?;
        let mut alignment = raw(columns, expected);
        right_align(&mut alignment, &profiles, &config(), None)?;

        assert_eq!(alignment.columns[1].query_base, '-');
        Ok(())
    }

    #[test]
    fn circular_canonicalization_stops_at_origin_seam() -> Result<()> {
        let profiles = one_hot("AAA");
        let columns = vec![
            column('A', 'A', Some(0), Some(2)),
            column('-', 'A', None, Some(3)),
            column('A', 'A', Some(1), Some(4)),
            column('A', 'A', Some(2), Some(5)),
        ];
        let expected = score(&columns, &profiles, &config())?;
        let mut alignment = raw(columns, expected);
        right_align(&mut alignment, &profiles, &config(), Some(4))?;

        assert_eq!(alignment.columns[1].query_base, '-');
        assert_eq!(alignment.columns[1].reference_index, Some(3));
        Ok(())
    }

    #[test]
    fn canonicalization_preserves_metrics_for_equivalent_homopolymer_shift() -> Result<()> {
        let profiles = one_hot("CAAAG");
        let columns = vec![
            column('C', 'C', Some(0), Some(0)),
            column('-', 'A', None, Some(1)),
            column('A', 'A', Some(1), Some(2)),
            column('A', 'A', Some(2), Some(3)),
            column('A', 'A', Some(3), Some(4)),
            column('G', 'G', Some(4), Some(5)),
        ];
        let expected = score(&columns, &profiles, &config())?;
        let mut alignment = raw(columns, expected);
        let before = alignment.metrics.clone();
        right_align(&mut alignment, &profiles, &config(), None)?;

        assert_eq!(alignment.metrics.exact_matches, before.exact_matches);
        assert_eq!(alignment.metrics.mismatches, before.mismatches);
        assert_eq!(alignment.metrics.gap_opens, before.gap_opens);
        assert_eq!(alignment.metrics.callable_columns, before.callable_columns);
        assert_eq!(
            alignment.metrics.callable_identity,
            before.callable_identity
        );
        assert_eq!(
            alignment.metrics.unresolved_query_bases,
            before.unresolved_query_bases
        );
        assert_eq!(alignment.metrics.gap_opens, 1);
        Ok(())
    }
}
