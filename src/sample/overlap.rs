//! Pairwise overlap admission discovered from independently placed reads.

use std::collections::BTreeMap;

use crate::config::SampleReconciliationConfig;
use crate::error::{Error, Result};
use crate::model::read_observation::ReadObservation;
use crate::model::sample_evidence::{OverlapExclusionReason, ReadOverlapEvidence};

#[derive(Debug, Clone, Copy)]
struct CoordinateObservation {
    reference_base: char,
    query_base: Option<char>,
}

/// Builds the sparse pairwise overlap graph used by later sample reconciliation.
pub(super) fn assess(
    reads: &[&ReadObservation],
    config: &SampleReconciliationConfig,
) -> Result<Vec<ReadOverlapEvidence>> {
    let coordinates = reads
        .iter()
        .map(|read| coordinates(read))
        .collect::<Result<Vec<_>>>()?;
    let mut overlaps = Vec::new();

    for left_read_index in 0..reads.len() {
        for right_read_index in left_read_index + 1..reads.len() {
            if let Some(overlap) = assess_pair(
                left_read_index,
                right_read_index,
                &coordinates[left_read_index],
                &coordinates[right_read_index],
                config,
            )? {
                overlaps.push(overlap);
            }
        }
    }

    Ok(overlaps)
}

fn coordinates(read: &ReadObservation) -> Result<BTreeMap<usize, CoordinateObservation>> {
    let mut coordinates = BTreeMap::new();
    for column in &read.alignment.columns {
        let Some(reference_index_0based) = column.reference_index_0based else {
            continue;
        };
        if column.reference_base == '-' {
            return Err(Error::Sample(
                "reference-coordinate overlap cannot contain an insertion reference base".into(),
            ));
        }
        let query_base = is_canonical(column.query_base).then_some(column.query_base);
        if coordinates
            .insert(
                reference_index_0based,
                CoordinateObservation {
                    reference_base: column.reference_base,
                    query_base,
                },
            )
            .is_some()
        {
            return Err(Error::Sample(format!(
                "read {} contains duplicate reference coordinate {}",
                read.input_sha256,
                reference_index_0based + 1
            )));
        }
    }
    Ok(coordinates)
}

fn assess_pair(
    left_read_index: usize,
    right_read_index: usize,
    left: &BTreeMap<usize, CoordinateObservation>,
    right: &BTreeMap<usize, CoordinateObservation>,
    config: &SampleReconciliationConfig,
) -> Result<Option<ReadOverlapEvidence>> {
    let mut shared_positions = 0usize;
    let mut comparable_bases = 0usize;
    let mut agreements = 0usize;
    let mut conflicts = 0usize;

    for (reference_index_0based, left_observation) in left {
        let Some(right_observation) = right.get(reference_index_0based) else {
            continue;
        };
        shared_positions += 1;
        if left_observation.reference_base != right_observation.reference_base {
            return Err(Error::Sample(format!(
                "reference base disagrees at position {}",
                reference_index_0based + 1
            )));
        }

        let (Some(left_base), Some(right_base)) =
            (left_observation.query_base, right_observation.query_base)
        else {
            continue;
        };
        comparable_bases += 1;
        if left_base == right_base {
            agreements += 1;
        } else {
            conflicts += 1;
        }
    }

    if shared_positions == 0 {
        return Ok(None);
    }

    let agreement = (comparable_bases > 0).then(|| agreements as f64 / comparable_bases as f64);
    let mut exclusion_reasons = Vec::new();
    if comparable_bases < config.minimum_comparable_bases {
        exclusion_reasons.push(OverlapExclusionReason::ComparableBasesBelowMinimum);
    }
    if agreement.is_some_and(|value| value < config.minimum_overlap_agreement) {
        exclusion_reasons.push(OverlapExclusionReason::AgreementBelowMinimum);
    }

    Ok(Some(ReadOverlapEvidence {
        left_read_index,
        right_read_index,
        shared_positions,
        comparable_bases,
        agreements,
        conflicts,
        agreement,
        eligible: exclusion_reasons.is_empty(),
        exclusion_reasons,
    }))
}

const fn is_canonical(base: char) -> bool {
    matches!(base, 'A' | 'C' | 'G' | 'T')
}

#[cfg(test)]
mod tests {
    use super::*;

    fn config(
        minimum_comparable_bases: usize,
        minimum_overlap_agreement: f64,
    ) -> SampleReconciliationConfig {
        SampleReconciliationConfig {
            minimum_comparable_bases,
            minimum_overlap_agreement,
        }
    }

    fn observation(reference_base: char, query_base: Option<char>) -> CoordinateObservation {
        CoordinateObservation {
            reference_base,
            query_base,
        }
    }

    fn coordinates(
        values: &[(usize, char, Option<char>)],
    ) -> BTreeMap<usize, CoordinateObservation> {
        values
            .iter()
            .map(|(position, reference, query)| (*position, observation(*reference, *query)))
            .collect()
    }

    #[test]
    fn omits_non_overlapping_pairs() -> Result<()> {
        let left = coordinates(&[(1, 'A', Some('A'))]);
        let right = coordinates(&[(2, 'C', Some('C'))]);

        assert!(assess_pair(0, 1, &left, &right, &config(1, 0.5))?.is_none());
        Ok(())
    }

    #[test]
    fn admits_sufficient_canonical_overlap_and_agreement() -> Result<()> {
        let left = coordinates(&[
            (1, 'A', Some('A')),
            (2, 'C', Some('C')),
            (3, 'G', Some('G')),
            (4, 'T', Some('T')),
        ]);
        let right = coordinates(&[
            (1, 'A', Some('A')),
            (2, 'C', Some('C')),
            (3, 'G', Some('G')),
            (4, 'T', Some('A')),
        ]);

        let overlap = assess_pair(0, 1, &left, &right, &config(4, 0.75))?
            .ok_or_else(|| Error::Sample("expected overlap".into()))?;
        assert_eq!(overlap.shared_positions, 4);
        assert_eq!(overlap.comparable_bases, 4);
        assert_eq!(overlap.agreements, 3);
        assert_eq!(overlap.conflicts, 1);
        assert_eq!(overlap.agreement, Some(0.75));
        assert!(overlap.eligible);
        assert!(overlap.exclusion_reasons.is_empty());
        Ok(())
    }

    #[test]
    fn excludes_gaps_and_unresolved_calls_from_nucleotide_agreement() -> Result<()> {
        let left = coordinates(&[(1, 'A', Some('A')), (2, 'C', None), (3, 'G', Some('G'))]);
        let right = coordinates(&[(1, 'A', Some('A')), (2, 'C', Some('C')), (3, 'G', None)]);

        let overlap = assess_pair(0, 1, &left, &right, &config(2, 0.5))?
            .ok_or_else(|| Error::Sample("expected overlap".into()))?;
        assert_eq!(overlap.shared_positions, 3);
        assert_eq!(overlap.comparable_bases, 1);
        assert_eq!(overlap.agreement, Some(1.0));
        assert!(!overlap.eligible);
        assert_eq!(
            overlap.exclusion_reasons,
            vec![OverlapExclusionReason::ComparableBasesBelowMinimum]
        );
        Ok(())
    }

    #[test]
    fn keeps_zero_comparable_overlap_explicit_without_synthetic_agreement() -> Result<()> {
        let left = coordinates(&[(1, 'A', None), (2, 'C', None)]);
        let right = coordinates(&[(1, 'A', Some('A')), (2, 'C', Some('C'))]);

        let overlap = assess_pair(0, 1, &left, &right, &config(1, 0.5))?
            .ok_or_else(|| Error::Sample("expected overlap".into()))?;
        assert_eq!(overlap.shared_positions, 2);
        assert_eq!(overlap.comparable_bases, 0);
        assert_eq!(overlap.agreements, 0);
        assert_eq!(overlap.conflicts, 0);
        assert_eq!(overlap.agreement, None);
        assert!(!overlap.eligible);
        assert_eq!(
            overlap.exclusion_reasons,
            vec![OverlapExclusionReason::ComparableBasesBelowMinimum]
        );
        Ok(())
    }

    #[test]
    fn reports_both_failed_overlap_rules_in_stable_order() -> Result<()> {
        let left = coordinates(&[(1, 'A', Some('A')), (2, 'C', Some('C'))]);
        let right = coordinates(&[(1, 'A', Some('G')), (2, 'C', Some('T'))]);

        let overlap = assess_pair(0, 1, &left, &right, &config(3, 0.5))?
            .ok_or_else(|| Error::Sample("expected overlap".into()))?;
        assert!(!overlap.eligible);
        assert_eq!(
            overlap.exclusion_reasons,
            vec![
                OverlapExclusionReason::ComparableBasesBelowMinimum,
                OverlapExclusionReason::AgreementBelowMinimum,
            ]
        );
        Ok(())
    }
}
