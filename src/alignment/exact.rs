//! Provably exact upper-bound alignment fast path.

use crate::alignment::traceback::RawAlignment;
use crate::config::AlignmentConfig;
use crate::model::locus_evidence::EvidenceProfile;

#[derive(Debug)]
pub(crate) enum UpperBoundPlacement {
    Unproven,
    Unattained,
    Unique(RawAlignment),
    Ambiguous,
}

pub(crate) fn classify(
    _query: &str,
    _profiles: &[Option<EvidenceProfile>],
    _reference: &str,
    _config: &AlignmentConfig,
    _modulo_length: Option<usize>,
) -> UpperBoundPlacement {
    UpperBoundPlacement::Unproven
}

#[cfg(test)]
mod tests {
    use crate::alignment::scoring::SCORE_SCALE;
    use crate::model::locus_evidence::EvidenceProfile;

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

    fn profiles(sequence: &str) -> Vec<Option<EvidenceProfile>> {
        sequence
            .bytes()
            .map(|base| {
                let weights = match base {
                    b'A' => [1.0, 0.0, 0.0, 0.0],
                    b'C' => [0.0, 1.0, 0.0, 0.0],
                    b'G' => [0.0, 0.0, 1.0, 0.0],
                    b'T' => [0.0, 0.0, 0.0, 1.0],
                    _ => [0.0; 4],
                };
                Some(EvidenceProfile { weights })
            })
            .collect()
    }

    fn unique(result: UpperBoundPlacement) -> RawAlignment {
        match result {
            UpperBoundPlacement::Unique(alignment) => alignment,
            other => panic!("expected unique proven placement, got {other:?}"),
        }
    }

    #[test]
    fn proves_unique_gapless_upper_bound_placement() {
        let query = "ACGT";
        let alignment = unique(classify(
            query,
            &profiles(query),
            "TTACGTGG",
            &config(),
            None,
        ));

        assert_eq!(alignment.score, 12 * SCORE_SCALE);
        assert_eq!(alignment.start_reference, 2);
        assert_eq!(alignment.end_reference, 6);
        assert_eq!(alignment.metrics.exact_matches, 4);
        assert_eq!(alignment.metrics.gap_opens, 0);
    }

    #[test]
    fn proof_follows_profile_optimum_not_primary_sequence() {
        let alignment = unique(classify(
            "TCGT",
            &profiles("ACGT"),
            "TTACGTGG",
            &config(),
            None,
        ));

        assert_eq!(alignment.score, 12 * SCORE_SCALE);
        assert_eq!(alignment.start_reference, 2);
        assert_eq!(alignment.metrics.exact_matches, 3);
        assert_eq!(alignment.metrics.mismatches, 1);
    }

    #[test]
    fn reports_multiple_upper_bound_placements_as_ambiguous() {
        assert!(matches!(
            classify("AAA", &profiles("AAA"), "AAAAA", &config(), None),
            UpperBoundPlacement::Ambiguous
        ));
    }

    #[test]
    fn distinguishes_unattained_bound_from_unprovable_bound() {
        assert!(matches!(
            classify("ACGT", &profiles("ACGT"), "TTTT", &config(), None),
            UpperBoundPlacement::Unattained
        ));

        let tied = vec![
            Some(EvidenceProfile {
                weights: [0.5, 0.5, 0.0, 0.0],
            }),
            profiles("C")[0],
        ];
        assert!(matches!(
            classify("AC", &tied, "AC", &config(), None),
            UpperBoundPlacement::Unproven
        ));

        let mut weak = config();
        weak.ambiguous_score = -5;
        weak.gap_extension_score = -1;
        let weak_profiles = vec![Some(EvidenceProfile {
            weights: [0.4, 0.3, 0.2, 0.1],
        })];
        assert!(matches!(
            classify("A", &weak_profiles, "A", &weak, None),
            UpperBoundPlacement::Unproven
        ));
    }

    #[test]
    fn proves_unique_circular_origin_crossing_placement() {
        let reference = "ACGT";
        let working_reference = format!("{reference}{reference}");
        let alignment = unique(classify(
            "GTAC",
            &profiles("GTAC"),
            &working_reference,
            &config(),
            Some(reference.len()),
        ));

        assert_eq!(alignment.start_reference, 2);
        assert_eq!(alignment.end_reference, 6);
        assert_eq!(alignment.metrics.exact_matches, 4);
    }
}
