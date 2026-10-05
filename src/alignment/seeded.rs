//! Score-bounded seeded alignment proof.

use crate::alignment::traceback::RawAlignment;
use crate::config::AlignmentConfig;
use crate::model::locus_evidence::EvidenceProfile;

#[derive(Debug)]
pub(crate) struct SeededProof {
    pub(crate) score: i64,
    pub(crate) placements: Vec<RawAlignment>,
}

pub(crate) fn classify(
    _query: &str,
    _profiles: &[Option<EvidenceProfile>],
    _reference: &str,
    _config: &AlignmentConfig,
    _modulo_length: Option<usize>,
) -> Option<SeededProof> {
    None
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

    #[test]
    fn proves_unique_one_substitution_placement() {
        let Some(proof) = classify("ACGT", &profiles("ACGT"), "TTACATGG", &config(), None) else {
            panic!("one-substitution placement should be provable");
        };

        assert_eq!(proof.score, 4 * SCORE_SCALE);
        assert_eq!(proof.placements.len(), 1);
        assert_eq!(proof.placements[0].start_reference, 2);
        assert_eq!(proof.placements[0].end_reference, 6);
        assert_eq!(proof.placements[0].metrics.exact_matches, 3);
        assert_eq!(proof.placements[0].metrics.mismatches, 1);
        assert_eq!(proof.placements[0].metrics.gap_opens, 0);
    }

    #[test]
    fn proves_all_equally_best_repeated_substitution_placements() {
        let Some(proof) = classify("ACGT", &profiles("ACGT"), "ACATGGACAT", &config(), None) else {
            panic!("repeated best placements should be completely enumerated");
        };

        assert_eq!(proof.score, 4 * SCORE_SCALE);
        assert_eq!(
            proof
                .placements
                .iter()
                .map(|alignment| alignment.start_reference)
                .collect::<Vec<_>>(),
            vec![0, 6]
        );
    }

    #[test]
    fn refuses_candidate_when_a_gap_can_rival_its_score() {
        let mut weak_gap = config();
        weak_gap.gap_open_score = -1;
        weak_gap.gap_extension_score = -1;

        assert!(classify("ACGT", &profiles("ACGT"), "TTACATGG", &weak_gap, None,).is_none());
    }

    #[test]
    fn refuses_candidate_outside_two_substitution_seed_bound() {
        assert!(classify("ACGTAC", &profiles("ACGTAC"), "TTATATTCGG", &config(), None,).is_none());
    }

    #[test]
    fn proves_one_substitution_across_circular_origin() {
        let reference = "ACGT";
        let working_reference = format!("{reference}{reference}");
        let Some(proof) = classify(
            "GTTC",
            &profiles("GTTC"),
            &working_reference,
            &config(),
            Some(reference.len()),
        ) else {
            panic!("origin-crossing substitution should be provable");
        };

        assert_eq!(proof.placements.len(), 1);
        assert_eq!(proof.placements[0].start_reference, 2);
        assert_eq!(proof.placements[0].end_reference, 6);
        assert_eq!(proof.placements[0].metrics.mismatches, 1);
    }
}
