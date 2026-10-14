//! Equivalence of the exact acceleration tiers with full-reference Gotoh.

use crate::alignment::scoring::SCORE_SCALE;
use dna_kernel::read_evidence::EvidenceProfile;

use crate::alignment::AlignmentConfig;
use crate::alignment::gotoh;
use crate::alignment::traceback::RawAlignment;

use super::*;

type TestResult<T = ()> = Result<T, Box<dyn std::error::Error>>;

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
fn over_limit_problem_falls_back_to_authoritative_gotoh() {
    let query = "A".repeat(2_000);
    let reference = "A".repeat(50_000);
    assert!(matches!(
        classify(&query, &profiles(&query), &reference, &config(), None),
        UpperBoundPlacement::Unproven
    ));
}

type ColumnSignature = (char, char, Option<usize>, Option<usize>);
type AlignmentSignature = (i64, usize, usize, Vec<ColumnSignature>);

fn alignment_signature(alignment: &RawAlignment) -> AlignmentSignature {
    (
        alignment.score,
        alignment.start_reference,
        alignment.end_reference,
        alignment
            .columns
            .iter()
            .map(|column| {
                (
                    column.query_base,
                    column.reference_base,
                    column.query_index,
                    column.reference_index,
                )
            })
            .collect(),
    )
}

fn assert_same_alignments(actual: &[RawAlignment], expected: &[RawAlignment]) {
    let actual = actual.iter().map(alignment_signature).collect::<Vec<_>>();
    let expected = expected.iter().map(alignment_signature).collect::<Vec<_>>();
    assert_eq!(actual, expected);
}

fn require_pruned(value: Pruned, scenario: &str) -> TestResult<Vec<RawAlignment>> {
    match value {
        Pruned::Proven(placements) => Ok(placements),
        Pruned::Unproven { .. } => {
            Err(format!("expected exact seeded pruning for {scenario}").into())
        }
    }
}

#[test]
fn seeded_pruning_matches_gotoh_for_one_snv() -> TestResult {
    let reference_query = "ACGTCAGTACGATCGTACCTGAGTACGA";
    let mut query = reference_query.to_owned();
    query.replace_range(10..11, "T");
    let reference = format!("TTTT{reference_query}CCCC");
    let query_profiles = profiles(&query);

    let expected = gotoh::align(&query, &query_profiles, &reference, &config(), None)?;
    let actual = require_pruned(
        align_pruned(&query, &query_profiles, &reference, &config(), None)?,
        "one-SNV alignment",
    )?;

    assert_same_alignments(&actual, &expected);
    Ok(())
}

#[test]
fn seeded_pruning_matches_gotoh_for_one_insertion() -> TestResult {
    let reference_query = "ACGTCAGTACGATCGTACCTGAGTACGA";
    let query = format!("{}T{}", &reference_query[..12], &reference_query[12..]);
    let reference = format!("TTTT{reference_query}CCCC");
    let query_profiles = profiles(&query);

    let expected = gotoh::align(&query, &query_profiles, &reference, &config(), None)?;
    let actual = require_pruned(
        align_pruned(&query, &query_profiles, &reference, &config(), None)?,
        "one-insertion alignment",
    )?;

    assert_same_alignments(&actual, &expected);
    Ok(())
}

#[test]
fn seeded_pruning_preserves_rightmost_homopolymer_deletion() -> TestResult {
    let reference_query = "GCCAAAAGTTACGTCAGTACGATCGTAC";
    let query = reference_query.replacen("AAAA", "AAA", 1);
    let reference = format!("TTTT{reference_query}CCCC");
    let query_profiles = profiles(&query);

    let expected = gotoh::align(&query, &query_profiles, &reference, &config(), None)?;
    let actual = require_pruned(
        align_pruned(&query, &query_profiles, &reference, &config(), None)?,
        "homopolymer deletion alignment",
    )?;

    assert_same_alignments(&actual, &expected);
    let deleted = actual[0]
        .columns
        .iter()
        .find(|column| column.query_base == '-')
        .and_then(|column| column.reference_index);
    assert_eq!(deleted, Some(10));
    Ok(())
}

#[test]
fn seeded_pruning_preserves_distinct_equal_placements() -> TestResult {
    let motif = "ACGTCAGTACGATCGTACCTGAGTACGA";
    let mut query = motif.to_owned();
    query.replace_range(10..11, "T");
    let reference = format!("GG{motif}TT{motif}CC");
    let query_profiles = profiles(&query);

    let expected = gotoh::align(&query, &query_profiles, &reference, &config(), None)?;
    let actual = require_pruned(
        align_pruned(&query, &query_profiles, &reference, &config(), None)?,
        "repeated SNV alignment",
    )?;

    assert_eq!(expected.len(), 2);
    assert_same_alignments(&actual, &expected);
    Ok(())
}

#[test]
fn seeded_pruning_supports_circular_non_origin_windows() -> TestResult {
    let reference = "ACGTCAGTACGATCGTACCTGAGTACGATTTTGGGGCCCCAAAATTTT";
    let reference_query = &reference[8..36];
    let mut query = reference_query.to_owned();
    query.replace_range(10..11, "A");
    let working_reference = format!("{reference}{reference}");
    let query_profiles = profiles(&query);

    let expected = gotoh::align(
        &query,
        &query_profiles,
        &working_reference,
        &config(),
        Some(reference.len()),
    )?;
    let actual = require_pruned(
        align_pruned(
            &query,
            &query_profiles,
            &working_reference,
            &config(),
            Some(reference.len()),
        )?,
        "non-origin circular alignment",
    )?;

    assert_same_alignments(&actual, &expected);
    Ok(())
}

fn reverse_complement(sequence: &str) -> String {
    sequence
        .bytes()
        .rev()
        .map(|base| match base {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => 'N',
        })
        .collect()
}

#[test]
fn threshold_proof_can_exclude_opposite_orientation() -> TestResult {
    let reference_query = "ACGTCAGTACGATCGTACCTGAGTACGA";
    let mut query = reference_query.to_owned();
    query.replace_range(10..11, "T");
    let reference = format!("TTTT{reference_query}CCCC");
    let query_profiles = profiles(&query);
    let forward = require_pruned(
        align_pruned(&query, &query_profiles, &reference, &config(), None)?,
        "forward threshold source",
    )?;
    let threshold = forward[0].score;

    let reverse_query = reverse_complement(&query);
    let reverse_profiles = profiles(&reverse_query);
    let reverse = align_at_or_above(
        &reverse_query,
        &reverse_profiles,
        &reference,
        &config(),
        None,
        threshold,
    )?;

    assert!(reverse.is_some_and(|placements| placements.is_empty()));
    Ok(())
}

#[test]
fn seeded_pruning_matches_gotoh_across_circular_origin() -> TestResult {
    let reference = "ACGTCAGTACGATCGTACCTGAGTACGATTGGCCAATGCATGCAAGTTCCAGTCGATGGCATCTAGGACT";
    let length = reference.len();
    let working_reference = format!("{reference}{reference}");
    for start in [length - 20, length - 3, 2] {
        let mut query = working_reference[start..start + 30].to_owned();
        query.replace_range(12..13, if &query[12..13] == "A" { "C" } else { "A" });
        let query_profiles = profiles(&query);

        let expected = gotoh::align(
            &query,
            &query_profiles,
            &working_reference,
            &config(),
            Some(length),
        )?;
        let actual = require_pruned(
            align_pruned(
                &query,
                &query_profiles,
                &working_reference,
                &config(),
                Some(length),
            )?,
            "origin-adjacent circular alignment",
        )?;
        assert_same_alignments(&actual, &expected);
    }
    Ok(())
}

#[test]
fn seeded_pruning_matches_gotoh_with_unresolved_and_mixed_loci() -> TestResult {
    let reference_query = "ACGTCAGTACGATCGTACCTGAGTACGATTGGCCAATGCATG";
    let query = reference_query.to_owned();
    let mut query_profiles = profiles(&query);
    query_profiles[5] = None;
    query_profiles[17] = Some(EvidenceProfile {
        weights: [0.5, 0.5, 0.0, 0.0],
    });
    query_profiles[30] = Some(EvidenceProfile {
        weights: [0.3, 0.3, 0.2, 0.2],
    });
    let reference = format!("TTTT{reference_query}CCCC");

    let expected = gotoh::align(&query, &query_profiles, &reference, &config(), None)?;
    let actual = require_pruned(
        align_pruned(&query, &query_profiles, &reference, &config(), None)?,
        "alignment with unresolved and mixed loci",
    )?;

    assert_same_alignments(&actual, &expected);
    Ok(())
}

/// Deterministic xorshift generator for the randomized equivalence check.
struct Generator(u64);

impl Generator {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn below(&mut self, bound: usize) -> usize {
        usize::try_from(self.next() % u64::try_from(bound).unwrap_or(1)).unwrap_or(0)
    }
}

fn random_profile(generator: &mut Generator, base: u8) -> Option<EvidenceProfile> {
    let index = match base {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => return None,
    };
    let mut weights = [0.0; 4];
    match generator.below(10) {
        0 => return None,
        1 => {
            weights[index] = 0.5;
            weights[(index + 1) % 4] = 0.5;
        }
        2 => {
            let first = generator.below(1_000) as f64 / 1_000.0;
            weights[index] = first;
            weights[(index + 2) % 4] = 1.0 - first;
        }
        _ => weights[index] = 1.0,
    }
    Some(EvidenceProfile { weights })
}

/// Randomized reads with substitutions, indels, unresolved and mixed loci,
/// against linear and circular references with an `N`, must prune to exactly
/// the whole-reference Gotoh placements whenever a proof applies.
#[test]
fn seeded_pruning_matches_gotoh_on_random_reads() -> TestResult {
    let mut generator = Generator(0x9e37_79b9_7f4a_7c15);
    let mut pruned = 0;
    for case in 0..400 {
        let length = 200 + generator.below(200);
        let mut reference = (0..length)
            .map(|_| char::from(b"ACGT"[generator.below(4)]))
            .collect::<String>();
        if case % 3 == 0 {
            let position = generator.below(length);
            reference.replace_range(position..=position, "N");
        }
        let circular = case % 2 == 0;
        let working_reference = if circular {
            format!("{reference}{reference}")
        } else {
            reference.clone()
        };
        let modulo_length = circular.then_some(length);
        let span = 40 + generator.below(80);
        let start = if circular {
            generator.below(length)
        } else {
            generator.below(length - span)
        };
        let mut read = working_reference.as_bytes()[start..start + span].to_vec();
        for _ in 0..generator.below(4) {
            let position = generator.below(read.len());
            match generator.below(3) {
                0 => read[position] = b"ACGT"[generator.below(4)],
                1 => read.insert(position, b"ACGT"[generator.below(4)]),
                _ => {
                    read.remove(position);
                }
            }
        }
        let query_profiles = read
            .iter()
            .map(|base| random_profile(&mut generator, *base))
            .collect::<Vec<_>>();
        let query = read
            .iter()
            .zip(&query_profiles)
            .map(|(base, profile)| {
                if profile.is_some() {
                    char::from(*base)
                } else {
                    'N'
                }
            })
            .collect::<String>();

        let Ok(expected) = gotoh::align(
            &query,
            &query_profiles,
            &working_reference,
            &config(),
            modulo_length,
        ) else {
            continue;
        };
        if let Pruned::Proven(actual) = align_pruned(
            &query,
            &query_profiles,
            &working_reference,
            &config(),
            modulo_length,
        )? {
            pruned += 1;
            assert_eq!(actual.len(), expected.len(), "case {case}");
            if expected.len() == 1 {
                assert_same_alignments(&actual, &expected);
            }
        }
        if let Some(actual) = align_at_or_above(
            &query,
            &query_profiles,
            &working_reference,
            &config(),
            modulo_length,
            expected[0].score + 1,
        )? {
            assert!(actual.is_empty(), "case {case}");
        }
    }
    assert!(pruned > 200, "only {pruned} random cases were pruned");
    Ok(())
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
