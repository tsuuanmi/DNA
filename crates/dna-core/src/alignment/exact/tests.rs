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

fn upper_bound_on(
    query: &str,
    query_profiles: &[Option<EvidenceProfile>],
    reference: &str,
    config: &AlignmentConfig,
    modulo_length: Option<usize>,
) -> UpperBoundPlacement {
    let reference = WorkingReference::new(reference, modulo_length);
    Problem::new(query, query_profiles, &reference, config)
        .map_or(UpperBoundPlacement::Unproven, |problem| classify(&problem))
}

fn pruned_on(
    query: &str,
    query_profiles: &[Option<EvidenceProfile>],
    reference: &str,
    config: &AlignmentConfig,
    modulo_length: Option<usize>,
) -> dna_kernel::error::Result<Pruned> {
    let reference = WorkingReference::new(reference, modulo_length);
    match Problem::new(query, query_profiles, &reference, config) {
        Some(problem) => align_pruned(&problem),
        None => Ok(Pruned::Unproven { lower_bound: None }),
    }
}

fn at_or_above_on(
    query: &str,
    query_profiles: &[Option<EvidenceProfile>],
    reference: &str,
    config: &AlignmentConfig,
    modulo_length: Option<usize>,
    threshold: i64,
) -> dna_kernel::error::Result<Option<Vec<RawAlignment>>> {
    let reference = WorkingReference::new(reference, modulo_length);
    match Problem::new(query, query_profiles, &reference, config) {
        Some(problem) => align_at_or_above(&problem, threshold),
        None => Ok(None),
    }
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
    let alignment = unique(upper_bound_on(
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
    let alignment = unique(upper_bound_on(
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
        upper_bound_on("AAA", &profiles("AAA"), "AAAAA", &config(), None),
        UpperBoundPlacement::Ambiguous
    ));
}

#[test]
fn distinguishes_unattained_bound_from_unprovable_bound() {
    assert!(matches!(
        upper_bound_on("ACGT", &profiles("ACGT"), "TTTT", &config(), None),
        UpperBoundPlacement::Unattained
    ));

    let tied = vec![
        Some(EvidenceProfile {
            weights: [0.5, 0.5, 0.0, 0.0],
        }),
        profiles("C")[0],
    ];
    assert!(matches!(
        upper_bound_on("AC", &tied, "AC", &config(), None),
        UpperBoundPlacement::Unproven
    ));

    let mut weak = config();
    weak.ambiguous_score = -5;
    weak.gap_extension_score = -1;
    let weak_profiles = vec![Some(EvidenceProfile {
        weights: [0.4, 0.3, 0.2, 0.1],
    })];
    assert!(matches!(
        upper_bound_on("A", &weak_profiles, "A", &weak, None),
        UpperBoundPlacement::Unproven
    ));
}

#[test]
fn over_limit_problem_falls_back_to_authoritative_gotoh() {
    let query = "A".repeat(2_000);
    let reference = "A".repeat(50_000);
    assert!(matches!(
        upper_bound_on(&query, &profiles(&query), &reference, &config(), None),
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
        pruned_on(&query, &query_profiles, &reference, &config(), None)?,
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
        pruned_on(&query, &query_profiles, &reference, &config(), None)?,
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
        pruned_on(&query, &query_profiles, &reference, &config(), None)?,
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
        pruned_on(&query, &query_profiles, &reference, &config(), None)?,
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
        pruned_on(
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

/// The seed lower bound is the score of a valid alignment, so a complete
/// proof at that threshold that finds no placement is an internal fault, not a
/// reason to fall back.
#[test]
fn empty_proof_at_the_seed_lower_bound_is_inconsistent() {
    assert!(matches!(
        seeded::proven_at_lower_bound(Some(Vec::new()), 7),
        Err(dna_kernel::error::Error::Alignment(
            dna_kernel::error::AlignmentError::Inconsistent(_)
        ))
    ));
    assert!(matches!(
        seeded::proven_at_lower_bound(None, 7),
        Ok(Pruned::Unproven {
            lower_bound: Some(7)
        })
    ));
}

#[test]
fn threshold_proof_can_exclude_opposite_orientation() -> TestResult {
    let reference_query = "ACGTCAGTACGATCGTACCTGAGTACGA";
    let mut query = reference_query.to_owned();
    query.replace_range(10..11, "T");
    let reference = format!("TTTT{reference_query}CCCC");
    let query_profiles = profiles(&query);
    let forward = require_pruned(
        pruned_on(&query, &query_profiles, &reference, &config(), None)?,
        "forward threshold source",
    )?;
    let threshold = forward[0].score;

    let reverse_query = reverse_complement(&query);
    let reverse_profiles = profiles(&reverse_query);
    let reverse = at_or_above_on(
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
            pruned_on(
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

/// A read spanning the whole circle whose last probe seed lands just after the
/// origin needs a probe window reaching further back than one circle, which
/// must start at the doubled reference's beginning rather than underflow.
#[test]
fn seeded_pruning_clamps_a_window_reaching_past_the_doubled_reference_start() -> TestResult {
    let mut generator = Generator(0x05ee_d0fc_1c1e);
    let reference = (0..100)
        .map(|_| char::from(b"ACGT"[generator.below(4)]))
        .collect::<String>();
    let working_reference = format!("{reference}{reference}");
    // The last of four probe seeds starts at query locus 79 with a margin of
    // 25, so its window reaches 104 bases back from reference position 1.
    let start = 22;
    let mut query = working_reference[start..start + 100].to_owned();
    query.replace_range(10..11, if &query[10..11] == "A" { "C" } else { "A" });
    let query_profiles = profiles(&query);
    let seed = &query[79..95];
    assert_eq!(working_reference[..100].matches(seed).count(), 1);
    assert_eq!(working_reference.find(seed), Some(1));

    let expected = gotoh::align(
        &query,
        &query_profiles,
        &working_reference,
        &config(),
        Some(100),
    )?;
    let actual = require_pruned(
        pruned_on(
            &query,
            &query_profiles,
            &working_reference,
            &config(),
            Some(100),
        )?,
        "whole-circle rotation",
    )?;
    assert_same_alignments(&actual, &expected);
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
        pruned_on(&query, &query_profiles, &reference, &config(), None)?,
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
/// against linear and circular references with an `N` or a repeated unit,
/// must prune to exactly the whole-reference Gotoh placements whenever a proof
/// applies, at the seed lower bound and at explicit thresholds.
#[test]
fn seeded_pruning_matches_gotoh_on_random_reads() -> TestResult {
    const UNIT: usize = 130;
    let mut generator = Generator(0x9e37_79b9_7f4a_7c15);
    let mut pruned = 0;
    let mut ambiguous = 0;
    let mut thresholded = 0;
    for case in 0..400 {
        let repeated = case % 4 == 1;
        let length = 2 * UNIT + generator.below(140);
        let mut reference = (0..length)
            .map(|_| char::from(b"ACGT"[generator.below(4)]))
            .collect::<String>();
        if repeated {
            // The reference ends with a copy of its first unit, and the read
            // comes from that unit, so it has two equal placements.
            let unit = reference[..UNIT].to_owned();
            reference.replace_range(length - UNIT.., &unit);
        }
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
        let span = 40 + generator.below(UNIT - 50);
        let start = if repeated {
            generator.below(UNIT - span)
        } else if circular {
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
        if let Pruned::Proven(actual) = pruned_on(
            &query,
            &query_profiles,
            &working_reference,
            &config(),
            modulo_length,
        )? {
            pruned += 1;
            ambiguous += usize::from(expected.len() == 2);
            assert_same_alignments(&actual, &expected);
        }
        // At or below the optimum a proof yields exactly Gotoh's placements;
        // above it, a proof yields none.
        let optimum = expected[0].score;
        for threshold in [
            optimum - 5 * SCORE_SCALE,
            optimum - SCORE_SCALE,
            optimum - 1,
            optimum,
        ] {
            if let Some(actual) = at_or_above_on(
                &query,
                &query_profiles,
                &working_reference,
                &config(),
                modulo_length,
                threshold,
            )? {
                thresholded += 1;
                assert_same_alignments(&actual, &expected);
            }
        }
        for threshold in [optimum + 1, optimum + SCORE_SCALE] {
            if let Some(actual) = at_or_above_on(
                &query,
                &query_profiles,
                &working_reference,
                &config(),
                modulo_length,
                threshold,
            )? {
                assert!(actual.is_empty(), "case {case}");
            }
        }
    }
    assert!(
        ambiguous > 20,
        "only {ambiguous} pruned cases were ambiguous"
    );
    assert!(
        thresholded > 800,
        "only {thresholded} thresholds at or below the optimum were proven"
    );
    assert!(pruned > 200, "only {pruned} random cases were pruned");
    Ok(())
}

#[test]
fn proves_unique_circular_origin_crossing_placement() {
    let reference = "ACGT";
    let working_reference = format!("{reference}{reference}");
    let alignment = unique(upper_bound_on(
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
