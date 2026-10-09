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

fn require_pruned(
    value: Option<Vec<RawAlignment>>,
    scenario: &str,
) -> TestResult<Vec<RawAlignment>> {
    value.ok_or_else(|| format!("expected exact seeded pruning for {scenario}").into())
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
fn seeded_pruning_falls_back_for_circular_origin_crossing() -> TestResult {
    let reference = "ACGTCAGTACGATCGTACCTGAGTACGA";
    let mut query = format!("{}{}", &reference[18..], &reference[..18]);
    query.replace_range(5..6, "A");
    let working_reference = format!("{reference}{reference}");
    let query_profiles = profiles(&query);

    assert!(
        align_pruned(
            &query,
            &query_profiles,
            &working_reference,
            &config(),
            Some(reference.len()),
        )?
        .is_none()
    );
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
