//! Sample aggregation over synthetic read observations.

use crate::model::alignment::{
    Alignment, AlignmentColumn, AlignmentMetrics, Orientation, ReferenceSegment,
};
use crate::model::sample_evidence::SampleEvidence;
use crate::model::variant::{
    ObservedVariant, Variant, VariantCallMapping, VariantCallRole, VariantCallingResult,
    VariantExclusionReason, VariantKind,
};
use dna_kernel::read_evidence::{
    CallEvidence, CallMask, EvidenceProfile, EvidenceReason, MaskedAlignment, ReadEvidence, VetoSet,
};

use super::*;

type TestResult<T = ()> = std::result::Result<T, Box<dyn std::error::Error>>;

fn sample_config() -> SampleReconciliationConfig {
    SampleReconciliationConfig {
        minimum_comparable_bases: 1,
        minimum_overlap_agreement: 0.5,
    }
}

/// Evidence of a read whose every call is an unmasked `G` with one profile.
fn evidence(call_count: usize) -> ReadEvidence {
    let calls = (0..call_count).map(|_| call('G', Some(PROFILE))).collect();
    ReadEvidence::new(calls, 0..call_count, Vec::new())
        .unwrap_or_else(|error| panic!("valid synthetic evidence: {error}"))
}

const PROFILE: [f64; 4] = [0.1, 0.2, 0.3, 0.4];

fn call(base: char, profile: Option<[f64; 4]>) -> CallEvidence {
    CallEvidence {
        base,
        profile: profile.map(|weights| EvidenceProfile { weights }),
        mask: None,
        vetoes: VetoSet::default(),
    }
}

/// Runs aggregation over owned reads with no rejected read.
fn run(reads: &[CalledRead]) -> Result<SampleEvidence> {
    aggregate(&reads.iter().collect::<Vec<_>>(), &[], &sample_config())
}

fn observation(
    input_sha256: &str,
    reference_sha256: &str,
    configuration_sha256: &str,
    orientation: Orientation,
    columns: Vec<AlignmentColumn>,
    variants: Vec<Variant>,
) -> CalledRead {
    let call_count = columns
        .iter()
        .filter_map(|column| column.original_call_index_0based)
        .max()
        .map_or(0, |index| index + 1);
    let observed = variants
        .iter()
        .cloned()
        .map(|variant| ObservedVariant {
            variant,
            exclusion_reasons: Vec::new(),
        })
        .collect();
    let reference_segments = columns
        .iter()
        .filter_map(|column| column.reference_index_0based)
        .fold(None::<(usize, usize)>, |bounds, index| {
            Some(match bounds {
                Some((start, end)) => (start.min(index), end.max(index + 1)),
                None => (index, index + 1),
            })
        })
        .map_or_else(Vec::new, |(start_0based, end_0based_exclusive)| {
            vec![ReferenceSegment {
                start_0based,
                end_0based_exclusive,
            }]
        });
    CalledRead {
        input_name: format!("{input_sha256}.ab1"),
        input_sha256: input_sha256.into(),
        reference_sha256: reference_sha256.into(),
        configuration_sha256: configuration_sha256.into(),
        evidence: evidence(call_count),
        alignment: Alignment {
            orientation,
            score: 1,
            callable_segments: reference_segments.clone(),
            reference_segments,
            wraps_origin: false,
            metrics: AlignmentMetrics {
                exact_matches: 0,
                mismatches: 0,
                gap_opens: 0,
                callable_columns: 0,
                callable_identity: 0.0,
                unresolved_query_bases: 0,
                masked_query_bases: 0,
            },
            columns,
        },
        variants: VariantCallingResult {
            reported: variants,
            observed,
            excluded: Vec::new(),
        },
    }
}

fn column(
    query_base: char,
    reference_base: char,
    call: Option<usize>,
    reference: usize,
) -> AlignmentColumn {
    AlignmentColumn {
        query_base,
        reference_base,
        original_call_index_0based: call,
        reference_index_0based: Some(reference),
    }
}

fn snv(position_1based: usize, reference: &str, alternate: &str) -> Variant {
    Variant {
        contig: "reference".into(),
        position_1based,
        reference: reference.into(),
        alternate: alternate.into(),
        kind: VariantKind::Snv,
        calls: vec![VariantCallMapping {
            role: VariantCallRole::Supporting,
            call_index_0based: 0,
            reference_position_0based: Some(position_1based - 1),
        }],
    }
}

#[test]
fn orders_reads_once_and_factors_read_identity_from_evidence() -> TestResult {
    let forward = observation(
        "a",
        "reference",
        "config",
        Orientation::Forward,
        vec![column('G', 'A', Some(0), 72)],
        vec![snv(73, "A", "G")],
    );
    let reverse = observation(
        "b",
        "reference",
        "config",
        Orientation::Reverse,
        vec![column('G', 'A', Some(0), 72)],
        vec![snv(73, "A", "G")],
    );

    let evidence = run(&[reverse, forward])?;

    assert_eq!(evidence.reads[0].input_name, "a.ab1");
    assert_eq!(evidence.reads[1].input_name, "b.ab1");
    assert_eq!(evidence.coverage.len(), 1);
    assert_eq!(evidence.coverage[0].start_0based, 72);
    assert_eq!(evidence.coverage[0].end_0based_exclusive, 73);
    assert_eq!(evidence.coverage[0].read_depth, 2);
    assert_eq!(evidence.coverage[0].forward_depth, 1);
    assert_eq!(evidence.coverage[0].reverse_depth, 1);
    assert_eq!(evidence.overlaps.len(), 1);
    assert_eq!(evidence.overlaps[0].left_read_index, 0);
    assert_eq!(evidence.overlaps[0].right_read_index, 1);
    assert_eq!(evidence.overlaps[0].comparable_bases, 1);
    assert!(evidence.overlaps[0].eligible);
    assert_eq!(evidence.locus_differences.len(), 1);
    assert_eq!(evidence.locus_differences[0].support_topology.reads, 2);
    assert_eq!(
        evidence.locus_differences[0].support_topology.forward_reads,
        1
    );
    assert_eq!(
        evidence.locus_differences[0].support_topology.reverse_reads,
        1
    );
    assert_eq!(
        evidence.locus_differences[0]
            .support_topology
            .alternate_reads,
        2
    );
    assert_eq!(
        evidence.locus_differences[0]
            .support_topology
            .reference_reads,
        0
    );
    assert_eq!(
        evidence.locus_differences[0].support_topology.profile_reads,
        2
    );
    assert_eq!(
        evidence.locus_differences[0]
            .support_topology
            .profile_forward_reads,
        1
    );
    assert_eq!(
        evidence.locus_differences[0]
            .support_topology
            .profile_reverse_reads,
        1
    );
    let nucleotide_support = evidence.locus_differences[0].nucleotide_support;
    assert_eq!(nucleotide_support.contributors, 2);
    assert_eq!(nucleotide_support.forward_contributors, 1);
    assert_eq!(nucleotide_support.reverse_contributors, 1);
    assert_eq!(nucleotide_support.support, [0.5, 0.5, 0.5, 0.5]);
    assert_eq!(nucleotide_support.forward_support, [0.1, 0.2, 0.3, 0.4]);
    assert_eq!(nucleotide_support.reverse_support, [0.4, 0.3, 0.2, 0.1]);
    assert_eq!(
        nucleotide_support
            .mean_profile
            .map(|profile| profile.weights),
        Some([0.25, 0.25, 0.25, 0.25])
    );
    assert_eq!(
        nucleotide_support
            .forward_mean_profile
            .map(|profile| profile.weights),
        Some([0.1, 0.2, 0.3, 0.4])
    );
    assert_eq!(
        nucleotide_support
            .reverse_mean_profile
            .map(|profile| profile.weights),
        Some([0.4, 0.3, 0.2, 0.1])
    );
    let heterogeneity = nucleotide_support
        .heterogeneity
        .ok_or("total profile heterogeneity is missing")?;
    assert!((heterogeneity.within_profile_impurity - 0.7).abs() < 1e-12);
    assert!((heterogeneity.between_profile_dispersion - 0.05).abs() < 1e-12);
    assert!((heterogeneity.total - 0.75).abs() < 1e-12);
    assert_eq!(
        nucleotide_support
            .forward_heterogeneity
            .map(|geometry| geometry.between_profile_dispersion),
        Some(0.0)
    );
    assert_eq!(
        nucleotide_support
            .reverse_heterogeneity
            .map(|geometry| geometry.between_profile_dispersion),
        Some(0.0)
    );
    assert!(
        (nucleotide_support
            .directional_profile_distance
            .ok_or("directional profile distance is missing")?
            - 0.4)
            .abs()
            < 1e-12
    );
    assert_eq!(evidence.locus_differences[0].observations.len(), 2);
    assert_eq!(evidence.locus_differences[0].observations[0].read_index, 0);
    let forward_observation = &evidence.locus_differences[0].observations[0];
    assert_eq!(forward_observation.call_index_0based, Some(0));
    assert_eq!(
        forward_observation.profile.map(|profile| profile.weights),
        Some([0.1, 0.2, 0.3, 0.4])
    );
    assert_eq!(
        forward_observation.nucleotide_contribution,
        crate::model::sample_evidence::NucleotideContribution::Eligible
    );
    assert_eq!(evidence.locus_differences[0].observations[1].read_index, 1);
    let reverse_observation = &evidence.locus_differences[0].observations[1];
    assert_eq!(reverse_observation.call_index_0based, Some(0));
    assert_eq!(
        reverse_observation.profile.map(|profile| profile.weights),
        Some([0.4, 0.3, 0.2, 0.1])
    );
    assert_eq!(evidence.variants.len(), 1);
    assert_eq!(evidence.variants[0].support_topology.reads, 2);
    assert_eq!(evidence.variants[0].support_topology.eligible_reads, 2);
    assert_eq!(evidence.variants[0].support_topology.forward_reads, 1);
    assert_eq!(evidence.variants[0].support_topology.reverse_reads, 1);
    assert_eq!(
        evidence.variants[0].support_topology.eligible_forward_reads,
        1
    );
    assert_eq!(
        evidence.variants[0].support_topology.eligible_reverse_reads,
        1
    );
    assert_eq!(evidence.variants[0].support[0].read_index, 0);
    assert_eq!(
        evidence.variants[0].support[0].calls[0].call_index_0based,
        0
    );
    assert_eq!(evidence.variants[0].support[0].calls[0].base, 'G');
    assert_eq!(evidence.variants[0].support[1].read_index, 1);
    assert_eq!(
        evidence.variants[0].support[1].calls[0].call_index_0based,
        0
    );
    assert_eq!(evidence.variants[0].support[1].calls[0].base, 'C');
    Ok(())
}

#[test]
fn omits_reference_matches_but_preserves_non_reference_states() -> Result<()> {
    let mut read = observation(
        "a",
        "reference",
        "config",
        Orientation::Forward,
        vec![
            column('A', 'A', Some(0), 10),
            column('G', 'C', Some(1), 11),
            column('N', 'T', Some(2), 12),
            column('-', 'G', None, 13),
            column('C', 'C', Some(3), 14),
        ],
        Vec::new(),
    );
    read.evidence = read.evidence.with_call(2, call('N', None));

    let evidence = run(&[read])?;

    assert_eq!(evidence.locus_differences.len(), 3);
    assert_eq!(evidence.locus_differences[0].position_1based, 12);
    assert_eq!(evidence.locus_differences[0].support_topology.reads, 1);
    assert_eq!(
        evidence.locus_differences[0]
            .support_topology
            .alternate_reads,
        1
    );
    assert_eq!(
        evidence.locus_differences[1]
            .support_topology
            .unresolved_reads,
        1
    );
    assert_eq!(
        evidence.locus_differences[2]
            .support_topology
            .deletion_reads,
        1
    );
    assert_eq!(
        evidence.locus_differences[0].support_topology.profile_reads,
        1
    );
    assert_eq!(
        evidence.locus_differences[1].support_topology.profile_reads,
        0
    );
    assert_eq!(
        evidence.locus_differences[2].support_topology.profile_reads,
        0
    );
    assert_eq!(
        evidence.locus_differences[0]
            .nucleotide_support
            .contributors,
        1
    );
    assert_eq!(
        evidence.locus_differences[1]
            .nucleotide_support
            .contributors,
        0
    );
    assert_eq!(
        evidence.locus_differences[1].nucleotide_support.support,
        [0.0; 4]
    );
    assert_eq!(
        evidence.locus_differences[2]
            .nucleotide_support
            .contributors,
        0
    );
    assert_eq!(
        evidence.locus_differences[2].nucleotide_support.support,
        [0.0; 4]
    );
    assert!(
        evidence.locus_differences[1]
            .nucleotide_support
            .mean_profile
            .is_none()
    );
    assert!(
        evidence.locus_differences[2]
            .nucleotide_support
            .mean_profile
            .is_none()
    );
    assert!(
        evidence.locus_differences[1]
            .nucleotide_support
            .heterogeneity
            .is_none()
    );
    assert!(
        evidence.locus_differences[2]
            .nucleotide_support
            .heterogeneity
            .is_none()
    );
    assert!(
        evidence.locus_differences[1]
            .nucleotide_support
            .directional_profile_distance
            .is_none()
    );
    assert!(
        evidence.locus_differences[2]
            .nucleotide_support
            .directional_profile_distance
            .is_none()
    );
    assert_eq!(
        evidence.locus_differences[0].observations[0].state,
        crate::model::sample_evidence::LocusState::Alternate
    );
    assert_eq!(
        evidence.locus_differences[1].observations[0].state,
        crate::model::sample_evidence::LocusState::Unresolved
    );
    assert_eq!(
        evidence.locus_differences[1].observations[0].call_index_0based,
        Some(2)
    );
    assert!(
        evidence.locus_differences[1].observations[0]
            .profile
            .is_none()
    );
    assert_eq!(
        evidence.locus_differences[1].observations[0].nucleotide_contribution,
        crate::model::sample_evidence::NucleotideContribution::MissingProfile
    );
    assert_eq!(
        evidence.locus_differences[2].observations[0].state,
        crate::model::sample_evidence::LocusState::Deletion
    );
    assert!(
        evidence.locus_differences[2].observations[0]
            .call_index_0based
            .is_none()
    );
    assert!(
        evidence.locus_differences[2].observations[0]
            .profile
            .is_none()
    );
    assert_eq!(
        evidence.locus_differences[2].observations[0].nucleotide_contribution,
        crate::model::sample_evidence::NucleotideContribution::DeletionEvent
    );
    Ok(())
}

#[test]
fn unresolved_call_with_profile_remains_nucleotide_eligible() -> TestResult {
    let read = observation(
        "a",
        "reference",
        "config",
        Orientation::Forward,
        vec![column('N', 'T', Some(0), 12)],
        Vec::new(),
    );

    let evidence = run(&[read])?;

    assert_eq!(evidence.locus_differences.len(), 1);
    assert_eq!(
        evidence.locus_differences[0].observations[0].state,
        crate::model::sample_evidence::LocusState::Unresolved
    );
    assert_eq!(
        evidence.locus_differences[0].support_topology.profile_reads,
        1
    );
    assert_eq!(
        evidence.locus_differences[0].observations[0].nucleotide_contribution,
        crate::model::sample_evidence::NucleotideContribution::Eligible
    );
    assert_eq!(
        evidence.locus_differences[0]
            .nucleotide_support
            .contributors,
        1
    );
    assert_eq!(
        evidence.locus_differences[0].nucleotide_support.support,
        [0.1, 0.2, 0.3, 0.4]
    );
    assert_eq!(
        evidence.locus_differences[0]
            .nucleotide_support
            .mean_profile
            .map(|profile| profile.weights),
        Some([0.1, 0.2, 0.3, 0.4])
    );
    let heterogeneity = evidence.locus_differences[0]
        .nucleotide_support
        .heterogeneity
        .ok_or("single-read profile heterogeneity is missing")?;
    assert!((heterogeneity.within_profile_impurity - 0.7).abs() < 1e-12);
    assert_eq!(heterogeneity.between_profile_dispersion, 0.0);
    assert!((heterogeneity.total - 0.7).abs() < 1e-12);
    assert!(
        evidence.locus_differences[0]
            .nucleotide_support
            .directional_profile_distance
            .is_none()
    );
    Ok(())
}

#[test]
fn masked_calls_never_retain_a_locus_but_are_kept_where_another_read_differs() -> TestResult {
    let differing = observation(
        "a",
        "reference",
        "config",
        Orientation::Forward,
        vec![column('A', 'T', Some(0), 12), column('C', 'C', Some(1), 13)],
        Vec::new(),
    );
    let mut masked = observation(
        "b",
        "reference",
        "config",
        Orientation::Forward,
        vec![column('N', 'T', Some(0), 12), column('N', 'C', Some(1), 13)],
        Vec::new(),
    );
    for index in 0..2 {
        masked.evidence = masked.evidence.with_call(
            index,
            CallEvidence {
                mask: Some(CallMask {
                    alignment: MaskedAlignment::Anchoring,
                    reason: EvidenceReason::new("dephased_signal"),
                }),
                ..call('G', Some(PROFILE))
            },
        );
    }

    let evidence = run(&[differing, masked])?;

    assert_eq!(evidence.locus_differences.len(), 1);
    let locus = &evidence.locus_differences[0];
    assert_eq!(locus.position_1based, 13);
    assert_eq!(locus.support_topology.masked_reads, 1);
    assert_eq!(locus.support_topology.alternate_reads, 1);
    let masked = locus
        .observations
        .iter()
        .find(|observation| observation.state == crate::model::sample_evidence::LocusState::Masked)
        .ok_or("masked observation missing")?;
    assert_eq!(masked.base, Some('G'));
    assert_eq!(
        masked.nucleotide_contribution,
        crate::model::sample_evidence::NucleotideContribution::MaskedCall
    );
    assert_eq!(locus.nucleotide_support.contributors, 1);
    Ok(())
}

#[test]
fn differential_locus_retains_reference_support_from_overlapping_reads() -> Result<()> {
    let alternate = observation(
        "a",
        "reference",
        "config",
        Orientation::Forward,
        vec![column('G', 'A', Some(0), 72)],
        vec![snv(73, "A", "G")],
    );
    let reference = observation(
        "b",
        "reference",
        "config",
        Orientation::Reverse,
        vec![column('A', 'A', Some(0), 72)],
        Vec::new(),
    );

    let evidence = run(&[reference, alternate])?;

    assert_eq!(evidence.locus_differences.len(), 1);
    let topology = evidence.locus_differences[0].support_topology;
    assert_eq!(topology.reads, 2);
    assert_eq!(topology.forward_reads, 1);
    assert_eq!(topology.reverse_reads, 1);
    assert_eq!(topology.reference_reads, 1);
    assert_eq!(topology.alternate_reads, 1);
    assert_eq!(topology.unresolved_reads, 0);
    assert_eq!(topology.deletion_reads, 0);
    assert_eq!(topology.profile_reads, 2);
    assert_eq!(topology.profile_forward_reads, 1);
    assert_eq!(topology.profile_reverse_reads, 1);
    let observations = &evidence.locus_differences[0].observations;
    assert_eq!(observations.len(), 2);
    assert_eq!(observations[0].read_index, 0);
    assert_eq!(
        observations[0].state,
        crate::model::sample_evidence::LocusState::Alternate
    );
    assert_eq!(observations[1].read_index, 1);
    assert_eq!(
        observations[1].state,
        crate::model::sample_evidence::LocusState::Reference
    );
    assert_eq!(observations[1].call_index_0based, Some(0));
    Ok(())
}

#[test]
fn all_reference_overlap_needs_no_per_locus_records() -> Result<()> {
    let first = observation(
        "a",
        "reference",
        "config",
        Orientation::Forward,
        vec![column('A', 'A', Some(0), 10), column('C', 'C', Some(1), 11)],
        Vec::new(),
    );
    let second = observation(
        "b",
        "reference",
        "config",
        Orientation::Reverse,
        vec![column('A', 'A', Some(0), 10), column('C', 'C', Some(1), 11)],
        Vec::new(),
    );

    let evidence = run(&[first, second])?;

    assert!(evidence.locus_differences.is_empty());
    Ok(())
}

#[test]
fn records_callable_reads_that_do_not_support_a_variant_as_opposition() -> Result<()> {
    let variant = snv(73, "A", "G");
    let supporting = observation(
        "a",
        "reference",
        "config",
        Orientation::Forward,
        vec![column('G', 'A', Some(0), 72)],
        vec![variant],
    );
    let opposing = observation(
        "b",
        "reference",
        "config",
        Orientation::Reverse,
        vec![column('A', 'A', Some(0), 72)],
        Vec::new(),
    );
    let mut masked = observation(
        "c",
        "reference",
        "config",
        Orientation::Reverse,
        vec![column('A', 'A', Some(0), 72)],
        Vec::new(),
    );
    masked.alignment.callable_segments.clear();

    let evidence = run(&[supporting, opposing, masked])?;

    let opposition = &evidence.variants[0].opposition;
    assert_eq!(opposition.read_indices, [1]);
    assert_eq!((opposition.forward_reads, opposition.reverse_reads), (0, 1));
    Ok(())
}

#[test]
fn preserves_filtered_variant_observation_without_reporting_it() -> Result<()> {
    let variant = snv(73, "A", "G");
    let forward = observation(
        "a",
        "reference",
        "config",
        Orientation::Forward,
        vec![column('G', 'A', Some(0), 72)],
        vec![variant.clone()],
    );
    let mut reverse = observation(
        "b",
        "reference",
        "config",
        Orientation::Reverse,
        vec![column('G', 'A', Some(0), 72)],
        vec![variant],
    );
    reverse.variants.reported.clear();
    reverse.variants.observed[0].exclusion_reasons = vec![VariantExclusionReason::Evidence(
        EvidenceReason::new("peak_below_minimum"),
    )];

    let evidence = run(&[forward, reverse])?;

    assert_eq!(evidence.variants[0].support.len(), 2);
    assert_eq!(evidence.variants[0].support_topology.reads, 2);
    assert_eq!(evidence.variants[0].support_topology.eligible_reads, 1);
    assert_eq!(evidence.variants[0].support_topology.forward_reads, 1);
    assert_eq!(evidence.variants[0].support_topology.reverse_reads, 1);
    assert_eq!(
        evidence.variants[0].support_topology.eligible_forward_reads,
        1
    );
    assert_eq!(
        evidence.variants[0].support_topology.eligible_reverse_reads,
        0
    );
    assert!(evidence.variants[0].support[0].eligible);
    assert!(!evidence.variants[0].support[1].eligible);
    assert_eq!(
        evidence.variants[0].support[1].exclusion_reasons,
        vec![VariantExclusionReason::Evidence(EvidenceReason::new(
            "peak_below_minimum"
        ))]
    );
    assert_eq!(evidence.variants[0].support[0].calls[0].base, 'G');
    assert_eq!(
        evidence.variants[0].support[0].calls[0].call_index_0based,
        0
    );
    Ok(())
}

#[test]
fn rejects_an_aligned_call_missing_from_the_evidence() {
    let mut read = observation(
        "a",
        "reference",
        "config",
        Orientation::Forward,
        vec![column('G', 'A', Some(0), 72)],
        Vec::new(),
    );
    read.evidence = ReadEvidence::new(Vec::new(), 0..0, Vec::new())
        .unwrap_or_else(|error| panic!("valid empty evidence: {error}"));

    assert!(run(&[read]).is_err());
}

#[test]
fn rejects_duplicate_reference_coordinate_within_one_read() {
    let read = observation(
        "a",
        "reference",
        "config",
        Orientation::Forward,
        vec![column('G', 'A', Some(0), 72), column('A', 'A', Some(1), 72)],
        Vec::new(),
    );

    assert!(run(&[read]).is_err());
}

#[test]
fn rejects_duplicate_normalized_variant_identity_within_one_read() {
    let variant = snv(73, "A", "G");
    let read = observation(
        "a",
        "reference",
        "config",
        Orientation::Forward,
        vec![column('G', 'A', Some(0), 72)],
        vec![variant.clone(), variant],
    );

    assert!(run(&[read]).is_err());
}

#[test]
fn rejects_incompatible_or_duplicate_reads_even_when_renamed() {
    let first = observation(
        "a",
        "reference",
        "config",
        Orientation::Forward,
        vec![column('A', 'A', Some(0), 0)],
        Vec::new(),
    );
    let incompatible = observation(
        "b",
        "other-reference",
        "config",
        Orientation::Forward,
        vec![column('A', 'A', Some(0), 0)],
        Vec::new(),
    );
    assert!(run(&[first.clone(), incompatible]).is_err());

    let mut duplicate = observation(
        "a",
        "reference",
        "config",
        Orientation::Reverse,
        vec![column('A', 'A', Some(0), 0)],
        Vec::new(),
    );
    duplicate.input_name = "renamed-copy.ab1".into();
    assert!(run(&[first, duplicate]).is_err());
}
