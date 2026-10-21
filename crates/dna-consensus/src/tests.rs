use dna_align::{AlignmentConfig, RawAlignmentConfig, align_best};
use dna_call::RawVariantCallingConfig;
use dna_call::model::variant::VariantCallingResult;
use dna_call::read_call::RawCallConfig;
use dna_kernel::error::Result;
use dna_kernel::model::reference::ReferenceTopology;
use dna_kernel::read_evidence::{EvidenceReason, Placement, ReadEvidence};

use super::*;
use crate::model::consensus::{ConsensusRun, ConsensusSite, RunEvidence, SiteState};
use crate::observe::{Char, Trust};

mod n_reads;

const REFERENCE: &str = "TTGACGTCAGTACGATCGTACCTGAGTACGAGGACT";

fn reference() -> Reference {
    Reference {
        name: "ref".into(),
        sequence: REFERENCE.into(),
        topology: ReferenceTopology::Linear,
        sequence_sha256: String::new(),
    }
}

fn alignment() -> Result<AlignmentConfig> {
    RawAlignmentConfig {
        match_score: 3,
        mismatch_score: -5,
        ambiguous_score: 0,
        gap_open_score: -10,
        gap_extension_score: -4,
        minimum_callable_bases: 1,
        minimum_identity: 0.8,
        max_stretch: 64,
    }
    .validate()
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
            max_stretch: 64,
        },
        variant_calling: RawVariantCallingConfig {
            max_indel_length: 50,
            read_end_margin: 0,
        },
    }
    .validate()
}

/// A read of `sequence` placed on the reference, with the calls in
/// `masked` masked as dephased.
fn placed(name: &str, sequence: &str, masked: &[usize]) -> Result<CalledRead> {
    placed_with(name, sequence, masked, &[])
}

/// As [`placed`], with the calls in `unresolved` masked as unresolved.
fn placed_with(
    name: &str,
    sequence: &str,
    masked: &[usize],
    unresolved: &[usize],
) -> Result<CalledRead> {
    let mut evidence = ReadEvidence::clean(sequence);
    let masks = masked
        .iter()
        .map(|&index| {
            (
                Placement::Anchoring(EvidenceReason::new("dephased_signal")),
                index,
            )
        })
        .chain(unresolved.iter().map(|&index| {
            (
                Placement::Unresolved(EvidenceReason::new("unresolved_signal")),
                index,
            )
        }));
    for (placement, index) in masks {
        let call = evidence.calls()[index].masked(placement);
        evidence = evidence.with_call(index, call);
    }
    let alignment = align_best(&evidence, &reference(), &alignment()?)?;
    Ok(CalledRead {
        input_name: name.into(),
        input_sha256: name.into(),
        reference_sha256: String::new(),
        configuration_sha256: String::new(),
        evidence,
        alignment,
        variants: VariantCallingResult::fixture(Vec::new(), Vec::new()),
    })
}

/// The observations of a read `read_sequence` on `reference_sequence`,
/// both `ACGTACGTAC G CCCCC T AG` + `tail` in shape, with calls 19-25
/// masked as a segment that starts right after a repeat run, under a
/// read-end margin of 8.
fn exempt_observations(
    reference_sequence: &str,
    read_sequence: &str,
) -> std::result::Result<ReadObservations, Box<dyn std::error::Error>> {
    let reference = Reference {
        name: "ref".into(),
        sequence: reference_sequence.into(),
        topology: ReferenceTopology::Linear,
        sequence_sha256: String::new(),
    };
    let mut evidence = ReadEvidence::clean(read_sequence);
    for index in 19..26 {
        let mut call = evidence.calls()[index];
        call = call.masked(Placement::Unresolved(EvidenceReason::new(
            "post_homopolymer",
        )));
        evidence = evidence.with_call(index, call);
    }
    // The modality marks the C run (11-15) and its bounds as a run bounded
    // before the post-repeat phase loss at 19 (ADR-0076).
    for index in 10..=16 {
        let mut call = evidence.calls()[index];
        call.post_repeat.call = (11..=15).contains(&index);
        call.post_repeat.bound = index == 10 || index == 16;
        evidence = evidence.with_call(index, call);
    }
    let evidence = evidence.with_signal_edges(8);
    let alignment = align_best(&evidence, &reference, &alignment()?)?;
    let read = CalledRead {
        input_name: "r".into(),
        input_sha256: "r".into(),
        reference_sha256: String::new(),
        configuration_sha256: String::new(),
        evidence,
        alignment,
        variants: VariantCallingResult::fixture(Vec::new(), Vec::new()),
    };
    let settings = RawVariantCallingConfig {
        max_indel_length: 50,
        read_end_margin: 8,
    }
    .validate()?;
    Ok(observe(&read, &settings, reference.sequence.as_bytes(), 4))
}

// SRS-CONS-003
#[test]
fn an_exempt_run_trusts_its_matches_but_never_a_substituted_call()
-> std::result::Result<(), Box<dyn std::error::Error>> {
    // The read's C run 11-15 covers the reference A at 11 and is bounded
    // by G10 and T16, which match the reference.
    let tail = "AAAAAACGTACGTACGTACG";
    let observed = exempt_observations(
        &format!("ACGTACGTACGACCCCTAG{tail}"),
        &format!("ACGTACGTACGCCCCCTAG{tail}"),
    )?;
    let clean = |position: usize| observed.positions.get(&position).map(|site| site.clean);
    // The substituted call keeps `read_end`; the run's matches are clean;
    // the bound that reads the reference base is clean too, since it can
    // never make a variant.
    assert_eq!(clean(11), Some(false));
    assert!((12..=15).all(|position| clean(position) == Some(true)));
    assert_eq!(clean(16), Some(true));
    Ok(())
}

// SRS-CONS-003
#[test]
fn an_exempt_run_bound_on_another_base_keeps_its_trust()
-> std::result::Result<(), Box<dyn std::error::Error>> {
    // The read's bound at 16 reads A where the reference has T.
    let tail = "AAAAAACGTACGTACGTACG";
    let observed = exempt_observations(
        &format!("ACGTACGTACGCCCCCTAG{tail}"),
        &format!("ACGTACGTACGCCCCCAAG{tail}"),
    )?;
    let clean = |position: usize| observed.positions.get(&position).map(|site| site.clean);
    assert!((11..=15).all(|position| clean(position) == Some(true)));
    assert_eq!(clean(16), Some(false));
    Ok(())
}

// SRS-CONS-003
#[test]
fn each_placement_and_verdict_maps_to_one_trust_level()
-> std::result::Result<(), Box<dyn std::error::Error>> {
    use super::observe::{Trust, observe};
    // A read equal to the reference: call `i` observes position `i`.
    let mut read = placed_with("r", REFERENCE, &[5], &[8])?;
    let peak = EvidenceReason::new("peak_below_minimum");
    let mut vetoed = read.evidence.calls()[11];
    vetoed.identity.withheld.insert(peak)?;
    let mut base_only = read.evidence.calls()[14];
    base_only.bound = dna_kernel::read_evidence::Bound::Withheld;
    read.evidence = read.evidence.with_call(11, vetoed).with_call(14, base_only);
    let observed = observe(&read, &config()?.variant_calling, REFERENCE.as_bytes(), 4);
    let char_at = |position: usize| {
        observed
            .positions
            .get(&position)
            .and_then(|observation| observation.chars.first().copied())
    };
    let level = |position: usize| char_at(position).map(|char| char.trust);
    // Scored, trusted, nothing withheld: it shows its base.
    assert_eq!(level(2), Some(Trust::Clean));
    // Scored with a reason against its base: it counts towards its run.
    assert_eq!(level(11), Some(Trust::Unmasked));
    // Masked but anchoring: it shows a run's end.
    assert_eq!(level(5), Some(Trust::Anchoring));
    // Masked as unresolved: nothing.
    assert_eq!(level(8), Some(Trust::Nothing));
    // Vouched for by its base only: it shows its base, but neither counts
    // towards a run nor shows a run's end.
    let only = char_at(14).ok_or("no call at 14")?;
    assert!(only.clean() && !only.unmasked() && !only.resolving());
    Ok(())
}

// SRS-CONS-003
#[test]
fn each_use_of_a_call_reads_its_own_reasons() -> std::result::Result<(), Box<dyn std::error::Error>>
{
    use super::observe::observe;
    // `REFERENCE` with an A deleted after call 10 and a G inserted after
    // call 20 (0-based reference positions 11 and after 20).
    let mut sequence = REFERENCE.to_owned();
    sequence.insert(21, 'G');
    sequence.remove(11);
    let mut read = placed("r", &sequence, &[])?;
    let mixed = EvidenceReason::new("mixed_supporting_dna");
    let peak = EvidenceReason::new("peak_below_minimum");
    // Mixed signal on the inserted G withholds its base, not the
    // insertion; a weak peak beside the deletion withholds that call's
    // base, not the deletion.
    let mut inserted = read.evidence.calls()[20];
    inserted.identity.withheld.insert(mixed)?;
    let mut flank = read.evidence.calls()[10];
    flank.identity.withheld.insert(peak)?;
    flank.indel.withheld.insert(peak)?;
    read.evidence = read.evidence.with_call(20, inserted).with_call(10, flank);
    let observed = observe(&read, &config()?.variant_calling, REFERENCE.as_bytes(), 4);
    assert_eq!(
        observed.positions.get(&11).map(|deleted| deleted.clean),
        Some(true)
    );
    assert_eq!(
        observed.positions.get(&10).map(|base| base.clean),
        Some(false)
    );
    let junction = observed
        .junctions
        .values()
        .find(|junction| !junction.chars.is_empty());
    assert_eq!(junction.map(|junction| junction.clean), Some(true));
    Ok(())
}

// SRS-CONS-009
#[test]
fn an_indel_with_an_undecided_neighbour_stays_undecided() {
    let decision = |value: Option<&str>| Decision {
        value: value.map(Into::into),
        state: if value.is_some() {
            SiteState::Called
        } else {
            SiteState::Unresolved
        },
        supporting: Vec::new(),
        opposing: Vec::new(),
        uninformative: Vec::new(),
        runs: Vec::new(),
    };
    // 10 undecided, 11-12 deleted, 13 an A; then 20 a C, 21-22 deleted,
    // 23 a G.
    let mut decided = BTreeMap::from([
        (10, (10, decision(None))),
        (11, (12, decision(Some("")))),
        (13, (13, decision(Some("A")))),
        (20, (20, decision(Some("C")))),
        (21, (22, decision(Some("")))),
        (23, (23, decision(Some("G")))),
    ]);
    // 29 undecided, 30-31 a CC run lengthened to CCC that a read shows
    // in phase: its reads show the run's bounds themselves.
    let mut shown = decision(Some("CCC"));
    shown.runs = vec![run('C', 3, RunEvidence::InPhase)];
    decided.insert(29, (29, decision(None)));
    decided.insert(30, (31, shown));
    // 39 undecided, 40-41 a CC run the two-strand rule lengthened to CCC:
    // no read shows its bounds, so the bases beside it must place it, and
    // build runs this rule again after the two-strand rule.
    let mut two_strand = decision(Some("CCC"));
    two_strand.runs = vec![run('C', 3, RunEvidence::TwoStrand)];
    decided.insert(39, (39, decision(None)));
    decided.insert(40, (41, two_strand));
    unplace_unbounded_indels(&mut decided, 100);
    assert_eq!(decided[&11].1.state, SiteState::Unresolved);
    assert!(decided[&11].1.value.is_none());
    assert_eq!(decided[&21].1.value.as_deref(), Some(""));
    assert_eq!(decided[&30].1.value.as_deref(), Some("CCC"));
    assert_eq!(decided[&40].1.state, SiteState::Unresolved);
}

/// A tract decode of a run of `run` entered at call `call_start`, the
/// anchor `mark` after it: its first run `first` or one longer, and the
/// total `total` or one longer, each at shares 0.95 and 0.05.
fn decoded(
    run: char,
    mark: char,
    call_start: usize,
    first: usize,
    total: usize,
) -> dna_kernel::read_evidence::TractLengths {
    use dna_kernel::read_evidence::{LengthMarginal, Margin, StructureFit, TractLengths};
    let marginal = |length: usize| LengthMarginal {
        lengths: vec![length, length + 1],
        shares: vec![0.95, 0.05],
        errors: vec![0.01, 0.01],
        dominant: Some(length),
        next: Some(length + 1),
        z: Some(10.0),
        z_adjusted: Some(10.0),
        margin: Margin::Clear,
    };
    TractLengths {
        onset: call_start + total,
        call_start,
        extent: total,
        window: total,
        run,
        mark: Some(mark),
        candidates: vec![
            StructureFit {
                mark: None,
                residual: 0.02,
            },
            StructureFit {
                mark: Some(mark),
                residual: 0.01,
            },
        ],
        residual: 0.01,
        sigma: 0.01,
        autocorrelation: 0.0,
        first: marginal(first),
        total: marginal(total),
        cleaned: None,
    }
}

// SRS-CONS-013
#[test]
fn both_strands_decide_a_window_run_left_as_an_estimate()
-> std::result::Result<(), Box<dyn std::error::Error>> {
    // A stretch of C7 T C5 at 11-23 (1-based), the anchor T at 18, its
    // window by its own structure. Both strands decode C8 T C5: the
    // forward read enters run 1 at 11, the reverse read run 2 at 23, and
    // each derives the other run.
    let window_reference = Reference {
        name: "ref".into(),
        sequence: "TTGACGTCAGCCCCCCCTCCCCCGATCGTACCTGAGTACGAGGACT".into(),
        topology: ReferenceTopology::Linear,
        sequence_sha256: String::new(),
    };
    let read = |name: &str, sequence: &str, tract| -> Result<CalledRead> {
        let evidence = ReadEvidence::clean(sequence).with_tracts(vec![tract])?;
        let alignment = align_best(&evidence, &window_reference, &alignment()?)?;
        Ok(CalledRead {
            input_name: name.into(),
            input_sha256: name.into(),
            reference_sha256: String::new(),
            configuration_sha256: String::new(),
            evidence,
            alignment,
            variants: VariantCallingResult::fixture(Vec::new(), Vec::new()),
        })
    };
    let sequence = window_reference.sequence.as_str();
    let reverse_sequence = dna_kernel::model::nucleotide::reverse_complement(sequence);
    let forward = read("f", sequence, decoded('C', 'T', 10, 8, 14))?;
    let reverse = read(
        "r",
        &reverse_sequence,
        decoded('G', 'A', sequence.len() - 23, 5, 14),
    )?;
    let decided_with = |evidence: RunEvidence| {
        let decision = Decision {
            value: Some("CCCCCCCTCCCCC".into()),
            state: SiteState::Called,
            supporting: vec!["f".into(), "r".into()],
            opposing: Vec::new(),
            uninformative: Vec::new(),
            runs: vec![
                run('C', 7, evidence),
                run('T', 1, RunEvidence::InPhase),
                run('C', 5, RunEvidence::InPhase),
            ],
        };
        BTreeMap::from([(10, (22, decision))])
    };
    let settle =
        |evidence: RunEvidence| -> std::result::Result<Decision, Box<dyn std::error::Error>> {
            let mut decided = decided_with(evidence);
            let (reads, call, consensus) =
                ([&forward, &reverse], config()?, config::tests::config());
            two_strand(
                &mut decided,
                &reads,
                window_reference.sequence.as_bytes(),
                consensus.length_rule(),
                |first, last| {
                    read_runs_over(&reads, &window_reference, &call, &consensus, first, last)
                },
            );
            Ok(decided.remove(&10).ok_or("missing stretch")?.1)
        };
    // A phase-lost run 1 takes the strands' length; run 2, in phase,
    // keeps its own.
    let decision = settle(RunEvidence::PhaseLoss)?;
    assert_eq!(
        decision.runs,
        [
            run('C', 8, RunEvidence::TwoStrand),
            run('T', 1, RunEvidence::InPhase),
            run('C', 5, RunEvidence::InPhase)
        ]
    );
    assert_eq!(decision.value.as_deref(), Some("CCCCCCCCTCCCCC"));
    // A run whose ends two reads confirm keeps its length.
    let decision = settle(RunEvidence::ConfirmedEnd)?;
    assert_eq!(decision.runs[0], run('C', 7, RunEvidence::ConfirmedEnd));
    assert_eq!(decision.value.as_deref(), Some("CCCCCCCTCCCCC"));
    Ok(())
}

fn site_at(consensus: &Consensus, position: usize) -> Option<&ConsensusSite> {
    consensus
        .sites
        .iter()
        .find(|site| (site.first_0based..=site.last_0based).contains(&position))
}

// SRS-CONS-004, SRS-CONS-005
#[test]
fn agreeing_reads_give_their_haplotype() -> std::result::Result<(), Box<dyn std::error::Error>> {
    // A substitution at 10 (T>G) and a GG insertion after 20 in both reads.
    let haplotype = "TTGACGTCAGGACGATCGTACGGCTGAGTACGAGGACT";
    let forward = placed("f", haplotype, &[])?;
    let reverse = placed(
        "r",
        &dna_kernel::model::nucleotide::reverse_complement(haplotype),
        &[],
    )?;
    let consensus = build(
        &[&forward, &reverse],
        &reference(),
        &config()?,
        &config::tests::config(),
    );
    assert_eq!(consensus.segments.len(), 1);
    assert_eq!(consensus.segments[0].sequence, haplotype);
    // The substituted G joins the G before it: the stretch is that run.
    let substitution = site_at(&consensus, 10).ok_or("missing substitution site")?;
    assert_eq!(
        (substitution.first_0based, substitution.last_0based),
        (9, 10)
    );
    assert_eq!(substitution.call.as_deref(), Some("GG"));
    assert_eq!(substitution.state, SiteState::Called);
    assert_eq!(substitution.supporting, ["f", "r"]);
    assert_eq!(
        substitution.runs,
        [ConsensusRun {
            base: 'G',
            length: 2,
            evidence: RunEvidence::InPhase
        }]
    );
    // The insertion is decided over the whole CC run it lengthens.
    let insertion = site_at(&consensus, 20).ok_or("missing insertion site")?;
    assert_eq!((insertion.first_0based, insertion.last_0based), (20, 21));
    assert_eq!(insertion.call.as_deref(), Some("CGGC"));
    Ok(())
}

// SRS-CONS-005
#[test]
fn a_lone_read_against_a_reference_read_is_contested()
-> std::result::Result<(), Box<dyn std::error::Error>> {
    let variant = placed("v", "TTGACGTCAGGACGATCGTACCTGAGTACGAGGACT", &[])?;
    let plain = placed("p", REFERENCE, &[])?;
    let consensus = build(
        &[&variant, &plain],
        &reference(),
        &config()?,
        &config::tests::config(),
    );
    assert_eq!(consensus.segments[0].sequence, REFERENCE);
    let site = site_at(&consensus, 10).ok_or("missing contested site")?;
    assert_eq!(site.state, SiteState::Contested);
    assert_eq!(
        (site.supporting.as_slice(), site.opposing.as_slice()),
        (["p".to_owned()].as_slice(), ["v".to_owned()].as_slice())
    );
    assert_eq!(consensus.summary.contested_sites, 1);
    Ok(())
}

/// `sequence` placed with its call `index` withheld as mixed and, when
/// `reference_only`, asserted reference-only, its second base `second`
/// standing out.
fn confirming(
    name: &str,
    sequence: &str,
    index: usize,
    reference_only: bool,
) -> std::result::Result<CalledRead, Box<dyn std::error::Error>> {
    let mut read = placed(name, sequence, &[])?;
    let mut call = read.evidence.calls()[index];
    call.identity
        .withheld
        .insert(EvidenceReason::new("mixed_supporting_dna"))?;
    call.reference_only = reference_only;
    call.second = Some(dna_kernel::read_evidence::SecondBase {
        base: 'C',
        share: 0.3,
        stands_out: true,
    });
    read.evidence = read.evidence.with_call(index, call);
    Ok(read)
}

// SRS-CONS-010
#[test]
fn a_reference_only_call_confirms_the_reference_and_nothing_else()
-> std::result::Result<(), Box<dyn std::error::Error>> {
    let build_alone =
        |read: &CalledRead| -> std::result::Result<Consensus, Box<dyn std::error::Error>> {
            Ok(build(
                &[read],
                &reference(),
                &config()?,
                &config::tests::config(),
            ))
        };
    // Withheld, the call at 14 leaves its position undecided.
    let withheld = build_alone(&confirming("r", REFERENCE, 14, false)?)?;
    assert_eq!(&withheld.segments[0].sequence[13..16], "GNT");
    // Reference-only, it confirms the reference there, and the position
    // is marked as decided by that confirmation alone.
    let confirmed = build_alone(&confirming("r", REFERENCE, 14, true)?)?;
    assert_eq!(confirmed.segments[0].sequence, REFERENCE);
    assert_eq!(confirmed.reference_only, BTreeSet::from([14]));
    // A reference-only call reading another base than the reference's
    // decides nothing.
    let other = "TTGACGTCAGTACGTTCGTACCTGAGTACGAGGACT";
    let differing = build_alone(&confirming("r", other, 14, true)?)?;
    assert_eq!(&differing.segments[0].sequence[13..16], "GNT");
    assert!(differing.reference_only.is_empty());
    Ok(())
}

// SRS-CONS-010
#[test]
fn a_confirmation_beside_a_differing_site_is_withdrawn()
-> std::result::Result<(), Box<dyn std::error::Error>> {
    // The read substitutes G>T at 13, cleanly, and its call at 14 is a
    // reference-only A: the substitution's neighbour stays undecided, as
    // the read's verdict leaves it.
    let substituted = "TTGACGTCAGTACTATCGTACCTGAGTACGAGGACT";
    let read = confirming("r", substituted, 14, true)?;
    let consensus = build(&[&read], &reference(), &config()?, &config::tests::config());
    assert_eq!(&consensus.segments[0].sequence[12..16], "CTNT");
    assert!(consensus.reference_only.is_empty());
    Ok(())
}

// SRS-CONS-002
#[test]
fn reads_of_one_file_name_count_by_content() -> std::result::Result<(), Box<dyn std::error::Error>>
{
    // Three reads of one strand and one file name: two substitute G>T at
    // 13, one reads the reference. Each counts once, by its content.
    let substituted = "TTGACGTCAGTACTATCGTACCTGAGTACGAGGACT";
    let mut reads = vec![
        placed("sha-1", substituted, &[])?,
        placed("sha-2", substituted, &[])?,
        placed("sha-3", REFERENCE, &[])?,
    ];
    for read in &mut reads {
        read.input_name = "read.ab1".into();
    }
    let consensus = build(
        &reads.iter().collect::<Vec<_>>(),
        &reference(),
        &config()?,
        &config::tests::config(),
    );
    let site = consensus
        .sites
        .iter()
        .find(|site| site.first_0based <= 13 && 13 <= site.last_0based)
        .ok_or("the substitution must be a site")?;
    assert_eq!(
        site.call.as_deref().map(|call| call.contains('T')),
        Some(true)
    );
    assert_eq!(site.supporting, ["sha-1", "sha-2"]);
    assert_eq!(site.opposing, ["sha-3"]);
    Ok(())
}

/// `sequence` placed with the calls at `edges` at a signal edge and
/// reference-only: confirmations the modality's edge does not stop.
fn at_edges(
    name: &str,
    sequence: &str,
    edges: &[usize],
) -> std::result::Result<CalledRead, Box<dyn std::error::Error>> {
    let mut read = placed(name, sequence, &[])?;
    for &index in edges {
        let mut call = read.evidence.calls()[index];
        call.identity.edge = true;
        call.indel.edge = true;
        call.reference_only = true;
        read.evidence = read.evidence.with_call(index, call);
    }
    Ok(read)
}

// SRS-CONS-010
#[test]
fn a_reference_only_call_at_an_edge_confirms_clear_of_the_read_ends()
-> std::result::Result<(), Box<dyn std::error::Error>> {
    // At an edge, the call at 14 would leave its position undecided; as
    // a reference-only call it confirms the reference there.
    let read = at_edges("r", REFERENCE, &[14])?;
    let consensus = build(&[&read], &reference(), &config()?, &config::tests::config());
    assert!(consensus.reference_only.contains(&14));
    assert_eq!(consensus.segments[0].sequence, REFERENCE);
    Ok(())
}

// SRS-CONS-010
#[test]
fn a_vouched_call_within_the_end_margin_confirms_the_reference()
-> std::result::Result<(), Box<dyn std::error::Error>> {
    use super::observe::observe;
    // A read equal to the reference, its ends unvouched: calls 0-2 lie
    // within a margin of three calls of its start, and call `i` observes
    // position `i`.
    let mut read = placed("r", REFERENCE, &[])?;
    let mut edge = read.evidence.calls()[0];
    edge.identity.edge = true;
    let mut second = read.evidence.calls()[1];
    second.second = Some(dna_kernel::read_evidence::SecondBase {
        base: 'A',
        share: 0.3,
        stands_out: true,
    });
    read.evidence = read.evidence.with_call(0, edge).with_call(1, second);
    let settings = RawVariantCallingConfig {
        max_indel_length: 50,
        read_end_margin: 3,
    }
    .validate()?;
    let observed = observe(&read, &settings, REFERENCE.as_bytes(), 4);
    let site = |position: usize| {
        observed
            .positions
            .get(&position)
            .map(|site| (site.clean, site.confirms))
    };
    // The modality's edge verdict, or a second base standing out, leaves
    // the reference unconfirmed; the call at 2, vouched for its base,
    // confirms it. Beyond the margin, calls are clean on their own.
    assert_eq!(site(0), Some((false, false)));
    assert_eq!(site(1), Some((false, false)));
    assert_eq!(site(2), Some((true, true)));
    assert_eq!(site(3), Some((true, false)));
    Ok(())
}

// SRS-CONS-010
#[test]
fn a_confirmation_beside_any_observed_difference_is_withdrawn()
-> std::result::Result<(), Box<dyn std::error::Error>> {
    // A second read reads T at 16 (C>T), uncleanly: a shifted read can
    // show that, so the confirmation at 15 beside it is withdrawn, while
    // the one at 14 stands.
    let confirming_read = at_edges("r", REFERENCE, &[14, 15])?;
    let mut shifted = placed("s", "TTGACGTCAGTACGATTGTACCTGAGTACGAGGACT", &[])?;
    for index in [14, 15, 16] {
        let mut call = shifted.evidence.calls()[index];
        call.identity.edge = true;
        shifted.evidence = shifted.evidence.with_call(index, call);
    }
    let consensus = build(
        &[&confirming_read, &shifted],
        &reference(),
        &config()?,
        &config::tests::config(),
    );
    assert!(consensus.reference_only.contains(&14));
    assert!(!consensus.reference_only.contains(&15));
    Ok(())
}

// SRS-CONS-010
#[test]
fn a_stretch_is_decided_without_confirmations() {
    // T12C merges the tract into one stretch (7-16). A second read reads
    // the tract uncleanly but confirms the reference at 9: the stretch's
    // decision is the same with and without that confirmation. This
    // guards the outcome only: the stretch here does not fall back to its
    // position-by-position decision, where the exclusion acts.
    let merged = "GTCAAAACCCCCCCCCCATG";
    let observed = |confirm: bool| {
        let clean = rows("f", 0, TRACT, merged, &"c".repeat(TRACT.len()));
        let mut other = rows("g", 0, TRACT, TRACT, &"v".repeat(TRACT.len()));
        if confirm && let Some(observation) = other.positions.get_mut(&9) {
            observation.confirms = true;
            observation.clean = true;
        }
        vec![clean, other]
    };
    let decided = |confirm: bool| {
        decide_positions(&observed(confirm), TRACT.as_bytes(), 4)
            .get(&7)
            .map(|(last, decision)| (*last, decision.value.clone(), decision.state))
    };
    assert!(decided(false).is_some());
    assert_eq!(decided(true), decided(false));
}

// SRS-CONS-010
#[test]
fn a_confirmation_never_closes_an_unobserved_junction_in_a_long_run() {
    // The C run at 7-11 of the tract, five long. A confirmation at 9 with
    // no read observing the junction after 8 cleanly is withdrawn where
    // such a run is long enough to be a stretch, and stands otherwise.
    let confirmations = |seeded_run: usize| {
        let mut read = rows("r", 0, TRACT, TRACT, &"c".repeat(TRACT.len()));
        if let Some(observation) = read.positions.get_mut(&9) {
            observation.confirms = true;
        }
        if let Some(junction) = read.junctions.get_mut(&8) {
            junction.clean = false;
        }
        let mut observed = vec![read];
        withdraw_unobserved_run_junctions(&mut observed, TRACT.as_bytes(), seeded_run);
        observed[0].positions[&9].confirms
    };
    assert!(!confirmations(4));
    assert!(confirmations(6));
}

// SRS-CONS-010
#[test]
fn a_second_sequence_any_read_shows_withdraws_a_confirmation() {
    // Two reads at position 14: the first confirms the reference with a
    // second C; the second does not observe the position cleanly.
    let observations = |first_forward: bool, second: (bool, Option<char>, bool)| {
        let mut first = rows("a", 0, REFERENCE, REFERENCE, &"c".repeat(REFERENCE.len()));
        first.forward = first_forward;
        if let Some(observation) = first.positions.get_mut(&14) {
            observation.confirms = true;
        }
        first.seconds.insert(14, 'C');
        let mut other = rows("b", 0, REFERENCE, REFERENCE, &"c".repeat(REFERENCE.len()));
        other.forward = second.0;
        if let Some(base) = second.1 {
            other.seconds.insert(14, base);
        }
        if second.2 {
            other.unexplained.insert(14);
        }
        let mut observed = vec![first, other];
        hide_confirmations(&mut observed);
        observed[0].positions[&14].confirms
    };
    // Nothing on the other read: the confirmation stands.
    assert!(observations(true, (false, None, false)));
    // The other strand shows the same second base standing out.
    assert!(!observations(true, (false, Some('C'), false)));
    // The same strand doing so is not a second sequence by itself.
    assert!(observations(true, (true, Some('C'), false)));
    // A read of the confirming strand with a second peak its modality
    // does not explain.
    assert!(!observations(true, (true, None, true)));
    assert!(!observations(true, (true, Some('A'), true)));
}

// SRS-CONS-010
#[test]
fn an_unexplained_second_on_the_other_strand_alone_keeps_a_confirmation() {
    // Read `a` (forward) confirms the reference at 14 with a second C,
    // read `b` (reverse) shows an unexplained second peak there, and read
    // `c` (forward) optionally shows a second base standing out.
    let observations = |second: Option<char>, third: Option<char>| {
        let read = |name: &str, forward: bool| {
            let mut read = rows(name, 0, REFERENCE, REFERENCE, &"c".repeat(REFERENCE.len()));
            read.forward = forward;
            read
        };
        let mut first = read("a", true);
        if let Some(observation) = first.positions.get_mut(&14) {
            observation.confirms = true;
        }
        first.seconds.insert(14, 'C');
        let mut other = read("b", false);
        other.unexplained.insert(14);
        if let Some(base) = second {
            other.seconds.insert(14, base);
        }
        let mut same = read("c", true);
        if let Some(base) = third {
            same.seconds.insert(14, base);
        }
        let mut observed = vec![first, other, same];
        hide_confirmations(&mut observed);
        observed[0].positions[&14].confirms
    };
    // A second base the confirming strand does not show: one strand
    // alone, so the confirmation stands.
    assert!(observations(Some('A'), None));
    assert!(observations(Some('A'), Some('G')));
    // The confirming strand shows the same second base standing out.
    assert!(!observations(Some('A'), Some('A')));
    // The confirming call's own second base: the other-strand guard.
    assert!(!observations(Some('C'), None));
    // No second base standing out to compare: still withdrawn.
    assert!(!observations(None, None));
}

// SRS-CONS-012
#[test]
fn positions_only_unreliable_spans_observe_are_unreliable_only()
-> std::result::Result<(), Box<dyn std::error::Error>> {
    let reason = EvidenceReason::new("poor_start_signal");
    let mut poor = placed("p", REFERENCE, &[])?;
    poor.evidence =
        poor.evidence
            .with_unreliable(vec![dna_kernel::read_evidence::UnreliableSpan {
                calls: 0..5,
                reason,
            }])?;
    // Alone, the read observes its span's positions only unreliably.
    let alone = build(&[&poor], &reference(), &config()?, &config::tests::config());
    assert_eq!(
        alone.unreliable_only.into_iter().collect::<Vec<_>>(),
        (0..5)
            .map(|position| (position, reason))
            .collect::<Vec<_>>()
    );
    // A read that observes 3 and 4 cleanly leaves 0..=2.
    let late = placed("l", &REFERENCE[3..], &[])?;
    let both = build(
        &[&poor, &late],
        &reference(),
        &config()?,
        &config::tests::config(),
    );
    assert_eq!(
        both.unreliable_only.keys().copied().collect::<Vec<_>>(),
        [0, 1, 2]
    );
    Ok(())
}

// SRS-CONS-003, SRS-CONS-006
#[test]
fn masked_calls_never_decide() -> std::result::Result<(), Box<dyn std::error::Error>> {
    // The variant read's call at 10 is masked: the reference read decides.
    let variant = placed("v", "TTGACGTCAGGACGATCGTACCTGAGTACGAGGACT", &[10])?;
    let plain = placed("p", REFERENCE, &[])?;
    let consensus = build(
        &[&variant, &plain],
        &reference(),
        &config()?,
        &config::tests::config(),
    );
    assert_eq!(consensus.segments[0].sequence, REFERENCE);
    assert!(site_at(&consensus, 10).is_none());
    // Alone, a masked call leaves the position unresolved.
    let alone = build(
        &[&variant],
        &reference(),
        &config()?,
        &config::tests::config(),
    );
    assert_eq!(&alone.segments[0].sequence[8..12], "AGNA");
    assert_eq!(alone.summary.unresolved_positions, 1);
    Ok(())
}

// SRS-CONS-003
#[test]
fn no_masked_call_places_a_substitution_or_an_indel()
-> std::result::Result<(), Box<dyn std::error::Error>> {
    // A substitution (A14T), an insertion (a G after 14), and a deletion
    // (of 14), each read only through masked calls, anchoring or
    // unresolved: a masked call shows at most a run's end (ADR-0071), so
    // none is decided, beside a reference read or alone.
    let substituted = "TTGACGTCAGTACGTTCGTACCTGAGTACGAGGACT";
    let inserted = "TTGACGTCAGTACGAGTCGTACCTGAGTACGAGGACT";
    let deleted = "TTGACGTCAGTACGTCGTACCTGAGTACGAGGACT";
    let cases: [(&str, &[usize]); 3] = [
        (substituted, &[14]),
        (inserted, &[15]),
        (deleted, &[13, 14]),
    ];
    let plain = placed("p", REFERENCE, &[])?;
    for (sequence, masked) in cases {
        for read in [
            placed("v", sequence, masked)?,
            placed_with("v", sequence, &[], masked)?,
        ] {
            let both = build(
                &[&read, &plain],
                &reference(),
                &config()?,
                &config::tests::config(),
            );
            assert_eq!(both.segments[0].sequence, REFERENCE);
            let alone = build(&[&read], &reference(), &config()?, &config::tests::config());
            // Every segment reads its reference interval, or `N`.
            for segment in &alone.segments {
                let expected = &REFERENCE[segment.start_0based..=segment.end_0based];
                assert_eq!(segment.sequence.len(), expected.len(), "{sequence}");
                assert!(
                    segment
                        .sequence
                        .chars()
                        .zip(expected.chars())
                        .all(|(call, reference)| call == reference || call == 'N'),
                    "{sequence} gave {}",
                    segment.sequence
                );
            }
        }
    }
    Ok(())
}

// SRS-CONS-003
#[test]
fn an_insertion_into_a_run_the_read_does_not_bound_never_decides()
-> std::result::Result<(), Box<dyn std::error::Error>> {
    // A C inserted into the CC run at 20-21; the A before the run is call 19.
    let lengthened = "TTGACGTCAGTACGATCGTACCCTGAGTACGAGGACT";
    let bounded = build(
        &[&placed("b", lengthened, &[])?],
        &reference(),
        &config()?,
        &config::tests::config(),
    );
    assert_eq!(bounded.segments[0].sequence, lengthened);
    // Unresolved, the A no longer shows where the run starts, so the read
    // cannot tell the run's length and the reference length stays.
    let unbounded = build(
        &[&placed_with("u", lengthened, &[], &[19])?],
        &reference(),
        &config()?,
        &config::tests::config(),
    );
    let mut expected = REFERENCE.to_owned();
    expected.replace_range(19..20, "N");
    assert_eq!(unbounded.segments[0].sequence, expected);
    assert!(site_at(&unbounded, 20).is_none());
    Ok(())
}
/// A read observed over `reference` from 0-based `first`, given as two
/// alignment rows (`-` for an inserted or deleted base) and one kind per
/// column: `c` clean, `m` masked but anchoring, `v` unmasked but vetoed;
/// an `N` in the read row is an unresolved call.
fn rows(name: &str, first: usize, reference: &str, read: &str, kinds: &str) -> ReadObservations {
    use super::observe::Observation;
    let char_of = |base: char, kind: char, position: Option<usize>| {
        let canonical = base != 'N';
        Char {
            base,
            trust: match kind {
                _ if !canonical => Trust::Nothing,
                'c' => Trust::Clean,
                'v' => Trust::Unmasked,
                _ => Trust::Anchoring,
            },
            position,
            beside_unresolved: false,
            base_only: false,
        }
    };
    let mut positions = BTreeMap::new();
    let mut junctions = BTreeMap::new();
    let mut position = first;
    let mut previous: Option<usize> = None;
    let mut inserted: Vec<Char> = Vec::new();
    for ((reference, read), kind) in reference.chars().zip(read.chars()).zip(kinds.chars()) {
        if reference == '-' {
            inserted.push(char_of(read, kind, None));
            continue;
        }
        if let Some(previous) = previous {
            junctions.insert(
                previous,
                Observation {
                    value: SiteValue::Insertion(inserted.iter().map(|char| char.base).collect()),
                    confirms: false,
                    clean: inserted.iter().all(|char| char.clean()),
                    chars: std::mem::take(&mut inserted),
                },
            );
        }
        let observation = if read == '-' {
            Observation {
                value: SiteValue::Deletion,
                confirms: false,
                clean: true,
                chars: Vec::new(),
            }
        } else {
            let char = char_of(read, kind, Some(position));
            Observation {
                value: SiteValue::Base(read),
                confirms: false,
                clean: char.clean(),
                chars: vec![char],
            }
        };
        positions.insert(position, observation);
        previous = Some(position);
        position += 1;
    }
    ReadObservations {
        forward: true,
        seconds: BTreeMap::new(),
        unexplained: BTreeSet::new(),
        edges: BTreeMap::new(),
        id: name.into(),
        positions,
        junctions,
    }
}

/// A forward read of the reference and a reverse read that differs at 4,
/// both `edge` calls from their ends everywhere.
fn dissent(reverse: &str, kinds: &str, edge: usize) -> Vec<ReadObservations> {
    const SHORT: &str = "ACGTACGTAC";
    let mut forward = rows("f", 0, SHORT, SHORT, "cccccccccc");
    let mut opposing = rows("r", 0, SHORT, reverse, kinds);
    opposing.forward = false;
    for read in [&mut forward, &mut opposing] {
        read.edges = (0..SHORT.len()).map(|position| (position, edge)).collect();
    }
    vec![forward, opposing]
}

/// The decision of the site or stretch holding `position`.
fn decided_at(observed: &[ReadObservations], position: usize) -> Decision {
    let decided = decide_positions(observed, b"ACGTACGTAC", 4);
    decided
        .range(..=position)
        .next_back()
        .map_or_else(|| decide(&[], "A"), |(_, (_, decision))| decision.clone())
}

/// The decision of the stretch holding `position` once one-strand
/// insertions are contested (SRS-CONS-015).
fn contested_at(observed: &[ReadObservations], bases: &[u8], position: usize) -> Decision {
    let mut decided = decide_positions(observed, bases, 4);
    insertion::contest_one_strand(&mut decided, observed, bases, 4);
    decided
        .range(..=position)
        .next_back()
        .map_or_else(|| decide(&[], "A"), |(_, (_, decision))| decision.clone())
}

// SRS-CONS-015
#[test]
fn a_short_insertion_the_other_strand_reads_through_is_contested() {
    const SHORT: &str = "ACGATTCAGC";
    let reads = |aligned: &str, reverse: &str, kinds: &str| {
        let forward = rows("f", 0, "ACGATT-CAGC", "ACGATTTCAGC", "ccccccccccc");
        let mut other = rows("r", 0, aligned, reverse, kinds);
        other.forward = false;
        vec![forward, other]
    };
    // The reverse read reads TT C masked but anchoring: not clean, and
    // showing no run's length, but it reads the site through without the
    // inserted T.
    let through = reads(SHORT, SHORT, "cccaaaaccc");
    let before = decided_at_with(&through, SHORT.as_bytes(), 5);
    assert_eq!(
        (before.value.as_deref(), before.state),
        (Some("TTTC"), SiteState::Called)
    );
    let after = contested_at(&through, SHORT.as_bytes(), 5);
    assert_eq!(
        (after.value.as_deref(), after.state),
        (Some("TTC"), SiteState::Contested)
    );
    assert_eq!(
        (after.supporting, after.opposing),
        (vec!["r".to_owned()], vec!["f".to_owned()])
    );
    assert!(after.runs.is_empty());
    // A reverse read that shows the insertion, even masked, confirms it.
    let shows = reads("ACGATT-CAGC", "ACGATTTCAGC", "cccaaaaaccc");
    assert_eq!(
        contested_at(&shows, SHORT.as_bytes(), 5).value.as_deref(),
        Some("TTTC")
    );
    // An unresolved call there gives the reverse strand no view.
    let unresolved = reads(SHORT, "ACGATNCAGC", "cccaaaaccc");
    assert_eq!(
        contested_at(&unresolved, SHORT.as_bytes(), 5)
            .value
            .as_deref(),
        Some("TTTC")
    );
    // A read of the same strand reading it through does not contest it.
    let mut same_strand = reads(SHORT, SHORT, "cccaaaaccc");
    same_strand[1].forward = true;
    assert_eq!(
        contested_at(&same_strand, SHORT.as_bytes(), 5)
            .value
            .as_deref(),
        Some("TTTC")
    );
}

// SRS-CONS-015
#[test]
fn a_read_views_an_insertion_site_by_what_it_shows() {
    const SITE: &str = "ACGATTCAGC";
    let view = |first: usize, aligned: &str, read: &str, kinds: &str| {
        let observed = rows("r", first, aligned, read, kinds);
        insertion::view(&observed, 4, 6, "TTTC", SITE.len())
    };
    assert_eq!(
        view(0, "ACGATT-CAGC", "ACGATTTCAGC", "aaaaaaaaaaa"),
        InsertionView::Shows
    );
    assert_eq!(
        view(0, SITE, SITE, "cccaaaaccc"),
        InsertionView::ReadsThrough
    );
    // An unresolved call, or one outside the informative interval, shows
    // no base.
    assert_eq!(
        view(0, SITE, "ACGATNCAGC", "cccaaaaccc"),
        InsertionView::Unresolved
    );
    // Reaching the site but not the position after it.
    assert_eq!(
        view(0, "ACGATT", "ACGATT", "cccccc"),
        InsertionView::ReadEnd
    );
    // Ending before the site.
    assert_eq!(view(0, "ACG", "ACG", "ccc"), InsertionView::NotCovered);
}

// SRS-CONS-015
#[test]
fn an_insertion_into_a_long_run_is_left_to_the_run_rules() {
    const LONG: &str = "ACGTTTTTCAGC";
    let forward = rows("f", 0, "ACGTTTTT-CAGC", "ACGTTTTTTCAGC", "ccccccccccccc");
    let mut other = rows("r", 0, LONG, LONG, "cccaaaaaaccc");
    other.forward = false;
    let observed = vec![forward, other];
    assert_eq!(
        contested_at(&observed, LONG.as_bytes(), 5).value.as_deref(),
        Some("TTTTTTC")
    );
}

/// The decision of the site or stretch holding `position` over `bases`.
fn decided_at_with(observed: &[ReadObservations], bases: &[u8], position: usize) -> Decision {
    decide_positions(observed, bases, 4)
        .range(..=position)
        .next_back()
        .map_or_else(|| decide(&[], "A"), |(_, (_, decision))| decision.clone())
}

// SRS-CONS-011
#[test]
fn an_explained_second_on_the_other_strand_is_set_aside() {
    let mut observed = dissent("ACGTTCGTAC", "cccccccccc", 50);
    assert_eq!(decided_at(&observed, 4).opposing, ["r"]);
    // Its second base stands out, and the modality explains its peak.
    observed[1].seconds.insert(4, 'A');
    let decided = decide_positions(&observed, b"ACGTACGTAC", 4);
    let mut unexplained = observed.clone();
    unexplained[1].unexplained.insert(4);
    assert!(!set_aside_dissent(
        &mut unexplained,
        &decided,
        b"ACGTACGTAC",
        12,
        40
    ));
    let mut same_strand = observed.clone();
    same_strand[1].forward = true;
    assert!(!set_aside_dissent(
        &mut same_strand,
        &decided,
        b"ACGTACGTAC",
        12,
        40
    ));
    let mut near_end = observed.clone();
    near_end[0].edges.insert(4, 11);
    assert!(!set_aside_dissent(
        &mut near_end,
        &decided,
        b"ACGTACGTAC",
        12,
        40
    ));
    assert!(set_aside_dissent(
        &mut observed,
        &decided,
        b"ACGTACGTAC",
        12,
        40
    ));
    // Set aside, the read shows nothing at the stretch its substitution
    // made, and the position is decided alone.
    let after = decided_at(&observed, 4);
    assert_eq!(after.value.as_deref(), Some("A"));
    assert_eq!(after.state, SiteState::Called);
    assert!(after.opposing.is_empty());
}

// SRS-CONS-011
#[test]
fn a_run_length_dissent_from_a_called_variant_is_set_aside() {
    const RUN: &str = "ACATGGGGTACA";
    let reads = |opposing: &str| {
        let forward = rows("f", 0, "ACATGGGG-TACA", "ACATGGGGGTACA", "ccccccccccccc");
        let mut reverse = rows("r", 0, "ACATGGGG-TACA", "ACATGGGGGTACA", "ccccccccccccc");
        reverse.forward = false;
        let mut other = rows("o", 0, RUN, opposing, "cccccccccccc");
        // Its second base stands out, and the modality explains its peak.
        other.seconds.insert(8, 'G');
        let mut observed = vec![forward, reverse, other];
        for read in &mut observed {
            read.edges = (0..RUN.len()).map(|position| (position, 50)).collect();
        }
        observed
    };
    let set_aside = |observed: &mut Vec<ReadObservations>| {
        let decided = decide_positions(observed, RUN.as_bytes(), 4);
        set_aside_dissent(observed, &decided, RUN.as_bytes(), 12, 40)
    };
    let decided_at = |observed: &[ReadObservations]| {
        decide_positions(observed, RUN.as_bytes(), 4)
            .range(..=5)
            .next_back()
            .map(|(_, (_, decision))| decision.clone())
    };
    let mut length = reads(RUN);
    let before = decided_at(&length).map(|decision| (decision.value, decision.opposing));
    assert_eq!(
        before,
        Some((Some("GGGGGT".to_owned()), vec!["o".to_owned()]))
    );
    assert!(set_aside(&mut length));
    let after = decided_at(&length).map(|decision| (decision.value, decision.opposing));
    assert_eq!(after, Some((Some("GGGGGT".to_owned()), Vec::new())));
    // Without a second base standing out, the dissent stands.
    let mut plain = reads(RUN);
    plain[2].seconds.clear();
    assert!(!set_aside(&mut plain));
    // A read showing other runs differs in more than a length.
    let mut substituted = reads("ACATGGAGTACA");
    assert!(!set_aside(&mut substituted));
}

// SRS-CONS-011
#[test]
fn an_indel_near_its_read_end_is_set_aside() {
    let near = dissent("ACGT-CGTAC", "cccccccccc", 50);
    let decided = decide_positions(&near, b"ACGTACGTAC", 4);
    assert_eq!(decided_at(&near, 4).opposing, ["r"]);
    // Far from its ends, the deletion stands.
    assert!(!set_aside_dissent(
        &mut near.clone(),
        &decided,
        b"ACGTACGTAC",
        12,
        40
    ));
    let mut close = near.clone();
    close[1].edges.insert(4, 5);
    // The supporting read must itself lie that far from its ends.
    let mut short = close.clone();
    short[0].edges.insert(4, 39);
    assert!(!set_aside_dissent(
        &mut short,
        &decided,
        b"ACGTACGTAC",
        12,
        40
    ));
    assert!(set_aside_dissent(
        &mut close,
        &decided,
        b"ACGTACGTAC",
        12,
        40
    ));
    assert!(decided_at(&close, 4).opposing.is_empty());
    // A substitution near the end is no indel.
    let mut substituted = dissent("ACGTTCGTAC", "cccccccccc", 50);
    substituted[1].edges.insert(4, 5);
    let decided = decide_positions(&substituted, b"ACGTACGTAC", 4);
    assert!(!set_aside_dissent(
        &mut substituted,
        &decided,
        b"ACGTACGTAC",
        12,
        40
    ));
}

fn run(base: char, length: usize, evidence: RunEvidence) -> ConsensusRun {
    ConsensusRun {
        base,
        length,
        evidence,
    }
}

/// An A-tract, then the C-tract `CCCCCTCCCC` (0-based 7-16).
const TRACT: &str = "GTCAAAACCCCCTCCCCATG";

// SRS-CONS-004, SRS-CONS-008
#[test]
fn a_run_two_reads_enter_from_its_ends_takes_their_agreed_count() {
    // T12C merges the tract. The forward read enters it in phase and
    // loses phase at its end, on an unresolved call; the reverse read the
    // other way round. Both read every call of the run unmasked and count
    // ten.
    let forward = rows(
        "f",
        0,
        TRACT,
        "GTCAAAACCCCCCCCCCNTG",
        "cccccccccccccccccmmm",
    );
    let reverse = rows(
        "r",
        0,
        TRACT,
        "GTCAAANCCCCCCCCCCATG",
        "mmmmmmmccccccccccccc",
    );
    let decision = decide_stretch(&[forward.clone(), reverse.clone()], 7, 16, TRACT.as_bytes());
    assert_eq!(decision.value.as_deref(), Some("CCCCCCCCCC"));
    assert_eq!(decision.runs, [run('C', 10, RunEvidence::ConfirmedEnd)]);
    assert_eq!(decision.supporting, ["f", "r"]);
    // When the reverse read loses phase one call earlier, the counts
    // differ: the run takes the reference frame between its observed ends.
    let early = rows(
        "r",
        0,
        TRACT,
        "GTCAAANCCCCCCCCCCATG",
        "mmmmmmmmcccccccccccc",
    );
    let decision = decide_stretch(&[forward, early], 7, 16, TRACT.as_bytes());
    assert_eq!(decision.value.as_deref(), Some("CCCCCCCCCC"));
    assert_eq!(decision.runs, [run('C', 10, RunEvidence::ReferenceFrame)]);
    // When the forward read's exit call reads the run's base, masked, the
    // run goes on past its count: the count bounds the run from below.
    let going_on = rows(
        "f",
        0,
        TRACT,
        "GTCAAAACCCCCCCCCCCTG",
        "cccccccccccccccccmmm",
    );
    let decision = decide_stretch(&[going_on, reverse], 7, 16, TRACT.as_bytes());
    assert_eq!(decision.runs, [run('C', 10, RunEvidence::ReferenceFrame)]);
}

// SRS-CONS-004
#[test]
fn a_clean_end_holding_every_run_does_not_stand_for_a_read_with_more_runs() {
    // The tract CCCCC T CCCC (7-16). One read shows its run structure
    // with its C run masked. Another reads C at 7 and T at 11 unclean,
    // a C inserted after 11, then T and CCCC cleanly: its clean end, C T
    // CCCC, has the composition's runs, but the read reads two more runs
    // before them, so it shows neither where the first C run starts nor
    // its length.
    let structured = rows("s", 0, TRACT, TRACT, "cccccccmmmmmcccccccc");
    let longer = rows(
        "l",
        0,
        "GTCAAAACCCCC-TCCCCATG",
        "GTCAAAACCCCTCTCCCCATG",
        "cccccccvcccvccccccccc",
    );
    let decision = decide_stretch(&[structured, longer], 7, 16, TRACT.as_bytes());
    assert_ne!(decision.value.as_deref(), Some("CTCCCC"));
    assert!(
        decision.runs.first().is_none_or(|run| run.length != 1),
        "{decision:?}"
    );
}

// SRS-CONS-004
#[test]
fn a_read_reading_only_the_exit_cleanly_opposes_a_run_that_swallows_it() {
    // The forward read reads the whole stretch 3-17 as `AAAA` and eleven
    // C, its exit beyond 17 anchored: the tract swallows the A at 17. The
    // reverse read reads only the exit cleanly, `TCCCCA` at 12-17, in
    // phase: it contradicts the composition's last runs.
    let forward = rows(
        "f",
        0,
        TRACT,
        "GTCAAAACCCCCCCCCCCTG",
        "ccccccccccccccccccmm",
    );
    let reverse = rows("r", 0, TRACT, TRACT, "mmmmmmmmmmmmcccccccc");
    let decision = decide_stretch(&[forward.clone(), reverse], 3, 17, TRACT.as_bytes());
    assert_ne!(decision.value.as_deref(), Some("AAAACCCCCCCCCCC"));
    assert_ne!(decision.state, SiteState::Called);
    // Alone, the forward read decides the stretch.
    let alone = decide_stretch(&[forward], 3, 17, TRACT.as_bytes());
    assert_eq!(alone.value.as_deref(), Some("AAAACCCCCCCCCCC"));
}

// SRS-CONS-004, SRS-CONS-008
#[test]
fn a_read_losing_phase_in_a_run_does_not_decide_what_follows_it() {
    // The forward read reads an eleventh C in phase over the A at 17 before
    // it loses phase; the reverse read reads that A in phase.
    let forward = rows(
        "f",
        0,
        TRACT,
        "GTCAAAACCCCCCCCCCCCG",
        "ccccccccccccccccccmm",
    );
    let reverse = rows(
        "r",
        0,
        TRACT,
        "GTCAAACCCCCCCCCCCATG",
        "mmmmmmmccccccccccccc",
    );
    let decision = decide_stretch(&[forward, reverse], 7, 17, TRACT.as_bytes());
    assert_eq!(decision.value.as_deref(), Some("CCCCCCCCCCA"));
    assert_eq!(
        decision.runs,
        [
            run('C', 10, RunEvidence::ReferenceFrame),
            run('A', 1, RunEvidence::InPhase)
        ]
    );
}

/// A C7 run, T, and a C5 run (0-based 6-18), as in a poly-C window.
const WINDOW: &str = "GTCAAACCCCCCCTCCCCCGCTT";
const WINDOW_ROW: &str = "GTCAAACCCCCCC--TCCCCC-GCTT";

// SRS-CONS-004, SRS-CONS-008
#[test]
fn a_run_read_until_phase_loss_takes_its_in_phase_length() {
    // A C9 T C6 haplotype; the only read loses phase where the C9 run
    // starts, so no read locates that end.
    let reverse = rows(
        "r",
        0,
        WINDOW_ROW,
        "GTCAACCCCCCCCCCTCCCCCCGCTT",
        "mmmmmmcccccccccccccccccccc",
    );
    let decision = decide_stretch(&[reverse], 6, 18, WINDOW.as_bytes());
    assert_eq!(decision.value.as_deref(), Some("CCCCCCCCCTCCCCCC"));
    assert_eq!(
        decision.runs,
        [
            run('C', 9, RunEvidence::PhaseLoss),
            run('T', 1, RunEvidence::InPhase),
            run('C', 6, RunEvidence::InPhase)
        ]
    );
}

// SRS-CONS-004, SRS-CONS-008
#[test]
fn a_lone_unresolved_call_at_a_run_end_shows_an_anchored_end() {
    // The same read, its phase lost through one unresolved call before a
    // masked A: the run's length is known only as the dominant one.
    let reverse = rows(
        "r",
        0,
        WINDOW_ROW,
        "GTCAANCCCCCCCCCTCCCCCCGCTT",
        "mmmmmmcccccccccccccccccccc",
    );
    let decision = decide_stretch(&[reverse], 6, 18, WINDOW.as_bytes());
    assert_eq!(decision.value.as_deref(), Some("CCCCCCCCCTCCCCCC"));
    assert_eq!(decision.runs[0], run('C', 9, RunEvidence::AnchoredEnd));
}

// SRS-CONS-004
#[test]
fn a_disputed_run_no_read_bounds_seeds_a_stretch() {
    // A C inserted into the C7 run (0-based 6-12) through masked calls: the
    // read shows neither run end cleanly, yet it differs on the run.
    let inserting = || {
        rows(
            "i",
            0,
            "GTCAAACCCCCCC-TCCCCCGCTT",
            "GTCAAACCCCCCCCTCCCCCGCTT",
            "ccccccmmmmmmmmmccccccccc",
        )
    };
    assert_eq!(stretches(&[inserting()], WINDOW.as_bytes(), 4), [(6, 12)]);
    // A read that shows both ends of the run cleanly bounds it.
    let bounding = rows("b", 0, WINDOW, WINDOW, &"c".repeat(WINDOW.len()));
    assert!(stretches(&[inserting(), bounding], WINDOW.as_bytes(), 4).is_empty());
}

// SRS-CONS-007
#[test]
fn an_estimated_length_equal_to_the_reference_keeps_its_site() {
    let estimated = |evidence: RunEvidence| {
        let decision = Decision {
            value: Some(TRACT[7..=16].into()),
            state: SiteState::Called,
            supporting: vec!["f".into()],
            opposing: Vec::new(),
            uninformative: Vec::new(),
            runs: vec![
                run('C', 5, evidence),
                run('T', 1, RunEvidence::InPhase),
                run('C', 4, RunEvidence::InPhase),
            ],
        };
        let decided = BTreeMap::from([(7, (16, decision))]);
        let mut consensus = Consensus {
            segments: Vec::new(),
            sites: Vec::new(),
            summary: ConsensusSummary::default(),
            reference_only: BTreeSet::new(),
            unreliable_only: BTreeMap::new(),
        };
        assemble(&[7], &decided, TRACT.as_bytes(), 1, &mut consensus);
        consensus.sites.len()
    };
    assert_eq!(estimated(RunEvidence::ReferenceFrame), 1);
    assert_eq!(estimated(RunEvidence::InPhase), 0);
}

// SRS-CONS-004
#[test]
fn reads_that_place_a_run_boundary_differently_are_decided_by_position() {
    // One read deletes a T, the other substitutes it: their run lengths
    // never mix, and the read without an indel decides each position.
    let reference = "GTCAGCCCTTGACT";
    let deleting = rows("d", 0, reference, "GTCAGCCC-TGACT", "cccccccccccccc");
    let substituting = rows("s", 0, reference, "GTCAGCCCCTGACT", "cccccccccccccc");
    let decision = decide_stretch(&[deleting, substituting], 5, 9, reference.as_bytes());
    assert_eq!(decision.value.as_deref(), Some("CCCCT"));
    assert!(decision.runs.is_empty());
}
