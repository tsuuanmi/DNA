//! The consensus passes over more than one read per strand (ADR-0090): two
//! forward reads and one reverse read, and two of each. Each test pins what
//! the pass does today with N reads; where ADR-0090 leaves a step pending,
//! the test says so.

use super::*;

/// Three reads of `SHORT`-shaped rows, every position `edge` calls from the
/// read's ends: `forward` reads on the forward strand, the others reverse.
fn strands(reads: &[(&str, &str, &str, bool)], edge: usize) -> Vec<ReadObservations> {
    const SHORT: &str = "ACGTACGTAC";
    reads
        .iter()
        .map(|(name, read, kinds, forward)| {
            let mut observed = rows(name, 0, SHORT, read, kinds);
            observed.forward = *forward;
            observed.edges = (0..SHORT.len()).map(|position| (position, edge)).collect();
            observed
        })
        .collect()
}

const CLEAN: &str = "cccccccccc";
const SHORT_BASES: &[u8] = b"ACGTACGTAC";

/// Sets aside the dissent of `observed` once, and the decision at 4 after.
fn set_aside_once(observed: &mut [ReadObservations]) -> (bool, Decision) {
    let decided = decide_positions(observed, SHORT_BASES, 4);
    let changed = set_aside_dissent(observed, &decided, SHORT_BASES, 12, 40);
    (changed, decided_at(observed, 4))
}

// SRS-CONS-011, ADR-0090
#[test]
fn an_explained_second_is_set_aside_with_two_reads_on_the_supporting_strand() {
    // 2F + 1R: two forward reads of the reference, one reverse read with
    // an explained second at 4.
    let mut observed = strands(
        &[
            ("f1", "ACGTACGTAC", CLEAN, true),
            ("f2", "ACGTACGTAC", CLEAN, true),
            ("r", "ACGTTCGTAC", CLEAN, false),
        ],
        50,
    );
    observed[2].seconds.insert(4, 'A');
    assert_eq!(decided_at(&observed, 4).opposing, ["r"]);
    let (changed, after) = set_aside_once(&mut observed);
    assert!(changed);
    assert_eq!(after.value.as_deref(), Some("A"));
    assert!(after.opposing.is_empty());
    assert_eq!(after.supporting, ["f1", "f2"]);
}

// SRS-CONS-011, ADR-0090
#[test]
fn every_opposing_read_must_be_explained_for_any_to_be_set_aside() {
    // 2F + 2R: both reverse reads differ at 4.
    let reads = || {
        strands(
            &[
                ("f1", "ACGTACGTAC", CLEAN, true),
                ("f2", "ACGTACGTAC", CLEAN, true),
                ("r1", "ACGTTCGTAC", CLEAN, false),
                ("r2", "ACGTTCGTAC", CLEAN, false),
            ],
            50,
        )
    };
    // Two against two is a tie: contested, and nothing to set aside.
    let mut tied = reads();
    tied[2].seconds.insert(4, 'A');
    tied[3].seconds.insert(4, 'A');
    assert_eq!(decided_at(&tied, 4).state, SiteState::Contested);
    assert!(!set_aside_once(&mut tied).0);
    // With a third forward read the reference wins; both explained reverse
    // reads are set aside together.
    let third = strands(&[("f3", "ACGTACGTAC", CLEAN, true)], 50);
    let mut both = reads();
    both.extend(third.clone());
    both[2].seconds.insert(4, 'A');
    both[3].seconds.insert(4, 'A');
    let (changed, after) = set_aside_once(&mut both);
    assert!(changed);
    assert!(after.opposing.is_empty());
    // One opposing read with a peak its modality leaves unexplained keeps
    // both standing.
    let mut one = reads();
    one.extend(third);
    one[2].seconds.insert(4, 'A');
    one[3].seconds.insert(4, 'A');
    one[3].unexplained.insert(4);
    assert!(!set_aside_once(&mut one).0);
    assert_eq!(decided_at(&one, 4).opposing, ["r1", "r2"]);
}

// SRS-CONS-011, ADR-0090
#[test]
fn opposing_reads_on_both_strands_are_never_set_aside() {
    // 2F + 1R with one forward read and the reverse read opposing: no
    // supporting read lies on a strand other than every opposing read's.
    let mut observed = strands(
        &[
            ("f1", "ACGTACGTAC", CLEAN, true),
            ("f2", "ACGTTCGTAC", CLEAN, true),
            ("r1", "ACGTACGTAC", CLEAN, false),
            ("r2", "ACGTACGTAC", CLEAN, false),
            ("r3", "ACGTTCGTAC", CLEAN, false),
        ],
        50,
    );
    for read in [1, 4] {
        observed[read].seconds.insert(4, 'A');
    }
    assert_eq!(decided_at(&observed, 4).opposing, ["f2", "r3"]);
    assert!(!set_aside_once(&mut observed).0);
}

/// Reads of the insertion site `ACGATTCAGC`: forward reads that show
/// `TTTC` at 4-6 and reverse reads that read the site through.
fn insertion_reads(forward: usize, reverse: usize) -> Vec<ReadObservations> {
    const SITE: &str = "ACGATTCAGC";
    let mut observed = (1..=forward)
        .map(|index| {
            rows(
                &format!("f{index}"),
                0,
                "ACGATT-CAGC",
                "ACGATTTCAGC",
                "ccccccccccc",
            )
        })
        .collect::<Vec<_>>();
    for index in 1..=reverse {
        let mut read = rows(&format!("r{index}"), 0, SITE, SITE, "cccaaaaccc");
        read.forward = false;
        observed.push(read);
    }
    observed
}

// SRS-CONS-015, ADR-0090
#[test]
fn a_one_strand_insertion_is_contested_whatever_the_read_counts() {
    const SITE: &[u8] = b"ACGATTCAGC";
    for (forward, reverse) in [(2, 1), (2, 2), (3, 1)] {
        let observed = insertion_reads(forward, reverse);
        let after = contested_at(&observed, SITE, 5);
        assert_eq!(
            (after.value.as_deref(), after.state),
            (Some("TTC"), SiteState::Contested),
            "{forward}F + {reverse}R"
        );
        let names = |prefix: char, count: usize| {
            (1..=count)
                .map(|index| format!("{prefix}{index}"))
                .collect::<Vec<_>>()
        };
        assert_eq!(after.supporting, names('r', reverse));
        assert_eq!(after.opposing, names('f', forward));
    }
    // ADR-0090 pending: the pass reads strands, not counts, so one reverse
    // read reading through contests three forward reads as it does one.
}

// SRS-CONS-008, ADR-0090
#[test]
fn reads_entering_a_run_from_both_ends_confirm_it_in_any_number() {
    // As in `a_run_two_reads_enter_from_its_ends_takes_their_agreed_count`:
    // reads entering the merged tract from its left end, and from its right.
    let left = |name: &str| {
        rows(
            name,
            0,
            TRACT,
            "GTCAAAACCCCCCCCCCNTG",
            "cccccccccccccccccmmm",
        )
    };
    let right = |name: &str| {
        rows(
            name,
            0,
            TRACT,
            "GTCAAANCCCCCCCCCCATG",
            "mmmmmmmccccccccccccc",
        )
    };
    // 2 + 1 and 2 + 2 reads that agree on ten calls confirm the run.
    for observed in [
        vec![left("f1"), right("r1"), left("f2")],
        vec![left("f1"), right("r1"), left("f2"), right("r2")],
    ] {
        let decision = decide_stretch(&observed, 7, 16, TRACT.as_bytes());
        assert_eq!(decision.runs, [run('C', 10, RunEvidence::ConfirmedEnd)]);
        assert_eq!(decision.supporting.len(), observed.len());
    }
    // A read that loses phase inside the stretch gives no entry and does
    // not block the confirmation.
    let inside = rows(
        "r2",
        0,
        TRACT,
        "GTCAAANCCCCCCCCCCATG",
        "mmmmmmmmcccccccccccc",
    );
    let decision = decide_stretch(&[left("f1"), right("r1"), inside], 7, 16, TRACT.as_bytes());
    assert_eq!(decision.runs, [run('C', 10, RunEvidence::ConfirmedEnd)]);
    assert_eq!(decision.uninformative, ["r2"]);
}

// SRS-CONS-014, ADR-0090
#[test]
fn a_cross_strand_pair_needs_exactly_one_decode_per_strand() {
    use crate::length::{Run, Source, View, in_phase};
    use crate::noise::cross_strand;
    let view = |forward: bool, source: Source| View {
        source,
        forward,
        top: 8,
        next: Some(9),
        p: if forward { 0.8 } else { 0.7 },
        se: Some(0.03),
    };
    // A second forward read in phase is no decode: the pair stands.
    let views = [
        view(true, Source::Direct),
        in_phase(true, 8),
        view(false, Source::Derived),
    ];
    assert!(cross_strand(&views, Run::Before).is_some());
    // Two forward decodes and one reverse give none.
    let views = [
        view(true, Source::Direct),
        view(true, Source::Derived),
        view(false, Source::Direct),
    ];
    assert!(cross_strand(&views, Run::Before).is_none());
}

// SRS-CONS-013, ADR-0090
#[test]
fn the_two_strand_rule_takes_every_read_of_a_strand()
-> std::result::Result<(), Box<dyn std::error::Error>> {
    // As in `both_strands_decide_a_window_run_left_as_an_estimate`, with a
    // second forward read: C7 T C5 at 11-23, the anchor T at 18.
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
    let reverse = read(
        "r",
        &reverse_sequence,
        decoded('G', 'A', sequence.len() - 23, 5, 14),
    )?;
    let settle =
        |reads: &[&CalledRead]| -> std::result::Result<Decision, Box<dyn std::error::Error>> {
            let decision = Decision {
                value: Some("CCCCCCCTCCCCC".into()),
                state: SiteState::Called,
                supporting: Vec::new(),
                opposing: Vec::new(),
                uninformative: Vec::new(),
                runs: vec![
                    run('C', 7, RunEvidence::PhaseLoss),
                    run('T', 1, RunEvidence::InPhase),
                    run('C', 5, RunEvidence::InPhase),
                ],
            };
            let mut decided = BTreeMap::from([(10, (22, decision))]);
            let (call, consensus) = (config()?, config::tests::config());
            two_strand(
                &mut decided,
                reads,
                window_reference.sequence.as_bytes(),
                consensus.length_rule(),
                |first, last| {
                    read_runs_over(reads, &window_reference, &call, &consensus, first, last)
                },
            );
            Ok(decided.remove(&10).ok_or("missing stretch")?.1)
        };
    // 2F + 1R: both forward reads decode C8 T C5; the rule decides C8.
    let (f1, f2) = (
        read("f1", sequence, decoded('C', 'T', 10, 8, 14))?,
        read("f2", sequence, decoded('C', 'T', 10, 8, 14))?,
    );
    let decision = settle(&[&f1, &f2, &reverse])?;
    assert_eq!(decision.runs[0], run('C', 8, RunEvidence::TwoStrand));
    // A forward read that decodes C7 puts another length on top: every
    // view must agree, so the run keeps its estimate.
    let other = read("f2", sequence, decoded('C', 'T', 10, 7, 13))?;
    let decision = settle(&[&f1, &other, &reverse])?;
    assert_eq!(decision.runs[0], run('C', 7, RunEvidence::PhaseLoss));
    Ok(())
}

/// A small deterministic generator (xorshift64*), so the property test
/// needs no dependency and replays exactly.
struct Generator(u64);

impl Generator {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    fn below(&mut self, bound: usize) -> usize {
        usize::try_from(self.next() % bound as u64).unwrap_or(0)
    }
}

/// Two to five reads of `SHORT` on random strands, each with up to two
/// random substitutions, masked calls, explained or unexplained seconds,
/// and distances from its ends.
fn random_reads(generator: &mut Generator) -> Vec<ReadObservations> {
    const SHORT: &str = "ACGTACGTAC";
    let count = 2 + generator.below(4);
    (0..count)
        .map(|index| {
            let mut read = SHORT.chars().collect::<Vec<_>>();
            let mut kinds = vec!['c'; SHORT.len()];
            for _ in 0..generator.below(3) {
                let position = generator.below(SHORT.len());
                read[position] = ['A', 'C', 'G', 'T'][generator.below(4)];
                if generator.below(4) == 0 {
                    kinds[position] = 'a';
                }
            }
            let read = read.into_iter().collect::<String>();
            let kinds = kinds.into_iter().collect::<String>();
            let mut observed = rows(&format!("read{index}"), 0, SHORT, &read, &kinds);
            observed.forward = generator.below(2) == 0;
            for position in 0..SHORT.len() {
                observed.edges.insert(position, generator.below(60));
                if generator.below(5) == 0 {
                    observed
                        .seconds
                        .insert(position, ['A', 'C', 'G', 'T'][generator.below(4)]);
                }
                if generator.below(10) == 0 {
                    observed.unexplained.insert(position);
                }
            }
            observed
        })
        .collect()
}

/// The clean or trusted observations left: what a pass can withdraw.
fn standing(observed: &[ReadObservations]) -> usize {
    observed
        .iter()
        .flat_map(|read| read.positions.values().chain(read.junctions.values()))
        .filter(|observation| {
            observation.clean
                || observation
                    .chars
                    .iter()
                    .any(|char| char.trust > Trust::Nothing)
        })
        .count()
}

// SRS-CONS-011
#[test]
fn a_set_aside_keeps_each_decision_and_the_passes_end() {
    let mut generator = Generator(0x9E37_79B9_7F4A_7C15);
    let mut set_asides = 0;
    for case in 0..2_000 {
        let mut observed = random_reads(&mut generator);
        let mut decided = decide_positions(&observed, SHORT_BASES, 4);
        let mut rounds = 0;
        loop {
            let before = standing(&observed);
            let previous = decided.clone();
            let aside = set_aside_dissent(&mut observed, &decided, SHORT_BASES, 12, 40);
            if !aside {
                break;
            }
            // Each set-aside withdraws at least one observation, so the
            // passes end within as many rounds as there are observations.
            assert!(standing(&observed) < before, "case {case}");
            set_asides += 1;
            decided = decide_positions(&observed, SHORT_BASES, 4);
            // A read set aside no longer opposes, and a decided site keeps
            // its value.
            for (first, (last, decision)) in &previous {
                let Some((_, after)) = decided.get(first).filter(|(end, _)| end == last) else {
                    continue;
                };
                if decision.value.is_some() && !decision.opposing.is_empty() {
                    assert_eq!(after.value, decision.value, "case {case} at {first}");
                }
            }
            rounds += 1;
            assert!(rounds <= before, "case {case}: the passes do not end");
        }
    }
    // The cases exercise the pass (153 set-asides with this seed).
    assert!(set_asides >= 50, "only {set_asides} set-asides");
}

// SRS-CONS-005
#[test]
fn ties_of_three_or_more_sequences_decide_nothing() {
    let observations = |votes: &[(&'static str, &str)]| {
        votes
            .iter()
            .map(|(read, sequence)| (*read, Some((*sequence).to_owned())))
            .collect::<Vec<_>>()
    };
    // One read each for three sequences, the reference among them.
    let three = decide(&observations(&[("a", "T"), ("b", "G"), ("c", "C")]), "C");
    assert_eq!((three.value, three.state), (None, SiteState::Contested));
    assert_eq!(three.opposing.len(), 3);
    // Two against two against one: no strict majority.
    let two_two_one = decide(
        &observations(&[("a", "T"), ("b", "T"), ("c", "C"), ("d", "C"), ("e", "G")]),
        "C",
    );
    assert_eq!(
        (two_two_one.value, two_two_one.state),
        (None, SiteState::Contested)
    );
    // Three against three is no majority; two against one against one is.
    let three_three = decide(
        &observations(&[
            ("a", "T"),
            ("b", "T"),
            ("c", "T"),
            ("d", "C"),
            ("e", "C"),
            ("f", "C"),
        ]),
        "C",
    );
    assert_eq!(three_three.value, None);
    let two_one_one = decide(
        &observations(&[("a", "T"), ("b", "T"), ("c", "C"), ("d", "G")]),
        "C",
    );
    assert_eq!(
        (two_one_one.value.as_deref(), two_one_one.state),
        (Some("T"), SiteState::Called)
    );
    assert_eq!(two_one_one.opposing, ["c", "d"]);
}
