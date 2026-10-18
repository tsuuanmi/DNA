---
id: PROP-0008
type: proposal
status: accepted
owners: []
created: 2026-10-17
related-requirements: [SRS-VAR-011, SRS-VAR-013, SRS-VAR-016, SRS-CONS-003, SRS-CONS-004, SRS-CONS-008, SRS-ALN-015]
related-decisions: [ADR-0047, ADR-0069, ADR-0070, ADR-0071, ADR-0074, ADR-0076, ADR-0077]
implementation: []
---

# Proposal: One responsibility per pipeline step

## Problem

The core caller should answer one question: given a read placed on the
reference, what does that read assert? A correct read, correctly placed, must
then call correctly, whatever its sequencing technology. Today `dna-core` also
interprets how a Sanger trace degrades, which makes its rules hard to reason
about and lets signal behaviour decide sequence content.

A read-only inventory of `dna-core` and `dna-kernel` at `f6b979e` found:

- **Four definitions of a trusted call, computed three times:**
  - the callable segments (unmasked and canonical);
  - the read-end margin with the ADR-0076 exemption and the support-veto
    scope;
  - the consensus's "clean", which rejects every veto whatever its scope and
    so disagrees with variant calling on `mixed_supporting_dna`;
  - overlap agreement, which ignores masks.
- **Sanger meaning carried by kernel types:**
  - `VetoScope`, `SupportVeto`, and `VetoSet`;
  - `MaskOnset::AfterRun`, which exists only for the ADR-0076 exemption;
  - the stage errors of ABIF decoding, signal processing, and callability.
- **A model of how a Sanger read loses phase, inside the consensus:**
  - run lengths no read bounds, estimated as `phase_loss` or
    `reference_frame`;
  - the handling of an `N` at a run's exit.
- **Run lengths asserted without bounds.** Outside a consensus stretch, a run
  that no trusted read bounds is decided position by position, which asserts
  its length. This is how DNA writes `303-309` at the reference length in about
  600 development-panel samples where the truth has `309.1C` or `309.2C`.

## Goals

1. Each step owns one question, and no step re-derives another's answer.
2. The modality alone decides how far each call can be trusted, and says so
   in modality-neutral verdicts. The core never interprets a modality label.
3. The core keeps one rule for content: **a read asserts only what it shows
   whole**, and a run length needs both of its ends shown.
4. An assumption or an estimate is never reported as a call.
5. The crate boundaries enforce the steps.

## Non-goals

- Changing the signal processing itself. Phase recovery, denoising, and
  length estimation remain later work in `dna-sanger`.
- Heteroplasmy (IUPAC) output.
- The laboratory's notation conventions, which stay in `dna-normalize` after
  the consensus (ADR-0077 §5).

## Evidence and research

Development panel v5, scratch builds never merged. The outputs are under
`results/architecture/`.

**R1: how today's consensus decides, and how often it is wrong.** By the
evidence for a run's length, the share of runs whose call differs from the
truth. These are upper bounds, since a notation-only difference counts:

| Run evidence | Runs | Differ |
|---|---:|---:|
| both ends shown in phase (`in_phase`) | 230,811 | 0.34 % |
| reference frame assumed (`reference_frame`) | 2,413 | 1.1 % |
| one read's clean start joined to another's clean end, run by run | 3,330 | 2.5 % |
| the same join, as a whole sequence | 204 | 30 % |
| an end shown by a masked anchoring call (`anchored_end`) | 3,288 | 3.8 % |
| phase-loss estimate (`phase_loss`) | 779 | 7.2 % |
| outside any stretch, no trusted read bounds the run (`303-309`) | 614 | 97 % |

About 99 % of the positions are decided by whole-sequence votes of reads that
show every run's ends.

**R2: a consensus that votes per position over normalized per-read variants.**
This was the simpler design: each read asserts its normalized variants, and
the sample takes each position's majority. It was measured in four forms. The
best (an anchoring flank counts as a run end, and each read is normalized
before the vote) raised the `haplotype_lod` not-concordant rate from 0.111 to
0.152: 1,302 samples got worse and 606 better. The failure is structural:

- reads describe one length-variant sequence through different alignments and
  assert different subsets of positions;
- a per-position vote then assembles a sequence that no read shows (for
  example `310C 311.1T 311.2C 311.3C`);
- per-read normalization cannot fully align the descriptions, since a read
  that does not read a whole window keeps its normalized form (SRS-NOM-010).

The vote unit must therefore be a stretch's whole sequence, as today. R2 also
showed that an end shown by a masked anchoring call is load-bearing. Without
it, `315.1C` is lost in about 5,100 samples, because no read shows both ends
of `311-315` in phase.

**R3: inventory.** `results/architecture/r3-inventory.md` lists each item with
its file and line, its class (Sanger meaning, neutral read logic, consensus
logic, vocabulary), and a proposed home.

## Proposed design

### Steps and owners

| Step | Crate | Question it owns |
|---|---|---|
| modality | `dna-sanger` | What did the instrument read, and how far can each call be trusted? |
| placement | `dna-align` | Where does the read lie on the reference? One 3'-most placement for equivalent indels (ADR-0047). |
| read assertion | `dna-call` | What does this placed read assert: its variants, and the positions it shows? |
| sample evidence | `dna-sample` | What do the reads show together, for review? |
| consensus | `dna-consensus` | Which sequence do the reads agree on, stretch by stretch? |
| notation | `dna-normalize` | How does the laboratory write it? |
| export | `dna-empop` | In which format does it leave? |
| pipeline | `dna` (facade) | Which steps run, in which order, for each workflow. |

`dna-sanger` depends only on `dna-kernel`. The core crates form a chain
(`dna-kernel` ← `dna-align` ← `dna-call` ← `dna-sample`, `dna-consensus`),
and only the facade composes the steps (ADR-0070).

### The modality boundary: per-call verdicts

`ReadEvidence` becomes version 2. Each call carries:

- its base and profile, as now;
- an **alignment instruction**: scored, anchoring (it places the read but
  shows no base), or unresolved;
- three **verdicts**, each trusted or withheld with a reason label reported
  verbatim:
  - **identity**: the call shows its base;
  - **insertion**: an inserted call shows that base was inserted;
  - **bound**: the call shows where a run ends.

The modality computes them once. It phrases the
ADR-0076 exemption on the read's own base runs, which it has before
alignment: a run's resolving calls on each side are trusted as bounds, and the
run's own calls stay trusted for a length even inside the signal-edge margin
when both bounds are trusted. Deletions then rest on the flanking calls' bound
verdicts, and substitutions on identity. This matches today's exemption
wherever the read's own runs are the aligned runs. It may differ in tandem
repeats, so step 4 measures exactly this and expects no difference; if one
remains, the core keeps a rule that maps the read's runs to the aligned ones.

Its inputs are everything it knows: masks, peak and quality floors, mixed
channels, a phase loss after a repeat run, and its own trace edges. The kernel loses `VetoScope`, `SupportVeto`, `VetoSet`,
`MaskOnset`, and the Sanger stage errors. The core never branches on a label.
There is one definition of trust, and the consensus, variant calling, sample
evidence, and the evaluation harness all read it.

### The core rule

A read asserts:

- the base of every position it shows with an identity verdict;
- an insertion or deletion only when the calls it rests on are trusted for it;
- a run's length only when it shows both of the run's ends. An end is shown by
  a call trusted as a bound that reads another base; a masked anchoring call
  can be such a bound (ADR-0074 `anchored_end`, which R2 found load-bearing).

Positions a read does not assert are not callable in that read. The neutral
rules that stay in the core:

- `run_boundary` (both ends);
- `unresolved_flank`;
- `placement_conflict`;
- an excluded call is not callable (SRS-VAR-019);
- the repeat-equivalent judgement of an insertion's calls (SRS-VAR-020);
- a placement margin, counted as today, around a physical end the modality
  does not vouch for (SRS-VAR-015), because an indel there cannot be placed
  by any technology. An `N` is handled by `unresolved_flank` (SRS-VAR-018),
  not by a margin.

The signal-edge margin around the trim ends and masked calls is the
modality's, expressed as verdicts. For reviewed sequences, which have no
masked calls and vouched ends, the outcome does not change.

### Consensus

- The vote unit stays a stretch's whole sequence, decided by run structure
  (ADR-0074).
- A run's length is decided only by reads that show both its ends. A
  `phase_loss` or `reference_frame` length, which no read bounds, stays in the
  consensus document with its label for review, and `call` treats it as not
  callable (option (i), approved). Deleting the estimators, and with them the
  handling of an `N` at a run's exit, is a later decision.
- A run that no trusted read bounds, while some read disagrees on it, seeds a
  stretch, so the position-by-position path never asserts a run length (the
  narrow trigger measured for `303-309`).
- `call` reports the consensus's decisions without placing the reads again
  (PROP-0006).

### What does not change

- The output contracts, except where a step's output changes as described.
- The laboratory notation and the EMPOP export.

## Sequencing

Each step below is its own change, measured on development panel v5 against
the main revision before it, with its plugin versions raised:

0. **Callable yield in the harness**: the callable positions inside the
   truth's analysed intervals, per sample and region. Without it, the yield
   cost of steps 2 and 3 cannot be measured.
1. **Independent fixes**, already approved:
   - `CORE_REASONS` gains `run_boundary`;
   - the normalization trade gap;
   - neutral wording in the core's requirements.
2. **PROP-0006**: `call` from the consensus's decisions. It is implemented on
   a paused branch; its last measurement was `haplotype_lod` 14 better, 0
   worse. It includes two rules in this proposal's spirit:
   - a length change of the run at a segment's edge is excluded, since no
     read shows its outer end;
   - record rules are kept, with no `run_boundary` waiver, which measured
     worse.
3. **Consensus lengths**:
   - `phase_loss` and `reference_frame` lengths are not callable in `call`
     (option (i), approved);
   - add the stretch seed for a disputed unbounded run.
4. **Per-call verdicts**: `ReadEvidence` v2, modality trust moved to
   `dna-sanger`, one trust definition in the core. The measurement must show
   identical outputs wherever the semantics do not change. The intended
   changes:
   - the consensus reads vetoes by their scope, as variant calling does
     (`mixed_supporting_dna`);
   - overlap agreement reads the unified trust instead of ignoring masks, so
     admitted read pairs and coverage may change;
   - the evaluation harness maps the verdict reasons verbatim, so its cause
     evidence labels may change while statuses do not;
   - `dna::error` wraps a `dna-sanger` error type for decoding, signal, and
     callability failures: a Breaking change for Rust callers, with CLI
     messages unchanged;
   - no difference in the ADR-0076 exemption (see the modality boundary).
5. **Crate split**: `dna-align`, `dna-call`, `dna-sample`, `dna-consensus`.
   Code moves plus the minimal public types; outputs stay identical:
   - a read-only public view of the alignment columns in `dna-align`;
   - trust as a `dna-call` output type, so the consensus no longer imports
     the private `ReadEligibility` and `CoreConfig`;
   - the consensus plugin declaring the contract it reads (`CalledRead`).

## Alternatives considered

- **A per-position variant consensus after per-read normalization.** Measured
  (R2) and rejected: it assembles sequences no read shows.
- **Strict bounds without anchoring ends.** Measured (R2): it loses `315.1C`
  across the panel.
- **A Sanger crate that wires the whole flow.** Rejected: the modality crate
  would depend on the core (ADR-0070), and each future modality would repeat
  the flow. The facade wires the steps.
- **Splitting the crates first.** Rejected: it would fix today's unclear
  boundaries in place. The semantics change first, then the code moves.

## Validation plan

- Each step's measurement on development panel v5:
  - every measure;
  - class by cause;
  - callable positions lost per sample, by region;
  - for steps 4 and 5, identical outputs apart from the stated changes.
- Unit tests:
  - the verdict mapping in `dna-sanger`;
  - the core rule on synthetic reads (bounded, anchored, and unbounded runs);
  - the stretch seed.
- The module-layer validator encodes the new crate map.

## Rollout and rollback

Steps land one at a time and each can be reverted alone. Steps 1 to 3 change
calls, and step 4 changes the evidence contract. Step 5 moves code only.

## Observability and operations

The verdict reasons are reported verbatim in the sample evidence, so a
reviewer sees why a call was withheld. The harness reads the same verdicts
instead of recomputing them.

## Security and privacy

No change: all processing stays local.

## Compatibility and migration

`ReadEvidence` v2 and any changed result contracts follow the versioning
policy (unreleased: revised in place, recorded as Breaking). No migration of
stored results.

## Decision

2026-10-17, user: accepted as written. In particular:

- the read-end margin is split: the core keeps it around `N` and around
  unvouched ends, and the modality's verdicts cover its signal edges;
- `phase_loss` and `reference_frame` follow option (i); deleting the
  estimators is a later decision;
- the steps run in the order 0 to 5.

2026-10-17, user: step 3 cuts only the length of an estimated run, its 3'
end, and keeps its composition callable; this narrows option (i) after the
measured variants below.

An ADR records the architecture when step 4 or 5 lands.

## Implementation

- Step 0: callable yield in the harness (#105).
- Step 1: the independent fixes (#103).
- Step 2: PROP-0006, `call` from the consensus's decisions (#106).
- Step 3: consensus lengths (`consensus` 5, `decided` 2):
  - a reference homopolymer of at least four bases that no read shows
    cleanly with both neighbours, while some read differs on it, seeds a
    stretch (SRS-CONS-004). Runs of two or three were measured and left out:
    they caught mixed signal at a heteroplasmic base (`195Y`, `16093Y`);
  - a run of length evidence `phase_loss` or `reference_frame` keeps its site
    even where its length equals the reference (SRS-CONS-007), so the label
    reaches `call`;
  - `call` of the consensus document leaves such a run's length uncallable:
    its 3'-most position, the deleted positions beside it included, where a
    change of its length is written, leaves the callable segments, and a
    change there is excluded as `run_boundary` (SRS-VAR-021). The run's other
    positions stay callable, since the reads show its composition: this
    narrows option (i), which cut the whole run;
  - a site's own alignment costs a gap extension 2, so a length change in
    each of two runs stays two indels rather than one gap across their
    boundary and a substitution.

  Measured on development panel v5 against main (discordant share of all
  samples, then samples better and worse):

  | Measure | Main | Step 3 | Better | Worse |
  |---|---:|---:|---:|---:|
  | primary | 0.0910 | 0.0863 | 81 | 5 |
  | strict | 0.2146 | 0.2476 | 8 | 549 |
  | haplotype | 0.1251 | 0.0847 | 666 | 4 |
  | `haplotype_lod` | 0.1098 | 0.0688 | 677 | 5 |

  - Most of the `haplotype_lod` gain is the HVS-II `303-309` length that the
    consensus asserted unlabelled at the reference length.
  - The strict cost is coverage: the truth's `309` insertions after a
    `303-309` run whose length is an estimate are no longer callable.
  - Callable yield: HV1 0.9895 to 0.9892, HV2 0.9898 to 0.9894, HV3 0.9899
    to 0.9901.
  - Variants measured on the way (`haplotype_lod` better and worse, strict
    better and worse):

    | Variant | `haplotype_lod` | strict |
    |---|---:|---:|
    | seed and label, no cut (option (ii)) | 82 / 24 | 64 / 4 |
    | cut the whole run (option (i)) | 687 / 4 | 32 / 2,073 |
    | cut the whole run, `phase_loss` only | 129 / 3 | 24 / 548 |
    | cut the 3' end (chosen) | 677 / 5 | 8 / 549 |
    | cut from the shortest in-phase length to the 3' end | 677 / 5 | 8 / 549 |

    The whole-run cut also made about 1,840 samples' HVS-I `16183C` and
    `16189C` uncallable, although the reads show that composition. The
    shortest in-phase length almost always equals the decided one, so the
    last variant changed nothing and would have needed a new consensus
    field.
  - Left open: 41 samples whose primary status the whole-run cut fixed stay
    as before, with `310C`, `310DEL 311DEL`, or `310DEL`–`312DEL` just past
    the estimated `303-309` run. They are not new; they are the next item in
    that window.
- Step 4, part 1: signal-edge verdicts and the read start (`core` 9,
  `consensus` 6, `sanger` 4). The position-73 coverage loss moved here, since
  the core's margin around modality-masked calls was what kept strong
  leading calls untrusted:
  - the margin split. `CallEvidence` carries the modality's signal-edge
    verdict (`edge`), and the core keeps only the margin at a read end the
    modality does not vouch for (SRS-VAR-013, SRS-CALL-013). Before any rule
    change, development panel v5 was byte-identical to main in every
    harness output;
  - the leading-edge verdict. Before the callable span, a masked call with
    one qualifying channel and a primary peak at least the read's median is
    vouched for by its base only (SRS-CALL-014). Such a call supports a
    substitution but no insertion or deletion and no run's end, and where a
    neighbour is not an unmasked call trusted for more than its base and its
    base is the reference base next to it, it may sit one position off and
    places nothing (`unresolved_flank`, SRS-VAR-018, SRS-CONS-003); nor
    does it show a position of a reference homopolymer of at least four
    bases the read does not bound with calls trusted for more than their
    base, since that position would carry a run length no read shows;
  - each call's second base and its share (SRS-CALL-015), evidence only.

  Measured on development panel v5 against main (samples better and worse;
  the shift rule's effect in brackets, without it):

  | Measure | Better | Worse |
  |---|---:|---:|
  | primary | 3 | 103 (199) |
  | strict | 329 | 78 (142) |
  | `haplotype_lod` | 2 | 3 (39) |

  - 553 of the 745 samples whose truth `73G` DNA left uncallable now call
    it; 487 of them are `haplotype_lod`-concordant.
  - Of the primary and strict regressions, 100 and 74 are only new DNA
    alleles outside the truth's analysed interval (mostly `73G` where the
    reviewers start at `79`).
  - Leading calls vouched for by their base, against the truth inside its
    intervals: at truth variants 1,646 used and right, 8 set aside by the
    shift rule (`93G`, `16362C`, `16357C`, which become uncallable); at
    reference positions 12,778 used and right, 6 used and wrong, 193 set
    aside. HVS-I, HVS-II and HVS-III read starts all contribute.
  - Callable yield: HV2 0.9894 to 0.9910, HV3 0.9901 to 0.9905, HV1 as
    before.
  - Placement: no read changes orientation; 38,813 reads extend at their
    start, 472 move an end by a few positions, and 34 reads become ambiguous
    while 17 others resolve (51 samples).
  - The run rule turned one sample (`573`-`576`) from discordant to not
    callable and changed nothing else.
  - Left open (`haplotype_lod`, 3 samples): strong single-channel G calls at
    the HVS-I read start (`16025G`, `16026G`), probably dye blobs, a false
    variant class the rescue introduces; and a deletion placed differently
    at `16038`, a side effect of the earlier informative interval.
- Step 4, part 2a: `ReadEvidence` v2, with outputs identical. Each call
  carries its placement (scored, or masked as anchoring or unresolved with
  its reason), its verdicts on its base, on an indel, and on a run's end,
  each with the modality's reasons in the read's reason order and its
  signal edge, its mark in a bounded run before a post-repeat phase loss,
  and its second base. The kernel loses `CallMask`, `MaskOnset`, `VetoSet`,
  `SupportVeto`, and `VetoScope`. The Sanger plugin now finds the bounded
  runs of ADR-0076 on the read's own calls (SRS-CALL-016); the core keeps
  the alignment-dependent half. Development panel v5: every harness output
  byte-identical to main, cause evidence included; real traces identical.
- Step 4, part 2b: the intended changes (`core` 11, `consensus` 7), each
  measured on development panel v5 against main:
  - verdict reasons fail closed (a verdict refuses a reason it cannot hold),
    and a masked placement carries its reason itself: no output change;
  - the consensus reads reasons by scope, as variant calling does: a base
    its verdict on the base, an inserted base its verdict on an indel, and a
    call beside an insertion or deletion none (SRS-CONS-003);
  - an insertion or deletion on a call vouched for by its base only is
    excluded as `base_only`, a new core reason, rather than `read_end`: no
    harness change; on the 52 real samples, 44 such exclusions;
  - overlap agreement compares only bases the core trusts (SRS-SAMPLE-008):
    no harness change; on the 52 real samples, conflicts fall from 1,156 to
    11 and eligible pairs from 137 to 85, the lost pairs being short
    overlaps at read ends that the signal edges now leave below the
    comparable minimum;
  - review found that reading flanks by scope let a deletion be decided
    between undecided positions and then read as callable reference. Two
    fixes: every candidate excluded before it is observed (an indel no flank
    reads, an unresolved allele) leaves its positions undecided
    (SRS-VAR-019), and a site that changes the length stays undecided while
    a position beside it is undecided, unless reads show its runs' lengths
    (SRS-CONS-009). A scan of the panel for decided differences that are
    neither reported nor uncallable finds none.

  Net against main: `haplotype_lod` 3 better and none worse, primary 2 and
  1, strict 1 and 4 (three `523DEL 524DEL` that the two strands read
  differently, now contested; one `310DEL 311DEL` decided from the read
  that bounds the run); callable yield HV2 0.9910 to 0.9909, HV3 0.9905 to
  0.9904; causes `dna` 1,212 to 1,178, `technology` 519 to 556, `lab` 370
  to 377.
