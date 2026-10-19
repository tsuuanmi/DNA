---
id: PROP-0010
type: proposal
status: draft
owners: []
created: 2026-10-19
related-requirements: [SRS-CALL-004, SRS-CALL-005, SRS-CONS-008, SRS-API-014]
related-decisions: [ADR-0069, ADR-0071, ADR-0074, ADR-0076, ADR-0080]
implementation: []
---

# Proposal: Decompose a length mixture into its sequences

## Problem

A homopolymer that carries length heteroplasmy is a population of molecules with different
run lengths. After such a run, a Sanger read shows the sum of those molecules. That sum is the
same downstream sequence several times, each copy offset by the difference in run length. The
Sanger plugin marks that stretch `dephased` after a repeat (SRS-CALL-004), and a masked call
yields no base (ADR-0076, `dna.review/v1` `reads[].masked`). The consensus then decides the
run's length from weaker evidence (SRS-CONS-008: `confirmed_end`, `reference_frame`,
`phase_loss`), or leaves it undecided.

At the HVS-II C-tract (303–315) this is the largest single block of reviewer work:

- `estimated_length` checks are in 925 of 10,697 current samples (main `b49ff87`);
- solving them alone would clear 4.6 % more samples;
- the length DNA estimates matches the reviewers' in 47 % of those items, and the main
  error is the reference `C7` against a recorded `C9` (`309.1C 309.2C`).

The laboratory's requirement is to separate the length-heteroplasmy region into its
sequences, then compare each sequence independently against the rCRS to characterise how
they differ.

## Goals

- From the read's own signal, decompose a dephased stretch after a repeat run into the
  sequences it holds. Each component is an offset of the read's dominant sequence, with its
  fraction.
- Publish each component as its own sequence, compared independently with the reference, so a
  reviewer sees both sequences and their notation (e.g. `309.1C 315.1C` and
  `309.1C 309.2C 315.1C`).
- Let the consensus decide a run's dominant length from the component evidence where both
  strands agree, and keep the check where they do not.
- Data, not defaults:
  - the stretch comes from callability;
  - the offsets are bounded by the read's own repeat run;
  - the number of components is chosen against the read's own noise level;
  - no new fixed threshold.

## Non-goals

- A heteroplasmy verdict, a quantified minor allele as a call, or IUPAC output (SRS-VAR-009).
  The minor sequence is shown, not called.
- Report conventions. For example, the HVS-I C-run length that 16189C creates stays
  unrecorded by the profile convention.
- Reading single-base heteroplasmy inside the dephased stretch. A call the components do not
  explain is flagged, not decided.
- Porting Tracy (see Alternatives).

## Evidence and research

These are research spikes, out of repo (ADR-0065/0066), on current samples (batch date
≥ 2026-01), with reviewer truth only for scoring. Aggregates only.

**The tail is the read's own sequence, shifted.** At each mixed call after the HVS-II
C-tract, ask whether the second base equals the primary of the call d positions away. This
test is reference-free.

| truth C-tract length (303–315) | reads | winning offset | agreement ≥ 0.7 |
|---|---:|---|---:|
| 13 (as rCRS) | 18 | none consistent | 11 % (no mixture) |
| 15 (e.g. `309.1C 315.1C`) | 324 | d = −1 in 302 | 92 % |
| 16 | 243 | d = −1 in 143, +1 in 69 | 73 % |

**A component fit recovers the lengths.**

- Fit y_i = Σ_d f_d · onehot(S[i − d]), f_d ≥ 0, by non-negative least squares on the
  normalised four-channel heights at each call. S is the read's dominant sequence.
- On 800 HV2F reads, the dominant component's lengths match the truth for both runs (303–309
  and 311–315) in 98.8 % of reads.
- The dominant fraction has a median of 0.84. The minor sits at d = ±1 with fraction 0.1–0.3:
  a truth `C8` carries a minor `C9`.
- At HVS-I 16184–16193, with 16189 T, it is 99.9 %, essentially without a mixture.

**Each run from the strand that enters it in phase.** On 1,500 samples:

- the run 303–309 from HV2F and the run 311–315 from HV3R give both lengths right in 99.7 %;
- the strands agree on 303–309 in 1,494 samples, 4 of them wrong;
- the strands agree on 311–315 in 1,498 samples, all right.

A reverse read that reads 303–309 after the mixed exit alone is less reliable (92.6 % in one
draw), so the in-phase strand matters.

**Tracy** (`tracy decompose`) separates two alleles offset by an indel. It does so by:

- locating one breakpoint where the best-minus-second signal jumps;
- phasing the primary against the reference;
- testing every shift, counting calls the reference does not explain, and taking a
  MAD-bounded minimum;
- fitting allele fractions on a grid;
- aligning each allele to the reference separately.

The model and the per-allele comparison are what this proposal keeps.

## Proposed design

Increments are merged one at a time, each measured on its own.

### Increment 1: components as evidence, no change to any call

- **dna-sanger**, new `decompose` module, reference-free.
  - It works on every `dephased` segment that starts after a repeat run (`after_repeat`).
  - It takes the read's dominant sequence S over the segment, as its primary calls, plus the
    run before it.
  - It fits non-negative fractions over offsets bounded by the length of that repeat run.
  - It adds an offset only while that lowers the residual by more than the read's own noise
    level, as measured on its in-phase calls by signal processing.
  - Output per segment:
    - its components (offset, fraction);
    - the residual;
    - the calls the components leave unexplained.
- **dna-kernel**, the `ReadEvidence` contract.
  - A segment-level assertion `PhaseMixture { call_start, call_end, components: [{offset,
    fraction}], unexplained: [call] }`.
  - It is generic: any modality that sees a superposition of offsets can assert it. The core
    reads offsets and fractions, never Sanger labels.
- **dna-align.** Each component's sequence over the window is the dominant sequence with the
  run length shifted by its offset. Each is compared with the reference by the same aligner
  and placement rules, so each has its own edits.
- **dna**, review (`dna.review/v1`, new rules version). For each admitted read, publish its
  components:
  - span;
  - fraction;
  - sequence;
  - the component's calls in the profile notation.

  This is display only. No item, label or `clear` changes.

### Increment 2: the consensus decides the dominant length

- **dna-consensus.**
  - A new run-length evidence, `dominant_component`, slots into SRS-CONS-008 after
    `confirmed_end` and before `reference_frame`.
  - A run's length is the one the dominant component shows in the reads that enter that run
    in phase.
  - It is decided only when the two strands agree. Otherwise the run stays as today
    (`estimated_length` remains a check).
- **dna-call / decided.rs.** It reads the evidence label only through the published run
  evidence, as for `confirmed_end` today.
- **dna**, review. `estimated_length` no longer lists a run decided by `dominant_component`
  on both strands. The components stay visible.
- The guarantee that a masked call yields no base or indel (#147) is amended once, in the
  ADR: inside an asserted mixture, the dominant component may show a run's length. Bases
  downstream are not read in this increment.

### Increment 3 (later, separate decision): read the tail

The dominant component also reads the bases after the run (316–340, 16194–16200). Deciding
them is a separate change with its own gate.

## Alternatives considered

- **Port Tracy.** Rejected:
  - fixed parameters (peak ratio 0.33, window 25, jump 0.25, floor 10, 5·MAD, grid 0.01);
  - a single breakpoint per trace, where callability already locates the stretch;
  - phasing against the reference, where at mtDNA C-tracts neither length need be rCRS;
  - two alleles only.
- **Run-base exit survival curve.** Read the length distribution from the fraction of the
  run's base past its end (2026-10-17 research). It is kept as a cross-check: it gives a
  distribution but no sequences to compare with the reference.
- **Tune the existing rules** (`reference_frame`, `phase_loss`). These are the 47 % today.
  They decide without the signal that carries the length.
- **Leave it to reviewers.** This keeps the largest workload block.

## Validation plan

- **Unit tests** on synthetic mixtures:
  - a single component;
  - a 0.8/0.2 offset ±1;
  - three components;
  - a substitution inside the tail flagged as unexplained;
  - no mixture after a run.
- **Increment 1.** dna-eval (panel v5) outputs are byte-identical to main apart from the new
  review fields. Component lengths are reported against truth per region and strand.
- **Increment 2.** The gate on panel v5:
  - unreviewed and non-het do not rise on any slice;
  - strict and lod are listed block by block;
  - HVS-II length accuracy where both strands agree is ≥ the spike's 99.7 %.
- **The minor component.** Fraction and offset agree between the two strands of a sample.
  The truth records only the dominant length, so the minor is validated by strand agreement,
  not by truth.

## Rollout and rollback

Each increment is a PR with its plugin versions raised (sanger, consensus, review). Rollback
is a revert. Increment 2's rule is a new evidence step and can be withdrawn without touching
increment 1.

## Observability and operations

Logs record aggregate counts per read: mixture segments, components, unexplained calls. Per-call
arrays are not logged (SRS-CALL-006). Review and consensus outputs carry the evidence.

## Security and privacy

No new inputs or dependencies. Research outputs stay in `/tmp` and `results/`, with
aggregates only in docs.

## Compatibility and migration

`ReadEvidence`, `dna.review/v1` and `dna.consensus/v1` are revised in place before the
first release, as for ADR-0080.

## Decision

Pending.

## Implementation

Increment 1 is implemented ([ADR-0082](../decisions/adr/0082-phase-mixtures.md),
SRS-CALL-023, SRS-VAR-022), with these refinements of the design above:

- The fit covers every `dephased` segment that a repeat run opens, not only
  `after_repeat` ones, whose longer repeat missed the HVS-II C-tract. The
  opening run is the longest run near the segment's start, and the stretch
  ends where a run as long starts.
- A mixture is `explained` only when its residual is at most the read's own
  in-phase residual at `decompose_unexplained_quantile`. Only an explained
  mixture's components are compared with the reference; an unexplained one
  keeps its offsets and fractions.
- A component is the whole read with its run changed, so it places like the
  read. Its span covers the tract (the runs of the opening run's base one
  call apart before it in read order) and the stretch, may cross the origin,
  and its edits include the bases either side of it.
- The opening run starts at or after the trim interval's start.
- The review fields are a separate change on top of this one.

On panel v5, current slice, an explained component's entered run at 303–315
is the truth's in 95.6 % of HV2F reads (1,857 of 1,943), against the spike's
98.8 %, and in 99.3 % of HV3R reads (1,345 of 1,355). The gap is the dominant sequence being the read's primary calls,
which a heavily lagging copy can mislead; decoding it from the mixture is a
measured follow-up. HVS-I has no truth for its 16189C run length
(`hvs1_c_length`).
