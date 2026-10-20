---
id: PROP-0012
status: draft
created: 2026-10-20
summary: "a library crate `dna-trace` re-processes the trace windows the Sanger plugin flags after a repeat run: it models the read's own base signatures and leakage, decodes the length mixture, and gives each run's components with fractions and a margin; the consensus decides a run's length only when the decode is clear and both strands agree, and otherwise flags it with the components, so no length is passed unflagged on a weak signal."
---

# Proposal: Trace signal engineering in its own crate

## Problem

The user's measure (2026-10-20) is non-heteroplasmy errors that pass with no
flag. Heteroplasmy stays optional: detecting it is a bonus, and an unflagged
one counts as a pass. A block where the truth reads `N` goes to the reviewer
re-check list. On development panel v5's current data with the overnight
stack (#175, #176), 475 blocks are unflagged errors: 366 heteroplasmy, 73
non-het, and 36 truth-`N`. The non-het blocks in this proposal's scope are:

- 30 in the HVS-II polyC window, at 309: a wrong length passed with no check;
- 16 at 573 (`CC/C`): the 2026-10-18 adjudication found most of them truth
  errors, and they are on the reviewer re-check list;
- possibly some of the 12 HVS-I `reads_show_reference` blocks around 16193,
  not yet adjudicated.

The largest block of review work left on current data is homopolymer length,
where the trace is a superposition of molecules whose runs differ in length.
On development panel v5 with the overnight stack (#175, #176), 2,871 of 10,695
current samples are not clear. 2,013 of them are held by a single zone. The
largest zones are:

- the HVS-II polyC window (303–315): 665 samples;
- the end of HVS-III (568–576, the same class at 568–573): 149;
- the HVS-I window (16180s): 195, mostly real heteroplasmy or real errors.

[[PROP-0010]] increment 1 ([[ADR-0082]]) decomposes such a stretch, as
evidence only, inside `dna-sanger`. Two later measurements set the limits of
that approach:

- **A one-strand rule is unsafe.** A run length decided from one read's
  dominant component (draft #174) changed 41 runs on the panel, all on one
  strand. It raised samples with an unreviewed difference (current 440 → 442,
  older 174 → 175).
- **Two strands decide only when the decode is clear.** The reference-free
  two-strand decoder spike covers 1,058 hard samples with usable truth (both
  strands in 979). "Hard" means DNA's 303–315 was wrong or undecided on
  `main` on 2026-10-19; the easy set is 500 samples DNA got right. The spike
  decoded only HVS-II, not 16184–16193 (HVS-I) or 568–573 (HVS-III).
  - **run 311–315:** the reverse read enters it in phase. Its decoded length
    is right in 90.8 % of hard reads, and in 99.7 % when the strands agree.
  - **run 303–309, per read:** primary calls are right in 9.8 % (HV2F) and
    20.9 % (HV3R). The decoder's dominant length is right in 55.8 % and
    72.4 %.
  - **run 303–309, both strands:** they agree in 78 %. With a decisive
    margin scaled to each read's own fit error (z), coverage and accuracy
    trade off:

    | z | Coverage | Accuracy |
    |---|---|---|
    | ≥ 2 | 57 % | 96.0 % |
    | ≥ 4 | 38 % | 98.3 % |
    | ≥ 6 | 11 % | 99.1 % |

  - **What remains wrong** is almost all near-ties of a real length
    heteroplasmy, where the recorded length is a reporting convention.
- **The signal is already processed.** DNA decodes the instrument's analysed
  channels (`DATA.9-12`, `PLOC.2`), which already carry baseline, crosstalk,
  and mobility correction. Re-cleaning them adds little. What pays is
  model-side: an empirical per-read signature for each base, and a null for
  the share a neighbour's peak tail leaks into the next call.

## Goals

The user's goal (2026-10-20): unflagged non-het errors in auto-passed zones
must never rise and should fall; the number of auto-passes comes second. The
order of preference:

1. auto-passed and right;
2. not auto-passed, but flagged with what the reviewer needs;
3. auto-passed, wrong, and unflagged: far worse than either.

So:

- Move trace signal engineering into its own crate, `dna-trace`, testable on
  synthetic traces and reusable by any dye-terminator modality.
- Decode a length mixture into per-run components: each run's length in each
  molecule, its signal fraction, the margin of the decision, and whether the
  components explain the stretch.
- Decide a run's length only when the decode is clear and both strands
  agree. Otherwise publish a flag with the components, so the reviewer
  decides in seconds: "dominant C8 ~70 %, minor C9 ~30 %, strands agree".
- Take every level from the read's or the run's own data, or from validated
  configuration with provenance. No fixed cut.

## Non-goals

- A reporting convention, such as the shortest length whose fraction is at
  least 0.6 of the dominant's (the "major" rule). It is the lab's decision,
  not the decoder's. The decoder reports components; the user decides
  whether any convention is applied, and where (`dna-normalize`).
- Heteroplasmy levels. Fractions are signal shares: stutter and
  heteroplasmy cannot be separated in Sanger signal (SRS-VAR-009).
- Raw-channel processing (`DATA.1-4`) or re-cleaning the analysed channels
  of the whole read. That is stage 2, only if research shows it pays.
- Reading the tail after a mixture base by base ([[PROP-0010]] increment 3).

## Evidence and research

- [[PROP-0010]] increment 1 and [[ADR-0082]]: the decomposition, its
  explained gate, and the per-component comparison.
- Draft #174: a one-strand dominant-component run length, rejected on the
  gate.
- The two-strand decoder spike (Researcher scratchpad, `decode/`, decoder v2:
  joint NNLS over molecule lengths with exact DP downstream), on 1,077 hard
  and 500 easy current samples; the numbers are under Problem.
- A literature and code survey (Researcher scratchpad, `signal-lit/`):
  - Every CRISPR mixed-trace tool fits shifted copies of an empirical control
    trace (TIDE, ICE, DECODR, PtWAVE).
  - PtWAVE chooses the number of components by BIC, without a threshold, and
    reports that Lasso loses real low components.
  - Forensic practice (SWGDAM 2024, Berger 2011, Sturk-Andreaggi 2020, EMPOP):
    - a length heteroplasmy must be seen on both strands;
    - the dominant length is read from a marker after the tract;
    - labs legitimately differ in the length they record;
    - shares are signal, not heteroplasmy levels.
- Planned research experiments, measured per strand and on the two-strand
  rule:
  - E1: empirical signature and leakage null;
  - E2: run-length proposals chosen by BIC;
  - E3: run-level context templates;
  - E4: height, area, or local normalisation;
  - E5: bootstrap intervals and a near-equal mark.

## Proposed design

### Crate boundary

`dna-trace` is a library crate with no plugin descriptor
([[ADR-0084]]). It depends on `dna-kernel` only, and `dna-sanger` depends on
it.

- **Input:** plain data, not ABIF:
  - the four analysed channels;
  - the call positions;
  - the read's primary calls;
  - the call ranges the modality flags.
- **Output:** typed results the modality turns into kernel assertions. The
  core never reads `dna-trace`.

ABIF decoding, basecalling, callability, and stage orchestration stay in
`dna-sanger`.

### Stage 1: only the flagged windows

`dna-trace` re-processes only the stretches the Sanger plugin flags:

- dephased or `post_homopolymer` segments after a repeat run;
- the run tract before them (the window [[ADR-0082]] already fits).

Outside them nothing changes, so each increment's risk stays local and
measurable per zone.

### The contract the core reads

The kernel `PhaseMixture` ([[ADR-0082]]) stays the one assertion, extended
rather than duplicated. Per run of its tract it gains:

- the run's length in each component, with that component's fraction;
- the margin of the dominant length against the next, scaled to the read's
  own fit error;
- whether the top two lengths are near-equal.

The core still sees no signal array.

### The decision and the flag

The consensus combines the two strands' assertions per run:

- **Decided** (a new run-length evidence) only when every read that decodes
  the run is explained, clear by its margin, and agrees, and reads of both
  strands decode it. The margin is configuration with provenance, or derived
  from the run's data.
- **Otherwise a check.** A review label carries, per read: the strand, the
  dominant length and fraction, the minor lengths and fractions, the margin,
  the near-equal mark, and explained or not. Then whether the strands agree,
  and whether the profile's convention records that run's length. A length
  heteroplasmy always stays flagged; it is never decided away.

### Increments

Each increment is one PR, measured on development panel v5 and on the
readiness panel `results/test-panel/readiness-v2.tsv` at a logged use. That
panel holds 20,140 samples: 12,731 from 148 batch-disjoint new batches
(2026-09/10, `fresh_batch`) and 7,409 from dev-batch remainders
(`dev_batch_remainder`), reported separately. It is used for aggregates
only, never inspected per sample while tuning, and each use is logged in
`results/test-panel/readiness-log.md`. A result that drives a rule change
turns it into development data.

1. **Move.** `dna-sanger/src/decompose` moves into `dna-trace`, unchanged.
   The panel must be byte-identical. An optional debug dump of the window's
   channels and fit serves research; it is not a contract.
2. **Model.** Empirical per-read base signatures and a leakage null replace
   the one-hot model in the fit (E1). The number of components is chosen by
   BIC (E2), so the explained gate's quantile is the only key left. Each
   read's per-run lengths, fractions, margin, and near-equal mark go on
   `PhaseMixture`. Evidence only.
3. **Flag.** The review shows the per-run components and strand agreement
   for every mixture run. It is also the first increment that reduces
   unflagged errors: a run that passes with no check today becomes a check
   whenever its decode is not "clear, and both strands agree with DNA's
   call", that is, when the strands disagree, the top two are near-equal,
   or only one strand decodes it. Expect more checks. The gate wants
   unflagged non-het blocks in the HVS-II polyC window to fall from 30, with
   precision unchanged. Reported: runs newly flagged, and
   how many of them hold a truth difference. No decision changes yet.
4. **Decide when clear.** The consensus decides a run's length on "clear and
   both strands agree", with its own run evidence. Run 311–315, which the
   reverse read enters in phase, comes first. Its expected gain is small and
   must be measured: standalone 311–315 checks are about 70 of the 665
   samples held only by the HVS-II window on current data, and most checks
   span 303–309. Its value is mainly safety and the mechanism. Run 303–309
   waits for the blind re-check of its confident disagreements
   (results/review-recheck-2026-10-19). Before this increment extends to
   16184–16193 (HVS-I) or 568–573 (HVS-III), research measures the decode
   there.
5. **Read start.** Research on current data (2026-10-20) found no rule that
   decides more at a read's start without error: non-reference calls inside
   a poor start are right far less often than outside, and even clean-looking
   ones there are often one-position shifts. So the outcome is the label. The
   modality asserts the read's poor start as an unreliable span from its own
   data (SRS-CALL-024): to its first run of `start_identity_margin` sound
   calls. The review labels an item covered only by such spans `poor_signal`,
   and a difference there is a check.

## Alternatives considered

- **Clean the whole read before basecalling.** Every call of every read
  would move at once, and a partial solve would cost the whole pipeline. It
  is kept for stage 2.
- **Keep the decoder in `dna-sanger`.** It works, but the user asked for its
  own crate. A crate also tests on synthetic traces and serves another
  dye-terminator modality.
- **Decide from one strand.** Measured in #174 and rejected.
- **Decide by a fitted reporting convention (the "major" rule).** The
  decoder must not encode reviewer habit; the user decides any convention.
- **Name it `dna-signal`.** "Signal" is the run-level product (`dna signal`,
  `dna.signal/v1`, [[PROP-0009]]). The user chose `dna-trace`.

## Validation plan

- **Primary gate.** Unflagged non-het error blocks per zone kind (start, end,
  polyC window, after polyC, interior) and region. They must never rise on
  any period and should fall. Heteroplasmy escapes are reported but not
  gated, and truth-`N` blocks go to the re-check list. Samples with an unreviewed or non-het
  difference must not rise. Clear samples are secondary. The instrument is
  `dna-eval`'s per-zone table (#177): for each kind and region, the zones
  that pass, and those that pass with an unflagged block, het and non-het.
- **Decisions.** Every decided run is reported by bucket: two strands
  agreeing and clear, near-equal, one strand, disagreeing. Each bucket gets
  its count and its accuracy against truth. A bucket is decided only where
  its accuracy, set on development panel v5, holds on readiness-v2 at a
  logged use. The hold-out (validation-v3) stays reserved for `approve`.
- **Flags.** Each flagged run is reported with the components shown, so the
  reviewer workload is measured, not assumed.
- **Unit tests** on synthetic traces for every model step:
  - one component;
  - two and three components;
  - a leakage tail with no mixture;
  - a near-equal pair;
  - no mixture after a run.
- **Increment 1** is byte-identical on the panel.

## Rollout and rollback

One PR per increment. Each raises the `sanger` plugin version and, from
increment 4, the `consensus` version. A revert rolls an increment back.
Increments 2 and 3 change no decision.

## Observability and operations

Logs carry aggregate counts per read: windows re-processed, components,
explained, and decided or flagged runs. The debug dump is opt-in and never
written by default.

## Security and privacy

No new input or dependency outside the workspace. Research outputs stay in
scratch and `results/`.

## Compatibility and migration

`ReadEvidence`, `dna.read_evidence/v1`, `dna.consensus/v1`, and
`dna.review/v1` are revised in place before the first release, as for
[[ADR-0082]]. This proposal takes over [[PROP-0010]]'s increments 2 and 3.
PROP-0010 stays the record of increment 1.

## Decision

Pending.

## Implementation

None yet.
