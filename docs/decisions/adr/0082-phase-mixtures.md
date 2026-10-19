# ADR-0082: A modality may assert a phase mixture, and the core compares its components

- **Status:** Accepted
- **Date:** 2026-10-19
- **Related decisions:** [ADR-0069](0069-plugin-first-modality-core-post-calling.md),
  [ADR-0071](0071-run-length-edits-need-resolved-run-ends.md),
  [ADR-0074](0074-run-structure-consensus.md),
  [ADR-0076](0076-bounded-run-read-end-exemption.md),
  [ADR-0080](0080-reference-only-calls.md)
- **Proposal:** [PROP-0010](../../proposals/0010-length-mixture-decomposition.md),
  increment 1

## Context

After a homopolymer that carries length heteroplasmy, a Sanger read shows the
sum of molecules whose runs differ in length: the same downstream sequence
several times, offset by whole calls. Callability marks the stretch
`dephased`, and a masked call gives no base (ADR-0076), so the consensus
decides the run's length from weaker evidence (ADR-0074) or leaves it
undecided. At the HVS-II C-tract this is the largest block of review work.
The laboratory asks for the region separated into its sequences, each
compared with the rCRS on its own.

The signal holds those sequences. A reference-free fit of the stretch as the
read's own calls plus offset copies recovers them, and needs no fixed cut:
the read's in-phase calls measure its own noise.

## Decision

**Contract** (`dna-kernel`, `ReadEvidence::mixtures`). A modality may assert,
per stretch of calls, a `PhaseMixture`:

- the stretch (`call_start..call_end`) and the calls of the run that opens it
  (`run_calls`);
- its components, each an offset of the read's own calls in whole calls
  (positive: a longer run, lagging; negative: a shorter one) with its share of
  the fitted signal, offset 0 first;
- the residual the components leave, the calls they leave unexplained, and
  whether they explain the stretch (`explained`).

The assertion is generic: any modality that sees a superposition of offset
copies may make it, and the core reads offsets and fractions, never a
modality's labels. It decides nothing. A fraction is a share of the signal,
not a heteroplasmy level, and appears in no calling output (SRS-VAR-009).
The kernel validates the stretch against the read's calls, offset 0 first,
fractions in `[0, 1]`, offsets no longer than the run, and unexplained calls
inside the stretch.

**Sanger** (SRS-CALL-023) asserts a mixture for every dephased segment that a
repeat run opens: the longest run of equal primary calls near the segment's
start, the stretch ending where a run as long starts. A non-negative
least-squares fit adds offsets greedily while each lowers the residual by
more than the read's own in-phase noise. The mixture is explained when its
residual does not exceed the read's in-phase residual at
`signal_processing.decompose_unexplained_quantile`: the copies then fit the
stretch as well as the read's clean calls fit their own bases.

**Core** (`dna-call`, SRS-VAR-022) compares each component of an explained
mixture with the reference on its own: its sequence placed and called as a
read, by the same aligner and record rules, on the read's strand. Its span is
in read order and may cross the origin. An unexplained mixture keeps its
offsets and fractions, with no sequence read from it. No call, record,
consensus, or review output changes. `SampleAnalysis::components` exposes the
result, and `dna-eval` measures it (`components.tsv`).

## Alternatives

- **Gate on the median per-call residual.** More mixtures pass (HV2F 71 %
  against 60 %), but the dominant run is right less often among them (93.1 %
  against 95.6 %), and a long noisy tail passes more often.
- **Use the profile's argmax as the dominant sequence.** On the same panel,
  it showed no gain over the primary calls at either HVS-II run (measured
  before the haplotype metric, so as a comparison only).
- **Decode the dominant sequence from the mixture** (subtract the copies,
  re-pick each call's base, refit). It can restore a base a lagging copy
  hides (`T310` called `C`) but never a call count. It is left as a measured
  follow-up.
- **Fit only segments flagged `after_repeat`.** That flag uses a longer
  repeat than the HVS-II C-tract opens with, and missed it.
- **One breakpoint and two alleles, one of them the reference** (Tracy's
  `decompose`). This was rejected: the reference is not one of the
  sequences, and three-length mixtures exist.

## Consequences

Development panel v5, current slice (batch date from 2026-01), against
`main` 3de20b3:

- `checks.tsv`, `scores.tsv`, and `windows.tsv` are byte-identical.
- Explained mixtures: HV2F 3,644 of 7,594, HV3R 4,135 of 5,806, HV1F 1,390
  of 5,458, HV1R 831 of 3,900. Of the explained ones, 2,337 hold more than
  one component.
- The run an explained component enters in phase in the HVS-II C-tract
  (303–315) is the truth's in 1,857 of 1,943 HV2F reads (95.6 %) and 1,345
  of 1,355 HV3R reads (99.3 %).
- In 275 explained mixtures (2.7 %) no component is placed: without the
  read's profiles, the aligner finds equally scoring placements for it. They
  are logged and keep their offsets and fractions.
- At HVS-I the truth does not record the length of the C run that 16189C
  creates (convention `hvs1_c_length`), and samples with T16189 rarely show a
  mixture there, so HVS-I gives no measure.
- The dominant sequence is the read's primary calls. Where the copies lag
  heavily, a primary reads the lagging copy's base, and the mixture is then
  unexplained or carries the error.

Plugin `sanger` 13; configuration schema 17. `ReadEvidence` is revised in
place before its first release. Increment 2 of PROP-0010 (the consensus
deciding a run's length from the dominant component) is a separate decision.
