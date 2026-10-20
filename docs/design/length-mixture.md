---
implements: [SRS-API, SRS-CALL, SRS-VAR]
decided-by: [ADR-0074, ADR-0082, ADR-0084]
---

# Length mixtures

A homopolymer that carries length heteroplasmy is a population of molecules
whose runs differ in length. After it, a Sanger read shows their sum: the
same downstream sequence several times, each copy offset by a difference in
run length. Callability marks that stretch `dephased` and masks its calls
(SRS-CALL-004). `dna-trace`'s `decompose` module ([[ADR-0084]]), which the
Sanger plugin calls, asserts what the stretch
holds as a phase mixture (SRS-CALL-023, ADR-0082), and `dna-call`'s
`components` module compares each of its sequences with the reference on its
own (SRS-VAR-022). Neither decides a call.

## Model

Over a stretch of calls `i`, with `y_i` the call's normalised A/C/G/T profile
and `S` the read's primary calls,

`y_i = Σ_d f_d · onehot(S[i − d])`, `f_d ≥ 0`.

Offset 0 is the read's dominant sequence. A positive offset is a copy that
lags it, the molecules whose run is `d` calls longer; a negative one leads it,
a shorter run. No offset is longer than the run that opens the stretch.

## Stretch

Every `dephased` segment of the read is fitted when a repeat run opens it:
the longest maximal run of two or more equal canonical primary calls whose end
lies within `callability.repeat_window_calls` calls of the segment's start,
the nearer on a tie, and which starts at or after the trim interval's start:
a dephased stretch in the trimmed start is start noise, too short to place. The longest run nearby is the one whose length varies
most; a nearer short run (the `AAA` before the HVS-II C-tract) would bound
the offsets too tightly. The stretch starts at that run and ends at the
segment's end, or where the next run of primary calls at least as long as
the opening one starts, a new source of mixture. Shorter runs inside it (the
`C6` after `T310`) are part of the tract the offsets explain.

## Fit

The fractions are a non-negative least-squares fit (Lawson–Hanson) of the
stretch's profiles. The offsets are chosen greedily: from offset 0 alone, the
offset that lowers the residual most is added while it lowers the residual by
more than the read's own noise, the root mean square of
`profile − onehot(primary)` over the read's in-phase single-channel calls. The
read sets every level, so nothing is a fixed cut.

- A call of the dephased segment whose own residual exceeds the
  `signal_processing.decompose_unexplained_quantile` (0.99) nearest-rank
  quantile of the in-phase residuals is unexplained: no sum of the copies
  accounts for it.
- The mixture is explained when its residual does not exceed that quantile:
  the copies fit the stretch as well as the read's in-phase calls fit their
  own bases. Otherwise the components are only the fit's best, and no
  sequence is read from them. A long noisy tail, where the primaries are
  mixture noise, fails here.
- A read with fewer in-phase single-channel calls than
  `callability.minimum_callable_calls` asserts no mixture.

## Components

Each component of an explained mixture is a sequence of its own: the read's
calls over its informative interval, with the opening run lengthened by its
own last `d` calls (or shortened by them), since the molecule whose run
differs carries the same sequence on either side. The whole read places it,
so a mixture near the read's start places as reliably as the read. It is
placed by
the same aligner as a read, on the strand the read was placed on, with
one-hot profiles of its bases, and called by the same record rules, so its
edits are its own. The aligner may split one sequence into different edits
(`C9TC6` as `309 C→CCCT` with `310 T→C`), so a consumer compares haplotypes,
not edit lists.

The component's span covers the tract and the stretch: the opening run, the
runs of its base one other call apart before it in read order, and the calls
after it. Read from its 3' end (HVS-III reverse), the HVS-II C-tract opens
its mixture with `C9`, the longest run, after `C6 T`, and the `C6` run is the
one that read enters in phase. The span is the first and last reference
positions the tract and stretch place on, in read order on the reference
strand; the first exceeds the last
when the stretch crosses the origin (a long HVS-I forward or HVS-III reverse
read). Its edits are those within the span and on the bases either side of
it, where an indel at the span's first base is anchored and the runs it
measures end. A component that cannot be placed is left out and logged.

## Limits

The dominant sequence is the read's primary calls. Where the copies lag
heavily, a primary can read the lagging copy's base: the `T` at 310 called
`C` (`C16` for `C9TC6`), or one more `C` at the exit of a long run. The fit
then fails the explained test or carries the error into its dominant
component. Decoding the dominant sequence from the mixture, subtracting the
copies at each call, could restore a base but never a call count; it is a
measured follow-up, not part of this design.

## Tract decode

The offset model above reads each dephased stretch as the read's own calls
shifted, so it inherits the read's primary calls where they are wrong. The
tract decode (SRS-CALL-025, [[PROP-0012]] increment 2) instead fits the
molecules themselves. It reads the analysed channels at each call's locus
sample, divided by the row sum, and no reference.

- **Where.** Before each call where callability sees the read leave phase
  (any segment not `in_phase`), the tract is the latest run of at least
  `trace.locate_min_run` equal primaries starting before it, extended back
  over a single-call mark when at least `trace.locate_min_before` run-base
  calls precede that mark (`CCCCCTCCCC` is one tract). There is no fallback:
  an onset without such a run yields no tract.
- **Model.** Over a window of calls from the tract's first, each molecule is
  `run^a mark run^b D`, or `run^L D` without a mark, with `D` the downstream
  sequence every molecule shares and the decode does not know. The fractions
  of all molecules over the candidate lengths, from `trace.span_below` under
  the run extent to `trace.span_above` over it, are fitted by NNLS. Given
  them, `D` is decoded exactly by dynamic programming over at most
  `trace.joint_lengths` consecutive lengths holding the most signal. The two
  steps alternate from a start per candidate length, at most
  `trace.iterations` rounds each, and the start with the lowest residual is
  kept. A read that ends within the candidates fits only the lengths its
  window holds: a longer molecule would show its run alone, like every other
  one past the read's end. A molecule's column is its total length's run plus its mark's row,
  so only two marginals are identifiable, the first run's and the total's,
  and only they are published.
- **Structure.** A single run and a mark of each other base are all fitted
  over the same calls: the longest candidate extent's window. The lowest
  residual wins, the earlier candidate on a tie, and every candidate's
  residual is published. On the research spike, no two came within a
  relative 2e-5.
- **Margins.** Each marginal gives its lengths' shares with standard errors
  from the fit. The margin of the dominant length over the next is
  deflated by the residual's lag-1 autocorrelation ρ, by
  √((1+ρ)/(1−ρ)), and below `trace.near_equal_z` the two lengths are
  near-equal. When two lengths show and the deflated margin is not finite
  (ρ ≥ 1, or a zero standard error), the margin is `undefined`, which is
  never clear: only a clear margin on both strands may later decide.
- **Cleaning.** The read's own base signatures (each base's mean profile
  over its in-phase calls before the tract) and a leakage null (the share a
  neighbour one or two calls away adds, measured away from runs) refit the
  chosen structure. That cleans the minor shares and does not change the
  dominant length.

The decode is evidence only, a second assertion beside `PhaseMixture`
([[ADR-0084]]). The steps are memoised by downstream sequence, and the NNLS
reuses its QR prefix, so the results equal a plain implementation's. The
research reference's vectors check every output.

## Flag (increment 3)

The review judges the runs of each anchored repeat window from the decodes
(SRS-API-020). It credits each decode to the run its first call is aligned
to, the run at the tract's start. A forward read gives run 1, the one it
enters first, and a reverse read gives run 2. With the anchor as its mark
and both marginals clear, a read also gives the next run's length from its
total.

A read that stays in phase through a run is ADR-0074's `in_phase` evidence.
It gives that run's length as well, read by the consensus's own
`read_runs_over`, so the review and the consensus share one meaning of "in
phase". Without that, a reverse read in phase through the C5 would leave
most runs to one strand.

A strand gives a run one clear length when all its reads agree with clear
margins. A run is confirmed when both strands give the same length and it
is DNA's call, or when one strand gives none and the other's one clear
length, DNA's call, is shown in phase by at least one of its reads. An
in-phase read shows the length directly; the other strand having no view
is absence of evidence, not a contradiction. A decode is a model of a
dephased tract, so a run that only decodes give still needs both strands.
Every other run a decode reaches is a check
(`length_unconfirmed`). Where the profile leaves the window's run lengths
unrecorded (HVS-I), the same view is context (`length_not_recorded`).

Increment 4 replaces "clear" with the data-derived decision rule, one
function the consensus uses to decide and the review uses to confirm.

## Evaluation

`dna-eval run` writes `components.tsv`: per component, the run clusters of the
reference it covers, measured as haplotypes against the truth (the run a read
enters in phase, and whether the whole cluster agrees), and `tracts.tsv`:
per tract decode, the reference position of its first call, its structure,
and its marginals' dominant lengths, margins, and margin states. See the
[evaluation harness](../validation/evaluation-harness.md#outputs).
