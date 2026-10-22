---
implements: [SRS-API, SRS-CONS]
decided-by: [ADR-0074]
covers: [src/report/length.rs, crates/dna-normalize/src/length_convention.rs]
---

# Review length labels

The run-length labels of [[reference/review|`dna.review/v1`]]: a run's
estimated length, the reads' tract decodes of an anchored repeat window
(SRS-API-020, [[PROP-0012]] increment 3), and the noted length of a profile
length tract. Their ranks and the item kinds they make are in the review
reference.

## `estimated_length`

`{runs: [{base, length, length_evidence}]}`: the site's runs whose length is
`phase_loss` or `reference_frame`, so their 3' end is not callable. A check.

- A run whose length a profile report convention leaves unrecorded is not
  listed (rules version 9): a `run_length` convention labelled
  `not_recorded` whose window the site overlaps, of the run's base, with the
  run at least its `minimum_run` long without the bases the site's call
  adds to the reference. An estimate lengthens a run only by adding bases,
  so a run that reaches the bound only by them may be one whose length the
  reviewers record, and stays listed. On the human profile this is the
  HVS-I C run of ten or more that 16189C creates.
- A run of a profile length tract moves to `tract_length` (rules version
  28, below): DNA reports no length there, so the note shows the reads'
  floor reading, unconfirmed, never the estimate.

## `length_unconfirmed`

`{runs: [{run, base, called, agreement, convention, floor, reads}]}` (rules
version 21): runs of an anchored repeat window (on the human profile, HVS-II
`C7 T C5`) that the reads' tract decodes reach but do not confirm. A check.

`run` is 1 before the anchor and 2 after it (1 alone when DNA calls one
run). `called` is DNA's length, null when DNA reports none: its call holds an
`N`, or (rules version 30) a position DNA must show to measure the run is
not callable, so its decided length is an estimate, kept only in
`estimated_length`. When DNA's calls in the window only insert or delete
the repeat base, none at the anchor, run 1 is measured from the position
before the window to the anchor and run 2 from the anchor to the position
after the window. Where the consensus shows the anchor (a consensus
segment holds it and more than reference confirmations observe it), the
anchor is DNA's decided base the runs lie around and need not be callable
itself: run 1 then ends before it and run 2 starts after it (rules version
35). An insertion excluded at run 1's estimated 3' end leaves both
positions beside its junction, the anchor among them, out of the callable
segments ([[requirements/variants|excluded variants]]), and before rules
version 35 that left run 2 unreported too. An anchor the consensus does not
show must be callable. Otherwise (one run, or a call that changes another
base of the window, such as `T310C`) every run needs every window position
callable, as rules version 27 required of all runs. Each read gives a length with its `source`:

- `decode`, the first-run marginal of the run the read enters first (a
  forward read enters run 1, a reverse read run 2);
- `derived`, the next run, from the read's total, when its mark is the anchor
  and both its marginals are clear;
- `in_phase`, a run the consensus reads in phase in that read (ADR-0074).

Each also carries its `share`, `next` length, and `margin`. `agreement` says
why the run is not confirmed, the first that applies: `structure_differs` (a
read sees two runs where DNA calls one, or the reverse), `disagree` (the
strands give different clear lengths), `not_clear` (a strand's lengths are
not all clear and equal), `one_strand` (one strand gives no length, and no
read of the other shows the run in phase), `undecided_call` (the reads give
one length, DNA none), `differs_from_call` (the reads agree on another length
than DNA's).

The reads give a run its length when both strands give it one clear length,
the same, or (rules version 26) when one strand gives none and one of the
other's reads shows its one clear length `in_phase`; decodes alone on one
strand need the other. A run whose length so read is DNA's reported length is
confirmed and not listed, as is a run the consensus decided by the
two-strand rule (`two_strand`, SRS-CONS-013, rules version 23) unless its
structure differs.

`convention` (rules version 28; rules version 24 read a near-equal pair) is
the length convention's reading, a suggestion and never DNA's call:
`{length, rule: "shortest_dominant"}`. With the `[normalization]` policy
configured, each read gives a reading, its in-phase length or the policy's
reading of its decode (the next-shorter length when its share is at least
`shortest_dominant_ratio` of the dominant's, ISFG 2014 3.3.2). When every
read gives one and they hold one adjacent pair `n`, `n + 1`, with some read
topping elsewhere than `n`, the run reads `n`. Null otherwise, without the
policy, for a run whose structure differs from DNA's, and in
`length_not_recorded`.

`floor` (rules version 27): for a run whose `called` is null, the least
length any read gives as its top or `next`, a suggestion only like
`convention`; null otherwise and in `length_not_recorded`.

## `length_not_recorded`

`{runs: [...]}` (rules version 21): the same runs, as `length_unconfirmed`
lists them, confirmed ones included, of an anchored window over which the
profile's report convention leaves the run lengths unrecorded (on the human
profile, HVS-I). Context only: the reviewers do not record those lengths,
and a check against them would follow a convention, not an error.

## `tract_length`

`{tract, ignored_in_search, runs: [{run, base, called, status, reading,
agreement, strands, in_phase, floor, convention, reads}]}` (rules version 28,
user decision 2026-10-21): runs of a profile length tract (`length_tracts`,
[[reference/profiles|profiles]]) whose length is length heteroplasmy, noted.
Context, never a check (ISFG, Parson et al. 2014, Rec. 10, 11, 16; SWGDAM
2024 5.1.1). `ignored_in_search` is the tract's search window, where EMPOP's
standard query disregards its length variants (the EMPOP section lists them,
[[reference/notation|notation]]): EMPOP's search behaviour, not DNA's
confidence in a call.

A run of the tract is noted when its dominant length is determined. `status`
says how, and `reading` is the length a reviewer would record:

- `dominant`: every read tops at DNA's reported length, and the run's
  structure is shown, by a read in phase or by decodes on both strands;
- `convention`: the reads hold one adjacent pair `n`, `n + 1` that the
  length convention reads as `n` (`reading`), DNA's length one of the pair,
  the structure shown as for `dominant`;
- `confirmed`: the reads confirm DNA's length (listed only for an
  `estimated_length` or `long_run` run the note takes over);
- `unconfirmed`: DNA reports no length (a position DNA must show to measure
  the run is not callable);
  `reading` is the reads' floor reading, never DNA's estimate.

Each run also gives the reads' `agreement` (as `length_unconfirmed` gives
it), each strand's one clear length (`strands`), the reads `in_phase`,
`floor`, `convention`, and every read's decode components: its first-run
and total length distributions (`lengths`, signal `shares`, never
heteroplasmy levels, and the margin `z`).

A run stays a check (`length_unconfirmed`, or its `estimated_length` or
`long_run`) when:

- its dominant length is undecided: reads top at lengths no convention
  resolves, among them an in-phase read of another length, or a near tie
  without a configured policy;
- the reads' structure differs from DNA's call (the anchor base read
  otherwise, as with `T310C`, or a mark where `16189C` leaves one run);
- its only views are decodes on one strand: a decode's structure choice
  adds a spurious mark or drops one in 6-13 % of simulated reads, and its
  margin is over-confident under noise (trace-simulator benchmark,
  2026-10-21).

The tract holds only the homopolymer count. On the human profile HVS-II's
tract is its first run (303-309), so `315.1C`, a phylogenetically stable
insertion, and every composition change (`T310C`, `16189C`) are called and
checked as usual. Where the profile records no length (HVS-I, the
`hvs1_c_length` convention), the note gives no `reading`, `floor`, or
`convention`: it notes the heteroplasmy and asserts no `16193.x` count. A
plain tract outside every anchored window (HVS-III) is judged as one run and
noted only where it takes over an `estimated_length` or `long_run` run.
