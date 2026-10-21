---
id: ADR-0085
status: accepted
date: 2026-10-20
family: [consensus]
related: [ADR-0073, ADR-0078, ADR-0080, ADR-0084, PROP-0012]
---

# ADR-0085: The core decides the analysed range at a region edge from the reads' unreliable spans

## Context

[[PROP-0012]] increment 5 gave every read a poor start: the modality asserts
an unreliable span from the read's first call to its first run of
`sanger_evidence.start_identity_margin` sound calls (SRS-CALL-024). The review
then labelled an item that only such spans cover `poor_signal`, and made a
difference there a check. Research on current data found no rule that decides
more inside a poor start without error.

At the HVS-II read start (positions 73 to 78) this left 602 zones needing a
check on current data. Reviewers already treat that stretch as poor data. On
current data their HVS-II interval starts at 79 rather than 73 in about 40 %
of these samples, and the `poor_signal` label reaches the reviewer's start in
229 of 234.

The user decided on 2026-10-20 to report this stretch as poor data, not
analyse it, and shorten the analysed range to match. Rules hold no target
positions (ADR-0021), so the decision is applied as a rule over the reads'
own evidence, not as "73-78".

## Decision

The core decides the analysed range at a reportable region's end from the
reads' unreliable spans (SRS-CONS-012).

- **Consensus.** The consensus publishes `unreliable_only`. These are the
  positions that some read's call covers inside its unreliable span, and that
  no read shows soundly outside its own span. A sound call is scored, trusted,
  and withheld for nothing on its base: the same test the `poor_signal` label
  uses. The fact is profile-free and is carried in `dna.consensus/v1`, so
  `dna call` of a consensus document reproduces it.
- **Call.** The `decided` call starts at each region's first and last
  position and walks inwards over the contiguous `unreliable_only` positions.
  It leaves out those up to the last one the consensus does not settle as the
  reference. A position is unsettled when it lies in a site that differs from
  the reference, is undecided, or is contested, or when no segment holds it.
  - The positions left out are not callable, and a change there is excluded
    with the span's reason (`poor_start_signal`).
  - The `unreliable_only` positions after the last unsettled one stay
    analysed.
  - Only a region's own ends count. An interior position, or one after a
    stretch no read covers, is never left out.
  - The variants document lists the intervals as `not_analysed`. The EMPOP
    query cuts them with the reason `not_analysed`.
- **Review.** The review lists the stretches as `not_analysed`, with their
  reason, and holds no item or zone over them. It also guards the rest of the
  region's edge zone. Every difference DNA calls within `review.edge_zone`
  positions of a region end it leaves out is a check, labelled
  `beside_not_analysed`. The reads there sit right after their poor start,
  where indels are often wrong.
- **Edge zone width.** The edge zone is the review's existing region-edge
  geometry, `review.edge_zone`. Its default is widened from 20 to 40, the
  width of the review workload's start and end zones. The `edge_deletion`
  label keeps its own zone, `review.deletion_zone`, at its earlier width of
  20 (amended below).

### Variants considered

All variants were measured on development panel v5, current slice, against
the same base.

- **A. Leave out every poor-only position at a region's edge**, skipping
  the positions no segment holds.
  - It cut 16024-16031 at the HVS-I start in 960 samples.
  - The consensus had decided those positions as the reference.
  - Reviewers analyse them: their interval starts at 16024 in 952 of the 960.
- **B. As A, but a position that reference confirmations decide
  (ADR-0080) stays analysed.**
  - It still shortened the HVS-I start in 523 samples.
  - The stretch skipped past positions no segment holds, so it could start
    inside a region. For example, it ended a region at a reverse internal
    primer's poor start.
- **C. As B, with only a region's own ends counting.**
  - It shortened the HVS-I start in 488 samples.
- **D. As C, but leaving out only the positions up to the last one not
  settled as the reference.**
  - It keeps a poor end that reads the reference throughout.
  - The HVS-I start was shortened in 104 samples, the same samples that
    previously needed a check there.
  - The check reduction was A's.
- **E. As D, with the soundness test shared with the review's `poor_signal`
  label (chosen).**
  - In D the consensus and the label disagreed about soundness. One case was
    a position covered at the end of another read under a bounded-run
    exemption.
  - That disagreement left undecided `poor_signal` checks next to a cut
    stretch.
  - With one test there is one implementation: the label reads
    `Consensus::unreliable_only`.

The guard was added after E was measured. On E, three non-het escapes had
appeared in passed start@HVS-II zones in 2025-08..11: an insertion 111.1A in
2544421 and 2547122, and a shift at 102 in 2548396. A check at the cut edge
had covered them before. The guard at `review.edge_zone` 20 left all three
escapes in place. At 40 it catches all three, and no zone's non-het escapes
rise.

### Edge deletions keep their own zone

Amended 2026-10-20. Until then `review.edge_zone` also set the zone of the
`edge_deletion` label: a deletion resting on a single read near a region's
first or last position. Widening it to 40 for `beside_not_analysed` widened
`edge_deletion` too, and its rationale does not carry over. The label exists
because a read's own edge loses or doubles a base. The common true
`106-111DEL` lies 33 positions from HVS-II's start, read by one HVS-II
forward read 40 to 70 calls from its end, and the wider zone made it a
check: 41 items on current data, DNA right in all 41.

`edge_deletion` now reads its own key, `review.deletion_zone`, at its
earlier default of 20, so no new number. `beside_not_analysed` keeps
`review.edge_zone` at 40.

A second variant was measured and rejected. It kept the zone at 40 and also
required the deleting read's flanking calls to lie within
`consensus.dissent_edge` calls of the read's informative end. It added no
key, but a deletion at the HVS-I end (2515360, `16355-16365`) fell outside
it and became an unreviewed non-het escape in a passed `end@HV1` zone: GATE
FAIL. It also had the review read a consensus key, which would let one key
move another plugin's rules.

Measured on development panel v5 against a fresh run of `feat/eval-diff`
ae98063 (`dna-eval diff`, GATE PASS):

| Period | Clear | Unreviewed | start@HV2 checks | interior@HV2 checks |
|---|---|---|---|---|
| current (2026-01 and later) | 75.99 → 76.06 % | 387 → 387 | 285 → 245 | 257 → 221 |
| earlier | 64.55 → 64.62 % | 282 → 282 | 354 → 346 | 348 → 340 |

- **Non-het escapes.** Non-het passed zones are unchanged in every zone and
  every period. No coverage or status of any measure changes.
- **Workload.** `window@HV2` checks fall 1,664 → 1,658. Positions to review
  fall 67,698 → 66,878 on current data. Twelve samples become clear.
- **Het.** One of them, 2526880 (2025), has an unreviewed `16086Y` in an
  `interior@HVS-I` zone that no check touched before either.
- **Configuration.** The configuration schema is 20. Plugin `review` 22.

## Consequences

Measured on development panel v5 against a fresh run of `main` 69ee8be,
with the re-check list `--recheck`. A period is a range of batch dates.

| Period | Clear | Unreviewed samples | With a non-het block | start@HV2 checks | start@HV1 checks |
|---|---|---|---|---|---|
| current (2026-01 and later) | 73.28 → 75.99 % | 437 → 437 | 38 → 38 | 603 → 285 | 431 → 316 |
| 2025-08..11 | 53.56 → 55.46 % | 151 → 148 | 45 → 42 | 607 → 282 | 128 → 98 |
| 2025-12 | 64.53 → 65.04 % | 34 → 34 | 5 → 5 | 38 → 38 | 28 → 16 |
| older | 74.52 → 75.28 % | 140 → 140 | 44 → 44 | 42 → 34 | 91 → 75 |

- **Non-het escapes.** Non-het passed zones are unchanged in every zone and
  every period.
- **Truth-N block.** One truth-N block moves into a passed start@HVS-I zone
  (LN_26_AB0085, 16042N). It is on the re-check list.
- **Positions to review.** On current data, positions to review fall from
  73,956 to 67,698, and zones needing a check from 4,391 to 3,985.
- **Edge-zone cost.** Widening `review.edge_zone` added `edge_deletion`
  checks: 57 → 101 samples on current data. It also caught three
  unreviewed samples in 2025-08..11. The amendment above gives
  `edge_deletion` its own zone at 20 again; on the later base the
  unreviewed samples are unchanged by that.
- **Coverage cost.** A truth allele inside a stretch that is not analysed
  becomes a coverage item: strict-discordant, and ignored by the primary
  measure. It is never dropped silently.
  - 73G: 117 samples on current data, 337 in 2025-08..11, and 9 older.
  - Also lost: 16362C (6 current, 2 in 2025-08..11), 93G, and 103A.
  - DNA had called each of these alleles correctly, in samples that were
    already checks.
  - In 2025-08..11, reviewers analysed from 73 in every sample shortened.
  - The user accepted this cost: the original decision stands.
- **Contracts.**
  - `dna.consensus/v1` gains `unreliable_only`.
  - `dna.variants/v1` gains `not_analysed`.
  - `dna.review/v1` gains `not_analysed`, `summary.positions.not_analysed`,
    and the label `beside_not_analysed`.
  - The EMPOP cut gains `not_analysed`.
  - All are revised in place before the first release.
- **Plugins.** `consensus` 17, `decided` 4, `review` 20, `empop` 2.

## Links

- Rule: SRS-CONS-012. SRS-CALL-024 allows the core this use of the span;
  [[ADR-0086]] adds one more, in placement.
- Increment: [[PROP-0012]] increment 5. Crate: [[ADR-0084]].
