---
id: ADR-0089
status: accepted
date: 2026-10-21
family: [consensus]
related: [ADR-0069, ADR-0071, ADR-0074, ADR-0080, ADR-0083]
---

# ADR-0089: Insertions need both strands

<!-- docs:begin lifecycle -->
- **Status:** Accepted
- **Date:** 2026-10-21
- **Related:** [ADR-0069](0069-plugin-first-modality-core-post-calling.md), [ADR-0071](0071-run-length-edits-need-resolved-run-ends.md), [ADR-0074](0074-run-structure-consensus.md), [ADR-0080](0080-reference-only-calls.md), [ADR-0083](0083-set-aside-dissent.md)
<!-- docs:end lifecycle -->

<!-- docs:begin backlinks -->
- **Governs:** [design/consensus](../../design/consensus.md), [engineering/improvement-method](../../engineering/improvement-method.md), [reference/consensus](../../reference/consensus.md), [reference/review](../../reference/review.md), [requirements/consensus](../../requirements/consensus.md)
<!-- docs:end backlinks -->

## Context

The user decided on 2026-10-21: "If both strands cover a site, an insertion
should be decided by consensus of both strands. If only one strand covers
it, the problem must be identified."

On `main`, one-read insertions in short runs at a read's start passed
without a check and were wrong against the truth: `82.1A` and `93.1A` at the
HVS-II forward start, and `16025.1T` at the HVS-I forward start. These were
the known non-het escapes. Only clean observations decide (SRS-CONS-003), so
a read of the other strand that reads the site through, but not cleanly
(masked but anchoring, or with a modality reason against its bases), was
uninformative. The insertion then rested on one strand, and the review's
`one_read` and `one_strand` labels are context only.

A census on development panel v5 and readiness-v2 classified every decided
site that inserts bases by the other strand's best view of it. The view
types are those of SRS-CONS-015. Current period (batch date ≥ 2026-01),
error rate against the truth where the truth covers the site:

| Where | Other strand | Development v5 | Readiness-v2 |
|---|---|---|---|
| Outside the repeat windows, short repeat | both strands support | 0/7 | 0/1 |
| | shows it, not cleanly | — | 0/1 |
| | reads through without it | — | 1/1 (a check) |
| | ends in the site | 0/2 | — |
| | covers no position | 1/5 | 2/4 |
| Outside the repeat windows, long repeat | both strands support | 0/14 | 0/26 |
| | shows it, not cleanly | 0/0 | 0/2 |
| | ends in the site, shows no base, or covers no position | 0/3 | 1/5 |
| Inside the repeat windows | both strands support | 95/7,774 | 166/12,731 |
| | shows it, not cleanly | 2/1,778 | 3/3,193 |
| | reads through without it | 14/533 | 21/613 |
| | ends in the site, shows no base, or covers no position | 5/535 | 3/734 |

Every wrong insertion outside the repeat windows rested on one strand. Over
all periods, a short-repeat insertion that the other strand read through
without was wrong in every case (`16025.1T` of `2547034`, `16049.1G` of
`2556347`). Inside the repeat windows the picture is reversed. When one
strand reads a C tract through without the insertion, the insertion is
still right in more than 96 % of sites, because that strand reads the tract
out of phase. The review's length labels judge those runs (SRS-API-020).
Of the one-strand insertions inside the windows that pass today without a
check, 3 of 1,037 on current development data are wrong.

## Decision

**Consensus (SRS-CONS-015).** A site whose call is its reference with one
block of bases inserted is a short insertion when:
- the site holds no reference homopolymer of `consensus.seeded_run` bases
  or more;
- the block adds to a tandem repeat of its own unit with fewer than
  `consensus.seeded_run` reference copies.

When such a site is called with all its clean support on one strand, the
consensus looks at every read of the other strand. A read reads the site
through when it covers the site and both neighbouring positions with calls
that each show a base, in phase or masked but anchoring, and does not read
the decided sequence. If some read does that and none reads the decided
sequence, the insertion is not decided:
- the site is decided as the reference, `contested`;
- the reads that read it through support the reference;
- the insertion's reads oppose it.

The review then lists it as a check (`strands_disagree`, `undecided`). A
long run or repeat stays with the run-structure rules (ADR-0074).

The rule reads only the core's trust levels and the reads' strands, never a
modality's reasons, positions, or profile windows (ADR-0069, ADR-0080).

**Review (rules version 29).** A decided insertion outside the profile's
repeat windows and length tracts whose supporting reads all lie on one
strand, where no read of the other strand reads the decided sequence, is a
check: `one_strand_insertion`. The label names the problem, the other
strand's best view by the consensus's own reading:
- `not_covered`: no read of it reaches the site;
- `read_end`: a read of it ends in the site;
- `unresolved`: a read of it shows no base there;
- `reads_through`: it reads the site through without the insertion. This
  happens only where the consensus leaves a long repeat to its run rules.

The label also gives `from_start`, the fewest calls from a supporting read's
first call to the site. The existing context labels give the per-read
reasons (`noisy_signal`, `phase_shift`, `read_edge`, `poor_signal`). Inside
the repeat windows and length tracts, the length labels already judge the
runs, and the label leaves them alone.

## Consequences

Measured against `main` 3ca9580e with the evidence cache and the frozen
re-check list. The `dna-eval diff` gate passed on both panels, with equal
`failed` counts:
- **Development v5, current period:** non-het unreviewed blocks 4 → 4,
  clear 8,071 → 8,066. The wrong `16049.1G` of `LN_26_AB0020` (a re-check
  sample) becomes a check.
- **Development v5, earlier period:** non-het unreviewed 13 → 12. The
  contested `16025.1T` of `2547034` is the artefact bench's known case,
  now fixed.
- **Readiness-v2, current period:** non-het unreviewed 2 → 1 (`16049.1G`
  of `LN_26_AB5540`), clear 14,119 → 14,117.
- **Readiness-v2, earlier period:** non-het unreviewed 14 → 11 (`93.1A`,
  `82.1A`, `16025.1T`).

The consensus rule changes two sites over both panels, both wrong
insertions. The review label makes 15 and 8 sites checks. Plugins: `consensus` 20, `review` 29.

## Alternatives

- **Contest every one-strand insertion the other strand reads through,
  long repeats included.** On development v5 it also contested two right
  insertions where a base splits a C5 run (`312T` with `315.1C`, `314T` with
  `315.1C`). Inside the windows it would contest more than 500 right C-tract
  lengths.
- **Make `one_strand` a check for insertions everywhere.** Inside the
  repeat windows it would add about 1,000 checks on current development
  data for 3 errors.
