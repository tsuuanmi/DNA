---
id: ADR-0083
status: accepted
date: 2026-10-19
family: [consensus]
related: [ADR-0062, ADR-0073, ADR-0074, ADR-0080]
---

# ADR-0083: Set aside dissent from a reference decision that the other strand contradicts

<!-- docs:begin lifecycle -->
- **Status:** Accepted
- **Date:** 2026-10-19
- **Related:** [ADR-0062](0062-read-callability.md), [ADR-0073](0073-sample-consensus.md), [ADR-0074](0074-run-structure-consensus.md), [ADR-0080](0080-reference-only-calls.md)
<!-- docs:end lifecycle -->

<!-- docs:begin backlinks -->
- **Governs:** [design/consensus](../../design/consensus.md), [requirements/consensus](../../requirements/consensus.md)
<!-- docs:end backlinks -->

## Context

A site where reads on opposite strands disagree is a review check
(`strands_disagree`). On development panel v5's current slice there were 305
such items in 257 samples, and DNA's decision matched the reviewers' in 291.
In 12 of the 14 others, the opposing read was the one that matched, so the
check catches real errors and stays as a rule.

Almost every right item has one opposing read against a decision of the
reference, often a contested one-against-one tie. Two mechanisms explain
most of them:
- the opposing call shows a second peak that its modality explains by the
  signal's own shape (a neighbour's tail, a dye blob, a phase shift, an
  off-centre peak, noise);
- the opposing read differs by an insertion or deletion close to its own
  end. Across the corpus, one-read indels that a read of the other strand
  opposes are frequent within about 35–40 aligned positions of a read's
  start or end, then flat at 0.1–0.2 per 10⁴.

Substitutions near a read end were also tried. They resolved a run mix-up
between reads (`2551583:16180`) wrongly, so they are not included.

## Decision

At a site or stretch decided as the reference, the consensus sets aside the
reads opposing it (SRS-CONS-011) when a supporting read on a strand none of
them is on observes every position cleanly, shows no second base standing
out there, and lies far enough from its own ends, and either:

1. **Explained second peak.** Every opposing read shows a second base
   standing out there (`SecondBase.stands_out`), and no second peak there is
   one its modality leaves unexplained (`second_unexplained`). The
   supporting read lies at least `variant_calling.read_end_margin` calls
   from its ends.
2. **Indel at a read end.** Every opposing read differs there only by
   insertions or deletions within `consensus.dissent_edge` calls (40) of its
   own ends. The supporting read lies at least that far from its ends.

The core reads only the kernel's assertions and the reads' informative
intervals, never a modality's labels (ADR-0069). A set-aside read shows
nothing at the site, inside a stretch as well, and the positions are decided
again in the loop that withdraws reference confirmations (ADR-0080). Only
reads against the reference are set aside, so the decided value never
changes: a contested site becomes called, and a called one loses its
dissent.

## Alternatives

- **Set aside a substitution near a read end too.** This was rejected: it
  hides a run mix-up between reads.
- **Measure the distance from the callable span.** The research used it. The
  core reads the informative interval, the one every read-end rule reads.
- **Let the review drop the label instead of the consensus.** This was
  rejected: the decision would still be contested and the dissent recorded,
  so the call and the review would disagree.

## Consequences

Development panel v5 against `main` 9e69f93:

| Period | Unreviewed | Non-het | Clear |
|---|---|---|---|
| current | 440 → 440 | 94 → 94 | 70.45 → 71.00 % |
| 2025-08..11 | 151 → 151 | 51 → 51 | 47.66 → 48.59 % |
| older | 174 → 174 | 84 → 84 | 68.52 → 69.41 % |

- 110 samples become clear.
- `strands_disagree` items on the current slice fall from 271 to 170. Every
  one removed is a site where DNA's decision is the reviewers'.
- No status changes under any measure. The run mix-up samples (`2551583`,
  `2551763`, `2551661`) keep every check.
- A set-aside read leaves the site's `opposing` list for `uninformative` in
  `dna.consensus/v1`.

Plugin `consensus` 13; configuration schema 18 (`consensus.dissent_edge`).
