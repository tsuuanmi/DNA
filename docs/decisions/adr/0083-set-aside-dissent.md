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
reads against the reference are set aside, apart from the run-length
dissent below, so the decided value never changes: a contested site becomes
called, and a called one loses its dissent.

### Run-length dissent from a called variant

Amended 2026-10-20. On the current slice of development panel v5, 25 samples
were held only by the `start@HV2` zone, each by one check at 73
(`strands_disagree`). In all 25 the reviewers' sequence is `73G` alone, and
the consensus decides the stretch 66-73 (`GGGGGGTA` in the reference) as
`GGGGGGTG`, supported by HV2F and HV3R. The only opposing read is HV1F, the
forward read that crosses the origin. It reads `73G` cleanly as well and
differs only in the length of the G run 66-71: G7 in 24 samples, and in one
(`2565539`) G6 against two reads of G7, a run length outside the reportable
region. HV1F reads the run 630-720 calls into the read, 25-153 calls (median
about 57) before the end of its informative interval and inside its in-phase
signal. In 24 samples HV1F shows a G second peak standing out at 72 that its
modality explains as a neighbour's tail. Rule 1 or 2 did not set it aside
because:

| Reason | Samples |
|---|---|
| The decision is not the reference (`73G` is inside the stretch) | 25 |
| Rule 2 reads HV1F's own `73G` as a substitution | 25 |
| Rule 2: the insertion is within `dissent_edge` calls of HV1F's end | 4 of 25 |
| HV2F cannot anchor: it is on HV1F's strand | 25 |
| Rule 1 would hold at a reference decision | 17 |
| Rule 1 fails: HV3R, the only anchor, shows a G second standing out at 72 | 8 |

Rule 1 therefore also applies at a site or stretch called by a strict
majority (`called`) as a value other than the reference, to opposing reads
whose dissent is a run's length alone: every call they show there is
resolved, and they show the decided runs in the same order. Every other
condition of rule 1 is unchanged. Rule 2 stays reference-only, and a
substitution against a called variant is never set aside. The remaining
supporting reads keep their majority, so the decided value does not change.

The 8 samples where HV3R shows the G second at 72 keep their check. A second
standing out on both strands at one position is the signature of a second
sequence, and the anchor's guard keeps it (`2558783`, `2566026`, `2566574`,
`2566845`, `2567577`, `2567589`, `2568024`, `LN_26_AC7273`).

Measured on development panel v5 against `main` 69e1471:

| Period | Unreviewed | Non-het | Clear |
|---|---|---|---|
| current | 440 → 440 | 94 → 94 | 73.16 → 73.31 % (+17) |
| 2025-08..11 | 151 → 151 | 51 → 51 | 53.38 → 53.41 % (+1) |
| older | 174 → 174 | 84 → 84 | 71.96 → 72.00 % (+1) |

No status of any measure changes, and no unreviewed difference appears or
disappears. Checks change in 21 samples. In 20 the check at 73 goes, which
clears 19 samples: the 17 current samples rule 1 reaches, `LN_25_AA7280`,
and `2525350`. The only other change is a check at 16221-16224
(`2537451`), where DNA's `16223T` is the reviewers' too. Plugin `consensus`
16.

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
