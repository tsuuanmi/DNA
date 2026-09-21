# ADR-0029: Align reference placement from basecall-independent evidence profiles

- **Status:** Accepted
- **Date:** 2026-09-19

## Context

ADR-0028 introduced one basecall-independent `EvidenceProfile` per PLOC-defined locus.
The existing Gotoh implementation still scored the retained primary sequence as
discrete A/C/G/T/N characters, so the richer evidence could not affect reference
placement or orientation.

Tracy demonstrates the value of profile-aware alignment, but its implementation
accumulates floating-point expected scores and casts them to integers before the DP.
DNA requires an explicit deterministic numeric contract and must not reintroduce
basecall membership as a hidden profile gate.

## Decision

DNA adopts a profile-aware semi-global Gotoh scorer for reference-guided
placement while preserving the existing affine-gap state machine, traceback,
resource bounds, circular-reference handling, and downstream primary-sequence
variant extraction.

### Query evidence

After QC end trimming, the forward query uses the retained `EvidenceProfile`
sequence from the same original call indexes. The reverse query reverses profile
order and complements channel weights:

~~~text
A <-> T
C <-> G
~~~

The retained primary sequence remains attached to traceback columns for call
mapping, public alignment metrics, and primary-sequence variant extraction. It no
longer determines substitution scores when a profile exists.

### Fixed-point substitution score

All DP score deltas use a fixed scale of 1024 integer units.

For canonical reference base `r` and profile support `p(r)`:

~~~text
u = round(p(r) * 1024)          # ties upward because p >= 0
score =
    u          * match_score
  + (1024-u)   * mismatch_score
~~~

The DP therefore stores only `i64` integer scores. It never compares
floating-point values.

A missing evidence profile or non-canonical reference base receives:

~~~text
1024 * ambiguous_score
~~~

All gap deltas are multiplied by the same 1024 scale, so a clean one-hot profile
preserves the relative scoring and path ordering of the previous discrete method.

### Tie semantics

Within Gotoh, existing deterministic predecessor and gap-extension tie ordering is
retained.

Between forward and reverse orientations, only the best fixed-point profile score
selects an orientation. An exact score tie is an explicit alignment error. Primary
basecall exact-match/mismatch counts do not break a profile-score tie.

Within the selected orientation, genuinely distinct equally scoring placements remain
an explicit error. ADR-0047 supersedes the earlier blanket treatment of modulo-distinct
ties: optimal placements that differ only by repeat-equivalent indel shifting are
canonicalized to one deterministic 3'/right-most mtDNA alignment before uniqueness is
assessed.

### Existing admission metrics

The public/internal `callable_columns`, `callable_identity`,
`unresolved_query_bases`, and exact/mismatch counts remain derived from the
retained primary sequence on the selected traceback.

The existing `minimum_callable_bases` and `minimum_identity` gates therefore
remain conservative primary-sequence admission checks after evidence-profile
placement. Changing those public metric semantics requires a separate decision and
contract review.

### Topology and gaps

Affine gaps retain:

~~~text
gap(k) = gap_open_score + k * gap_extension_score
~~~

with every term represented in the same 1024 fixed-point score units.

Circular references continue to align against the doubled reference, reject
tracebacks consuming more than one reference length, and project the selected
coordinates modulo the original reference.

## Consequences

- Mixed and unresolved primary calls can contribute measured nucleotide evidence
  to placement instead of collapsing to one discrete character score.
- Clean one-hot traces preserve the previous Gotoh score ordering, modulo the
  common 1024 scale.
- Alignment score magnitude changes internally and in operational logs; public JSON
  already omits alignment score, so no public schema change is required.
- Orientation ties become strictly evidence-score ties and cannot be broken by the
  primary basecall projection.
- Variant extraction remains primary-sequence based after placement; this ADR does
  not turn mixed signal directly into a variant call.

## Non-goals

This decision does not change basecalling, QC trimming, public alignment metric
definitions, gap-evidence modeling, sample consensus, mixed-indel detection,
allele decomposition, genotype, or heteroplasmy estimation.
