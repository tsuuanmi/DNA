# Overlaps and Differential Loci

Part of the canonical [sample evidence contract](../sample-evidence.md).

## Pairwise overlap admission

`overlaps[]` is the **pre-consensus reconciliation layer** learned from Tracy's
explicit minimum-overlap and minimum-match admission checks, adapted to DNA's
N-read reference-coordinate model.

DNA evaluates every unordered pair only after both reads have independently
completed placement. A pair is present when the two selected alignments share at
least one reference coordinate. Non-overlapping reads have no edge; this does not
reject either read and no canonical F/R partner is required.

Each overlap record contains:

- `left` / `right`: reviewer-facing read names from the SHA-sorted registry;
- `shared_positions`: reference coordinates covered by both reads;
- `comparable_bases`: shared coordinates where both aligned query symbols are
  canonical A/C/G/T;
- `agreements` / `conflicts`: equal versus unequal canonical base/base
  observations;
- optional `agreement = agreements / comparable_bases`;
- `eligible`: whether the edge meets the configured pre-consensus gates;
- `exclusion_reasons`: exact failed rules.

The configured defaults are Tracy-derived:

~~~text
minimum_comparable_bases = 25
minimum_overlap_agreement = 0.50
~~~

DNA intentionally does **not** copy a base-vs-gap scalar match rule.
Deletions and unresolved symbols remain part of shared coverage but do not enter
the nucleotide agreement denominator. Insertions and deletions remain explicit in
the existing locus/normalized-variant evidence so future consensus can treat gap
support as a separate evidence problem.

An eligible edge has no exclusion reasons. An ineligible edge records
`comparable_bases_below_minimum`, `agreement_below_minimum`, or both. This decision
does not erase a read, mutate placement, or change read-level variant eligibility.

## Sparse locus differences

`locus_differences[]` is the **alignment-observation layer**. It contains only
reference positions where at least one covering read is not a canonical reference
match. It answers:

> What did each covering read actually observe at this reference coordinate?

Each retained locus also exposes `support_topology` derived from the same
authoritative observations:

- total and forward/reverse read counts;
- reference/alternate/unresolved/deletion read counts;
- profile-bearing read count and its forward/reverse partition.

These counts are evidence topology only. They are not a vote, confidence score,
or claim of biological independence.

At a retained locus, `observations[]` contains every read that covers that locus,
including reads that agree with the reference. A call-backed observation contains:

- `read`: human-readable read name;
- `state`: `reference`, `alternate`, or `unresolved`;
- `base`: reference-oriented observed base;
- `quality`: the existing uncalibrated relative quality score;
- optional `profile`: normalized basecall-independent A/C/G/T evidence already
  projected to reference orientation;
- `in_noisy_region`: whether the source call falls inside the existing merged
  candidate-noisy call context.

A valid zero-signal call has no `profile`; DNA does not invent a called-base,
uniform, or reference-derived replacement. A deletion observation contains only
`read` and `state: "deletion"`; DNA does not fabricate deleted-base
quality, profile, or noisy-call context.

Dense all-reference positions are omitted. The compact default is explicit:

- inside a read's mapped `reference_segments`, absence from
  `locus_differences[]` means that read is a canonical reference match there;
- outside the read's mapped segments, the position is uncovered;
- at a retained differential locus, the explicit observations are authoritative.

This preserves reference support versus missing coverage without serializing
routine reference loci one by one.

Inserted query bases have no reference-coordinate locus and therefore do not
create a synthetic `locus_differences` entry.
