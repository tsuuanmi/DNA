# `src/sample/coverage.rs`

## Purpose

Derives compact run-length reference-coordinate coverage topology from the
selected post-trim segments of independently placed sample reads.

## Responsibilities

- Validate that each read segment is non-empty and that segments within one read
  do not overlap.
- Sweep all segment starts/ends in reference-coordinate order.
- Count total, forward-orientation, and reverse-orientation read depth.
- Emit maximal half-open intervals with constant depth tuples.
- Merge adjacent runs when all three depth values are identical.
- Represent circular origin-spanning reads through their existing two explicit
  linearized reference segments.

## Non-responsibilities

No nucleotide comparison, overlap eligibility, read rejection, quality
weighting, gap scoring, consensus calling, variant adjudication, genotype, or
heteroplasmy inference.

## Invariants

Every emitted interval is non-empty and has
`read_depth = forward_depth + reverse_depth >= 1`. Coverage counts every
independently placed read regardless of pairwise overlap eligibility.
Orientation depth describes selected alignment orientation only; it is not a
claim of biological strand independence.

## Dependencies

`model::sample_evidence`, alignment orientation/segments, deterministic
`BTreeMap` ordering, and typed sample errors.

## Tests

Unit tests cover tiled forward/reverse reads, overlapping depth, adjacent
same-depth run merging, origin-spanning circular segments, missing mapped
segments, and rejection of overlapping segments within one read.

## Status

Implemented.
