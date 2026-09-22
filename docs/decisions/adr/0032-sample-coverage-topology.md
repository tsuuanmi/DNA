# ADR-0032: Expose sample coverage topology before consensus

- **Status:** Accepted
- **Date:** 2026-09-19

## Context

DNA already processes every AB1 independently, derives its selected reference
placement, retains a deterministic sample read registry, and exposes Tracy-derived
pairwise overlap/admission evidence.

The next reference-guided consensus requirements need two pieces of evidence that
are not yet explicit at sample scope:

1. the local number of independently placed reads covering each reference
   coordinate;
2. the forward/reverse orientation composition of that local coverage.

Jumping directly from the overlap graph to a consensus vote would force a new
quality/admission policy before DNA has calibrated one. Tracy's multi-trace
assembly ultimately uses quality-blind character voting; DNA explicitly does
not adopt that rule.

## Decision

DNA adds a compact run-length encoded `coverage[]` topology to sample
evidence before implementing consensus.

### Source of coverage

Coverage derives only from the selected post-trim `reference_segments` already
stored for every independently placed read.

Filename, HV/amplicon label, primer label, CLI order, pairwise overlap
eligibility, nucleotide agreement, variant eligibility, and optional vendor
metadata do not change whether a selected read covers a reference coordinate.

A deletion remains reference-coordinate coverage because the read alignment
spans that reference position. An insertion has no additional reference
coordinate and therefore does not increase coverage outside the mapped segments.

### Coverage counts

For each reference coordinate, define:

~~~text
read_depth     = number of selected reads whose mapped segments cover the coordinate
forward_depth  = covering reads with forward selected orientation
reverse_depth  = covering reads with reverse selected orientation
~~~

Every covered coordinate satisfies:

~~~text
read_depth = forward_depth + reverse_depth
read_depth >= 1
~~~

These values are coverage topology only. They are not canonical-base support,
agreement, confidence, or consensus admission.

### Run-length encoding

DNA emits maximal 0-based half-open reference intervals over which all three
depth values are constant.

Adjacent intervals with identical depth tuples are merged even if the identity
of the covering read changes at the boundary. Read identity and exact per-read
segments remain available in `reads[]`; `coverage[]` intentionally summarizes
only the local denominator/topology.

For circular reads that cross the origin, the existing two reference segments
contribute independently near the end and beginning of the linearized reference.
Coverage runs never wrap implicitly across the JSON coordinate boundary.

### Relationship to overlap admission

`coverage[]` counts **all independently placed reads**. An ineligible pairwise
overlap edge does not remove either read from coverage.

This preserves high-conflict and disconnected evidence instead of silently
converting pairwise reconciliation policy into read rejection.

`overlaps[]` remains the pairwise admission layer and retains canonical-base
agreement/comparability evidence separately.

### Public contract

`dna.sample_evidence/v5` is replaced by
`dna.sample_evidence/v6`, adding required top-level `coverage[]`.

Each item contains:

~~~text
reference: { start, end }
read_depth
forward_depth
reverse_depth
~~~

Earlier v5 schema/example files are removed. No compatibility alias or duplicate
output is retained.

## Consequences

- Future consensus code has an explicit local coverage denominator without
  reconstructing it from read segments.
- Forward/reverse support topology is visible without claiming strand
  independence or biological confirmation.
- Sparse differential-locus and normalized-variant evidence keep their existing
  semantics.
- No quality weighting, majority vote, consensus sequence, sample variant
  verdict, genotype, or heteroplasmy estimate is introduced.
- A later consensus ADR must define local contributor eligibility and
  nucleotide/gap decision rules separately.
