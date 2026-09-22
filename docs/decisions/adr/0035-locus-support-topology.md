# ADR-0035: Factor differential-locus support topology before consensus

- **Status:** Accepted
- **Date:** 2026-09-19

## Context

ADR-0032 exposes local coverage depth, ADR-0033 factors normalized-variant support,
and ADR-0034 preserves reference-oriented nucleotide profiles at sample scope.

Sparse differential loci still require consumers to rescan observations to learn
the local evidence shape: how many covering reads are forward/reverse and how many
observe reference, alternate, unresolved, or deletion state.

That topology is useful pre-consensus evidence, but it is not itself a vote,
confidence score, or biological verdict.

## Decision

Every retained `LocusDifferenceEvidence` derives one `LocusSupportTopology` from
its authoritative observations and selected read orientations.

The topology records:

~~~text
reads
forward_reads
reverse_reads
reference_reads
alternate_reads
unresolved_reads
deletion_reads
~~~

with invariants:

~~~text
reads = forward_reads + reverse_reads
reads = reference_reads + alternate_reads + unresolved_reads + deletion_reads
~~~

Each read can contribute at most one observation to one retained reference
coordinate, as already enforced by sparse locus aggregation.

This topology was initially internal and observable through aggregate production logging. ADR-0053 later promotes the same topology into each public `dna.sample_evidence/v8` differential-locus record without changing its semantics.

## Consequences

- Future sample interpretation gets an explicit local denominator and conflict
  shape without rescanning raw read mappings.
- Reference support remains distinct from missing coverage.
- Unresolved and deletion evidence remain first-class rather than being folded
  into nucleotide support.
- Forward/reverse counts remain derived alignment dimensions and do not imply
  biological independence.
- No threshold, weighting, consensus base, sample-variant verdict, genotype, or
  heteroplasmy policy is introduced.

## Non-goals

This decision itself does not define local contribution eligibility, evidence weights, gap quality, consensus state, sample confidence, or a public schema change; ADR-0053 separately defines the later public projection.
