# ADR-0030: Add Tracy-derived pairwise overlap admission before sample consensus

- **Status:** Accepted
- **Date:** 2026-09-19

## Context

DNA already processes every trace independently into an immutable
`ReadObservation`, derives orientation and mapped reference segments from the
authoritative alignment, and aggregates arbitrary N-read sample evidence in
reference-coordinate and normalized-variant space.

The next Tracy lesson with high ROI is not consensus voting. Tracy's pairwise
`consensus` command now requires both a minimum aligned overlap and a minimum
primary-character match fraction before combining two traces (defaults 25 bases
and 0.5). Tracy's reference-guided assembly likewise admits independently
processed reads before progressive multi-trace reconciliation.

DNA needs the same explicit pre-consensus boundary without importing Tracy's
pair-first domain model, quality-blind multi-trace majority vote, or
nucleotide-like treatment of gaps.

## Decision

DNA adds deterministic pairwise overlap-admission evidence after every read
has completed independent reference placement.

### Configuration

Strict configuration schema version 5 adds:

~~~toml
[sample_reconciliation]
minimum_comparable_bases = 25
minimum_overlap_agreement = 0.50
~~~

Both settings are required. `minimum_comparable_bases` must be positive and
`minimum_overlap_agreement` must be finite in `(0, 1]`.

### Pair discovery

Reads are first sorted by input SHA-256 under the existing deterministic sample
registry rule. Every unordered pair is evaluated from its selected alignment
columns.

Filename, HV/amplicon label, primer label, declared direction, CLI order, and a
canonical F/R relationship do not participate in pair discovery.

A pair is materialized only when the two reads share at least one reference
coordinate. Non-overlapping reads therefore produce no edge and are not rejected
merely for lacking a partner.

### Nucleotide agreement

For a shared reference coordinate:

- `shared_positions` counts the coordinate whenever both reads cover it;
- the position is `comparable` only when both aligned query symbols are
  canonical A/C/G/T;
- equal canonical bases increment `agreements`;
- unequal canonical bases increment `conflicts`;
- unresolved query symbols and deletions remain shared coverage but do not enter
  the nucleotide agreement denominator.

When at least one comparable base exists:

~~~text
agreement = agreements / comparable_bases
~~~

No synthetic agreement value is created when the comparable denominator is zero.

This deliberately differs from forcing gaps into Tracy's nucleotide match
fraction. DNA already preserves deletion and normalized indel evidence as a
separate evidence type, and future consensus must remain explicitly gap-aware.

### Admission

A pair is eligible for later consensus reconciliation exactly when:

~~~text
comparable_bases >= minimum_comparable_bases
and
agreement >= minimum_overlap_agreement
~~~

Failed rules are retained in deterministic order as:

~~~text
comparable_bases_below_minimum
agreement_below_minimum
~~~

The pair decision does not erase either read, alter placement, mutate a
`ReadObservation`, or change read-level variant eligibility.

### Public evidence

`dna.sample_evidence/v4` adds a compact `overlaps[]` array. Each edge
references the existing reviewer-facing read names and exposes:

- shared reference positions;
- comparable canonical bases;
- agreements and conflicts;
- optional agreement fraction;
- eligibility;
- exact exclusion reasons.

The existing read registry, sparse differential loci, and normalized variant
support keep their semantics.

### No consensus in this decision

This ADR establishes evidence/admission needed by a future sample-consensus
method. It does not emit a consensus sequence, adjudicated sample variant,
majority-vote result, genotype, heteroplasmy estimate, or reference vote.

## Consequences

- DNA adopts Tracy's explicit overlap/agreement gate while keeping the
  N-read, reference-coordinate architecture already implemented.
- Cross-amplicon overlaps such as HV2F/HV3R are treated identically to canonical
  F/R overlaps when their mapped coordinates meet.
- Single reads and disconnected read groups remain representable sample
  evidence; missing a canonical partner is not a failure.
- Gap and indel conflicts remain available to their dedicated evidence layers
  instead of being hidden inside a base-match scalar.
- Strict configuration changes incompatibly from schema v4 to v5.
- Sample evidence changes incompatibly from v3 to v4; no compatibility output or
  alias is retained.
