# `src/variant_calling/mod.rs`

## Purpose

Converts the primary-sequence alignment into normalized, configured-eligible
primary-sequence variants while preserving original trace-call mappings.

## Responsibilities

- Define `call` as the module boundary and sequence extraction before filtering.
- Extract SNVs and ≤50 bp insertions/deletions, exclude unresolved alleles, map
  each variant to its original calls, left-normalize (or circular-canonicalize)
  indels, then apply configured region and supporting-signal eligibility.

## Non-responsibilities

No heteroplasmy inference, breakpoint detection, two-allele decomposition,
allelic fraction fitting, or serialization.

## Key types and functions

- `call(alignment, reference, calls, quality, config) -> Result<VariantCallingResult>`:
  the crate-level entry point.
- Child modules: `extract` (alignment-difference extraction), `mapping`
  (original-call/reference mappings), `normalize` (left/circular normalization),
  and `filter` (configured eligibility).

## Invariants and errors

Variants use validated alleles and explicit coordinates. Equivalent repetitive
indels normalize to the same `(contig, position, ref, alt)` tuple. Unresolved,
oversized, out-of-region, or low-support differences are excluded with one count
per candidate.

## Dependencies

- `config` for `VariantCallingConfig`.
- `model::alignment`, `model::basecalls`, `model::quality`, `model::reference`,
  `model::variant`.
- `error` for `Error`/`Result`.

## Apollo mapping

Primary alignment-difference subset of
`apollo/include/apollo/variant_calling/variant.h`.

## Requirements and decisions

ADR-0001, ADR-0003, ADR-0012; `SRS-VAR-001` through `SRS-VAR-006`.

## Tests

Unit tests in `normalize` cover left normalization, circular canonicalization,
and preservation of observed call mappings. The integration tests exercise variant
extraction through the pipeline.

## Status

Implemented.
