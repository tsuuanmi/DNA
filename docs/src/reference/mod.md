# `src/reference/mod.rs`

## Purpose

Loads and validates the reference FASTA and computes its identity.

## Responsibilities

- Re-export `load` as the module boundary.
- Accept one non-empty plain FASTA record, normalize supported DNA symbols, retain
  its identifier, and enforce the direct-alignment length limit.

## Non-responsibilities

No gzip, multi-contig selection, indexing, or alignment.

## Key types and functions

- `load(path, topology) -> Result<Reference>`: the public entry point,
  re-exported from `fasta`.
- Child module: `fasta` (one-record parsing and validation). Sequence identity is
  computed with the shared `checksum` module.

## Invariants and errors

Exactly one record, non-empty sequence, valid supported symbols, and at most
`MAX_REFERENCE_LENGTH` bases. Errors are typed and include path/record context.

## Dependencies

- `config` for `MAX_REFERENCE_LENGTH`.
- `model::reference` for `Reference` and `ReferenceTopology`.
- `error` for `Error`/`Result`.

## Apollo mapping

Focused subset of `apollo/include/apollo/fasta.h` and reference-loading command
logic.

## Requirements and decisions

ADR-0004; `SRS-IN-004`, `SRS-IN-005`, `SRS-IN-006`.

## Tests

Unit tests in `fasta` cover multiple-record rejection. The integration tests
exercise reference loading through the pipeline.

## Status

Implemented.
