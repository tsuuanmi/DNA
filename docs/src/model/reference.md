# `src/model/reference.rs`

## Purpose

Defines the validated single-record reference sequence and its topology.

## Responsibilities

- Represent the reference name, sequence, topology, and SHA-256 identities.

## Non-responsibilities

No parsing, validation, or alignment.

## Key types and functions

- `ReferenceTopology`: `Linear` or `Circular`.
- `Reference`: name, sequence, topology, and `sequence_sha256`.
- `Reference::len()`: reference length in bases.

## Invariants and errors

- The sequence is non-empty and at most `MAX_REFERENCE_LENGTH` bases (enforced at
  load time).
- `sequence_sha256` identifies the normalized, validated sequence.

## Dependencies

- `serde` for serialization/deserialization.

## Biological semantics

Topology distinguishes linear references from circular ones (e.g. the
mitochondrial rCRS genome). Circular references allow alignments that wrap the
origin and use circular-canonical variant normalization.

## Tests

No dedicated unit tests; behavior is exercised through `reference::fasta` and the
integration tests.

## Status

Implemented.
