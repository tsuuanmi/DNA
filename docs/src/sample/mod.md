# `src/sample/mod.rs`

## Purpose

Defines the sample-level scientific evidence aggregation boundary shared by production and validation tooling.

## Responsibilities

- Expose deterministic production aggregation of independently produced `ReadObservation[]` into compact `SampleEvidence`.
- Expose `CoveredLoci` internally for validation measurements, pairing dense all-covered locus evidence with the exact SHA-ordered read references used by observation indexes so diagnostic provenance resolves without recomputation.
- Keep validation/order, coverage topology, overlap admission, generic locus aggregation, call-signal projection, contribution policy, nucleotide support/geometry, and variant aggregation in focused child modules.

## Non-responsibilities

No input loading, logging, publication, threshold fitting, consensus projection, or metadata-driven pairing.

## Status

Implemented.
