# `src/report/json.rs`

## Purpose

Assembles compact `dna.analysis/v7` and provides deterministic JSON serialization.

## Responsibilities

- Validate completed read/reference identity consistency.
- Project read/trim, trace-integrity/signal-quality, selected post-trim alignment, normalized variants, and warnings.
- Pass the selected alignment orientation into variant-call projection so public base/peak evidence is reference-oriented.
- Serialize typed results deterministically with a trailing newline.

## Non-responsibilities

No input loading, scientific stage execution,  compatibility output, or filesystem publication.

## Invariants

- `schema_version` is `dna.analysis/v7`.
- Single-read analysis continues to omit the input filename.
- `reference_segments` originate from alignment of the retained post-trim sequence.

## Status

Implemented.
