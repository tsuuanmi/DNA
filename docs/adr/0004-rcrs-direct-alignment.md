# ADR-0004: Align Directly to Bundled rCRS in the MVP

- **Status:** Superseded by ADR-0010
- **Date:** 2026-08-22

## Context

The MVP reference is the 16,569 bp mitochondrial rCRS sequence. An FM-index and
candidate-region search add persistence, format, and anchoring decisions that are
not required to align a short reference.

## Options

1. Port SDSL-compatible `.fm9` indexing.
2. Adopt the previous Rust `.seq.gz` fast path.
3. Use direct semi-global Gotoh alignment with an explicit reference-size limit.

## Decision

Choose option 3 for one plain FASTA record up to 50,000 bases. Use deterministic
affine-gap Gotoh and bounded allocation. Do not add hardcoded HV rescue positions.

## Consequences

The MVP is simpler and has fewer dependencies. Runtime and memory scale with
query times reference length, so large or multi-contig references remain out of
scope. Indexing can be added later behind the same validated reference model.

## Supersession

A later ADR may add reference search/indexing after benchmarks show direct
alignment is insufficient for an approved use case.
