# ADR-0010: Align to the Circular rCRS Reference with Origin Wrapping

- **Status:** Accepted
- **Date:** 2026-08-22

## Context

The MVP reference is the 16,569 bp mitochondrial rCRS sequence, which is
circular: its last base is adjacent to its first. A query that spans the origin
must be alignable across that boundary. ADR-0004 chose direct semi-global Gotoh
alignment against a single FASTA record, but did not specify how circularity is
handled.

## Options

1. Treat the reference as linear and forbid alignments that cross the origin.
2. Duplicate the reference (concatenate it with itself) and align against the
   doubled sequence, then project the result back onto the circular reference.

## Decision

Choose option 2. For a circular reference, the working reference is the sequence
concatenated with itself, and the working length is the modulo length. The
selected alignment is projected back onto the reference: if the aligned span
crosses the origin it is split into two half-open reference segments and
`wraps_origin` is set to `true`. Indel normalization uses the circular-canonical
rotation (`circular_canonical`) rather than linear left-normalization.

## Consequences

Queries spanning the origin align correctly and the result records the wrap
explicitly via `reference_segments` and `wraps_origin`. The doubled reference
increases the alignment cell count, which is bounded by the compiled cell cap.
When the selected alignment crosses the origin, the compact output sets the
boolean `reference_origin_wrap` flag in the warning summary rather than emitting
a free-form info warning.

## Supersession

A future ADR may add reference search or indexing; circular handling must be
preserved behind the same validated reference model.
