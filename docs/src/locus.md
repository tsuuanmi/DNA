# `src/locus.rs`

## Purpose

Owns validated PLOC-defined locus-window geometry shared by scientific stages.

## Responsibilities

- Build one half-open sample window per validated PLOC locus from neighboring midpoints.
- Apply bounded first/last extrapolation.
- Keep locus geometry independent of basecall classification and signal interpretation.

## Non-responsibilities

No channel peak selection, base calling, signal normalization, quality scoring, alignment, or reporting.

## Invariants

Every window contains its PLOC coordinate, satisfies `start < end`, stays inside the trace sample range, and preserves trace order.

## Status

Implemented.
