# `src/pipeline/sample_reads.rs`

## Purpose

Runs the authoritative per-trace scientific path for multi-trace sample and validation operations.

## Responsibilities

- Process each loaded trace through `pipeline::observation`.
- Reuse one operation logger for nested basecall/signal/QC/alignment/variant stages.
- Preserve sample read start/completion markers and aggregate warning totals.
- Return completed `ReadObservation` values without aggregation or publication.

## Non-responsibilities

No sample reconciliation, validation truth, thresholding, reporting, or filesystem publication.

## Status

Implemented.
