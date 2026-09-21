# `src/signal_processing/integrity.rs`

## Purpose

Derives observation-only whole-trace integrity evidence from validated PLOC
locations, analyzed signed-16-bit channels, and already-derived `LocusEvidence`.

## Responsibilities

- Require one locus-evidence record per PLOC locus.
- Retain PLOC count and optional vendor PBAS/PCON counts.
- Summarize adjacent PLOC spacing as minimum, median, and maximum when at least
  two loci exist.
- Count analyzed channel samples exactly equal to `i16::MIN` or `i16::MAX`.
- Sum corrected A/C/G/T amplitudes at each refined locus event and emit the
  six-decimal maximum-to-median ratio across positive event totals.

## Non-responsibilities

No ABIF parsing, PLOC fallback, event discovery beyond the current PLOC-defined
method, artifact/dye-blob classification, thresholding, base calling, trimming,
alignment, variant filtering, or evidence mutation.

## Invariants

A vendor-series length mismatch is evidence, not a reason to add/remove loci.
Clipping is an exact representation-boundary observation. The event-signal ratio
is absent when no locus has positive corrected signal and is never treated as an
artifact probability or automatic filter.

## Dependencies

`model::trace`, `model::locus_evidence`, `model::signal`, shared signal
statistics, and typed errors.

## Tests

Unit tests cover opposite-direction vendor/PLOC cardinality mismatches, PLOC
spacing, exact positive/negative clipping, deterministic event-signal ratio, and
absence of artifact reclassification.

## Status

Implemented.
