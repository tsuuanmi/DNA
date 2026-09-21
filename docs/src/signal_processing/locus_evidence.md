# `src/signal_processing/locus_evidence.rs`

## Purpose

Derives one basecall-independent `LocusEvidence` record for every validated PLOC locus.

## Responsibilities

- Derive both the locus window and fixed-width local-statistics context from shared PLOC geometry, without consuming basecall records or rolling noisy-window records.
- Estimate per-channel local baseline and noise.
- Build the total non-negative baseline-corrected A+C+G+T signal series inside the shared locus window.
- Refine the event to the positive total-signal local maximum nearest PLOC; for equally distant candidates prefer greater total signal and then the lower sample coordinate.
- Fall back to the validated PLOC sample only when no positive total-signal local maximum exists.
- Derive raw co-located heights, corrected amplitudes, SNR, and normalized evidence profile directly from chromatogram channels.

## Non-responsibilities

No use of primary calls, ambiguity codes, selected basecall peaks, qualifying-channel membership, alignment, variant filtering, or biological mixture interpretation.

## Tests

Tests cover context placement, unresolved/basecall-independent loci, mixed-channel mass retention, neighboring-event rejection, equal-distance event tie handling, PLOC fallback, and zero-signal loci.

## Status

Implemented.
