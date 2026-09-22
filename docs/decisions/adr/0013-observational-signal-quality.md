# ADR-0013: Add Observational DNA-Quality Windows Before DNA Cleaning

- **Status:** Accepted
- **Date:** 2026-08-24

## Context

DNA retains instrument-analyzed `DATA.9`–`DATA.12` channels and PLOC loci but previously calculated no baseline/noise feature or internal noisy-region annotation. Generic smoothing or hard exclusion can remove genuine secondary peaks, and the repository has no approved truth-labeled corpus that can establish sensitivity.

The decoded channels are evidence tied to the input SHA-256. Replacing them in place would make reported peaks ambiguous and couple preprocessing, calling, QC, and variant behavior.

## Decision

Add `dna.windowed_snr/v1` after basecalling as a pure observation-only stage.

- Preserve the original channels and calls.
- Estimate rolling primary/secondary SNR in sample space over configurable 5–10-call windows.
- Merge overlapping or adjacent candidate-noisy windows into call/sample intervals only when a run contains at least two candidate windows.
- Expose the bounded window records and merged regions in `dna.analysis/v4`.
- Use strict configuration schema version 4, including the required minimum noisy-window run length.
- Do not change trimming, alignment, warning totals, or variant eligibility.
- Do not include disabled smoothing, compatibility aliases, or duplicate v3 output.

## Consequences

Results become larger than v3 but remain bounded by call count and still omit full channel arrays and non-signal per-call tables. Consumers can inspect candidate-noisy regions and variant call indexes without DNA asserting that the regions are uncallable.

The prior configuration schema and output contract are replaced without compatibility paths. Existing historical JSON remains readable as JSON but is not emitted by this version.

## Deferred decision

DNA cleaning, noise-aware base calling, calibrated quality, or variant exclusion requires a later ADR backed by approved truth labels and synthetic minor-peak preservation tests. Any processed signal must be a separate projection; it must never overwrite the decoded evidence.
