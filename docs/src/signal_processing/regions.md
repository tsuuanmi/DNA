# `src/signal_processing/regions.rs`

## Purpose

Coalesces rolling candidate-noisy windows into concise regions.

## Responsibilities

- Consume ordered signal windows.
- Union overlapping or directly adjacent candidate-noisy call/sample intervals
  only when the consecutive run contains the configured minimum number of windows.
- Suppress isolated candidate windows from merged noisy regions.
- Preserve the minimum primary SNR across each merged region.
- Leave regions separated across every clean gap.

## Non-responsibilities

No feature calculation, thresholding, gap filling, base calling, trimming, or variant eligibility.

## Key function

- `merge(windows, minimum_noisy_windows) -> Vec<NoisyRegion>`.

## Tests

Unit tests cover overlapping, adjacent, and clean-gap-separated windows.

## Status

Implemented.
