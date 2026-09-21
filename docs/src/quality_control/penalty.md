# `src/quality_control/penalty.rs`

## Purpose

Computes per-call ambiguity and peak-spacing penalties and the best contiguous
section.

## Responsibilities

- Compute a penalty per call from local ambiguity and peak-spacing deviation.
- Find the contiguous section of the configured best fraction with the minimum
  penalty sum.

## Non-responsibilities

No quality-score conversion or trimming.

## Key types and functions

- `PenaltyResult`: penalties, best-section bounds, and best average.
- `calculate(calls, window_size, best_fraction) -> Result<PenaltyResult>`: the
  entry point.

## Invariants and errors

- Calls must be non-empty and `window_size` positive; otherwise
  `Error::QualityControl`.
- The best section length is `max(1, floor(count * best_fraction))`, clamped to
  the call count.
- Penalty arithmetic is saturating; ambiguity-count overflow returns
  `Error::QualityControl`.

## Dependencies

- `model::basecalls` for `BaseCalls`.
- `error` for `Error`/`Result`.

## Biological semantics

Ambiguous calls and irregular peak spacing indicate lower-confidence regions.
The best section is the lowest-penalty contiguous run, used as the reference
baseline for end trimming.

## Tests

- `calculates_deterministic_ambiguity_and_spacing_penalties`: verifies penalties
  and best-span selection with the compact `BaseCalls` fixture, where ambiguity
  is stored per call and only the primary sequence is aggregated.

## Status

Implemented.
