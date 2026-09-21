# `src/quality_control/quality.rs`

## Purpose

Converts non-negative penalties into deterministic, uncalibrated bounded quality
scores.

## Responsibilities

- Map a penalty vector to a bounded score in `0..=maximum`, with the lowest
  penalty receiving the maximum score.

## Non-responsibilities

No penalty computation, trimming, or calibration.

## Key types and functions

- `relative_scores(penalties, maximum) -> Vec<u8>`: the conversion.

## Invariants and errors

- If the maximum penalty is non-positive, every call receives the maximum score.
- Otherwise each score is `floor(maximum * (1 - penalty / max_penalty))`, with the
  fraction manually clamped to `[0, 1]` so the resulting score stays in
  `[0, maximum]`.
- The function is total and never fails.

## Dependencies

None.

## Biological semantics

Scores are relative and uncalibrated: they rank calls within a read rather than
representing an absolute error probability. Compact v5 does not retain a redundant
`phred_calibrated` flag.

## Tests

- `zero_penalty_receives_maximum`: a zero penalty yields the maximum score.
- `worst_penalty_receives_zero`: the worst penalty yields zero.

## Status

Implemented.
