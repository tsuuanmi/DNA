# `src/quality_control/mod.rs`

## Purpose

Estimates relative base quality and selects a retained read interval.

## Responsibilities

- Re-export `analyze` as the module boundary.
- Compute bounded per-call quality, identify low-quality left/right tails, and
  return auditable trim boundaries.

## Non-responsibilities

No middle-read masking, denoising, ML scoring, or variant filtering.

## Key types and functions

- `analyze(trace, calls, config) -> Result<QualityControlResult>`: the public
  entry point, re-exported from `trim`.
- Child modules: `penalty` (per-position penalties and best section), `quality`
  (bounded score conversion), `trim` (end trimming and retained interval).

## Invariants and errors

One quality per call; trim bounds are within the call range and retain at least
`minimum_retained_bases`. Invalid or entirely unusable reads return
`Error::QualityControl`.

## Dependencies

- `config` for `QualityControlConfig`.
- `model::basecalls`, `model::quality`, `model::trace`.
- `error` for `Error`/`Result`.

## Apollo mapping

Quality helpers in `preprocessing/abif.h` and end-trimming logic in
`quality_control/trim.h`.

## Requirements and decisions

ADR-0001, ADR-0003; `SRS-QC-001` through `SRS-QC-004`.

## Tests

Unit tests in `quality` cover score boundaries. The integration tests exercise
the full quality-control path.

## Status

Implemented.
