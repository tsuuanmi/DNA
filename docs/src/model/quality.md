# `src/model/quality.rs`

## Purpose

Defines the relative quality scores and auditable end-trim bounds produced by
quality control.

## Responsibilities

- Represent per-call quality evidence and the retained interval.

## Non-responsibilities

No scoring, trimming, or algorithm logic.

## Key types and functions

- `CallQuality`: call index, penalty, uncalibrated relative quality score, and the
  boolean indicating whether optional vendor quality applies. The vendor score
  itself and a redundant calibration flag are not retained.
- `QualityControlResult`: per-call records, trim bounds, and the retained
  sequence.

## Invariants and errors

- One quality record per call.
- Trim bounds are within the call range and retain at least the configured
  minimum number of bases.

## Dependencies

No external dependencies; report serialization uses separate compact result types.

## Biological semantics

Quality scores are relative and uncalibrated; they rank calls by local ambiguity
and peak-spacing penalties. No `phred_calibrated` field is retained because the
method is always uncalibrated. The retained interval is the auditable read section
kept for alignment.

## Tests

No dedicated unit tests; behavior is exercised through `quality_control`.

## Status

Implemented.
