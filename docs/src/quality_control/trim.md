# `src/quality_control/trim.rs`

## Purpose

Scores calls and selects one auditable retained interval by trimming low-quality
ends.

## Responsibilities

- Compute penalties and relative quality scores for every call.
- Determine left and right trim bounds by walking outward from the best section
  until a window exceeds the stringency threshold.
- Build the per-call quality records and the retained sequence.

## Non-responsibilities

No middle-read masking, denoising, ML scoring, or variant filtering.

## Key types and functions

- `analyze(trace, calls, config) -> Result<QualityControlResult>`: the module
  entry point, re-exported from `mod.rs`.

## Invariants and errors

- The call count must be at least `minimum_retained_bases`; otherwise
  `Error::QualityControl`.
- The retained interval must be at least `minimum_retained_bases` long;
  otherwise `Error::QualityControl`.
- One quality record per call; trim bounds are within the call range.
- `vendor_quality_applies` is true only when PCON is present and PBAS agrees with
  the signal-derived call. The vendor score itself is not copied into
  `CallQuality` or compact v5 output.

## Dependencies

- `penalty` and `quality`.
- `config` for `QualityControlConfig`.
- `model::basecalls`, `model::quality`, `model::trace`.
- `error` for `Error`/`Result`.

## Biological semantics

Read ends are typically lower quality. Trimming removes low-quality tails while
retaining the best-supported section, so the aligned read reflects high-confidence
bases.

## Tests

- `retains_a_uniform_read_and_applies_matching_vendor_quality`: verifies full-span
  retention and the applicability boolean for present PCON whose PBAS agrees with
  the signal-derived primary call; it does not retain the PCON value.

## Status

Implemented.
