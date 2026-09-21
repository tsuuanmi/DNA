# `src/config/types.rs`

## Purpose

Defines the typed, validated configuration records consumed by the scientific
stages.

## Responsibilities

- Define the effective `Config` and the per-stage configuration structs
  (`ReferenceConfig`, `BasecallingConfig`, `DNAProcessingConfig`,
  `QualityControlConfig`, `AlignmentConfig`, `SampleReconciliationConfig`, `VariantCallingConfig`).
- Define the strict `RawConfig` deserialization shape with unknown-key rejection.
- Validate every value against the scientific contract and produce a typed
  `Config` with source identity.

## Non-responsibilities

No file I/O, environment handling, or algorithm execution.

## Key types and functions

- `Config`: complete effective configuration plus `source_path` and
  `source_sha256`.
- `RawConfig::validate(source_path, source_sha256) -> Result<Config>`: the
  validation entry point.
- `require_fraction` and `require_finite_range`: shared range validators.

## Invariants and errors

- `schema_version` must equal `5`; otherwise `Error::Config`.
- Fractions (`secondary_peak_ratio`, `best_section_fraction`, `minimum_identity`)
  must be finite and in `(0, 1]`.
- `window_size_bases` must be in `5..=10`, `minimum_primary_snr` must be
  finite and positive, and `minimum_noisy_windows` must be at least 2.
- `trim_window_size`, `max_relative_quality_score`, `minimum_retained_bases`,
  `minimum_callable_bases` and `sample_reconciliation.minimum_comparable_bases` must be positive.
- `trim_stringency` must be finite and in `[0, 9]`.
- `match_score` must be positive; `mismatch_score`, `gap_open_score`, and
  `gap_extension_score` must be negative.
- `sample_reconciliation.minimum_overlap_agreement` must be finite and in `(0, 1]`.
- `max_indel_length` must be in `1..=MAX_INDEL_LENGTH` and
  `minimum_peak_height` in `1..=MAX_PEAK_HEIGHT`.
- `relative_quality_threshold` must be below `max_relative_quality_score`.
- `regions` must be a non-empty list of inclusive `[start, end]` pairs within
  `1..=MAX_REFERENCE_LENGTH`.
- All raw structs use `#[serde(deny_unknown_fields)]`, so unknown keys are
  rejected at parse time.

## Dependencies

- `defaults` for `MAX_INDEL_LENGTH`, `MAX_PEAK_HEIGHT`, and
  `MAX_REFERENCE_LENGTH`.
- `model::reference` for `ReferenceTopology`.
- `error` for `Error`/`Result`.
- `serde` for serialization/deserialization.

## Biological semantics

The alignment scores encode the biological cost model: matches are rewarded,
mismatches and gaps are penalized, and ambiguous bases are scored neutrally. The
secondary-peak ratio and trim settings control how mixed positions and low-quality
read ends are handled.

## Tests

Unit tests cover schema version 5, required signal/reconciliation/filter fields, signal bounds,
overlap-admission bounds, list-of-lists parsing, threshold relationships, and invalid/empty/out-of-range
regions. End-to-end tests exercise the selected strict config.

## Status

Implemented.
