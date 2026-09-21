# `src/alignment/scoring.rs`

## Purpose

Owns deterministic fixed-point evidence-profile substitution scoring and packed
Gotoh state ordering.

## Responsibilities

- Define the `Match` / `Insertion` / `Deletion` state enum and traceback bit decoding.
- Scale every configured alignment delta into one fixed-point domain with denominator 1024.
- Score a canonical reference base from the normalized `EvidenceProfile` mass assigned to that base.
- Use the configured ambiguous score when the evidence profile is absent or the reference symbol is non-canonical.
- Preserve the negative-infinity sentinel during score addition.

## Scoring contract

For canonical reference base `r`:

~~~text
u = round(profile[r] * 1024)
substitution =
    u          * match_score
  + (1024-u)   * mismatch_score
~~~

Evidence weights are non-negative, so exact half-unit ties round upward under
Rust's round semantics. The dynamic program compares only integer `i64` scores.

All gap and ambiguous deltas are multiplied by the same 1024 scale. Therefore a
one-hot profile preserves the previous discrete match/mismatch ordering exactly
apart from the common scale factor.

## Non-responsibilities

No DP matrix traversal, traceback, orientation selection, QC admission, variant
extraction, or configuration loading.

## Key items

- `SCORE_SCALE = 1024`.
- `substitution(profile, reference, config) -> i64`.
- `scaled(delta) -> i64`.
- `is_canonical(base) -> bool`.
- `add(score, delta) -> i64`.

## Tests

Tests cover one-hot compatibility, mixed expected scores, missing/ambiguous
fallback semantics, and exact half-unit quantization.

## Traceability

ADR-0029; `SRS-ALN-002`, `SRS-ALN-003`, and `SRS-ALN-010`.

## Status

Implemented.
