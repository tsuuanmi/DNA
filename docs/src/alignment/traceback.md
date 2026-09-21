# `src/alignment/traceback.rs`

## Purpose

Decodes the packed Gotoh traceback into gapped alignment rows and computes
alignment metrics.

## Responsibilities

- Walk the packed traceback from a chosen endpoint back to the query start,
  producing ordered `RawColumn` records.
- Build internal gapped query and reference strings used to distinguish equal-score
  placements; operation-run serialization is not retained.
- Count exact matches, mismatches, gap opens, callable columns, and unresolved
  query bases, and derive callable identity.

## Non-responsibilities

No dynamic-programming scoring, strand selection, or variant extraction.

## Key types and functions

- `RawColumn`: one aligned column with query/reference bases and optional
  original indices.
- `RawAlignment`: score, reference span, aligned columns, and `AlignmentMetrics`;
  gapped-row and operation-run strings are not duplicated.
- `TracebackInput`: the packed trace, dimensions, endpoint, state, and score.
- `decode(input) -> Result<RawAlignment>`: performs the traceback walk.
- `metrics(columns) -> AlignmentMetrics`: derives the single authoritative exact/mismatch/gap/callable metric set and is reused after canonical gap shifts.
- `gap_open_count(columns) -> usize`: counts transitions into insertion/deletion runs without serializing an operation string.

## Invariants and errors

- Every traceback index is bounds-checked; overflow or out-of-range access
  returns `Error::Alignment`.
- Reaching a match or deletion state at reference column zero is an error.
- The ordered columns are the single alignment-row representation; each is a
  match, insertion (query base, reference `-`), or deletion (query `-`, reference
  base).
- Callable identity is `exact_matches / callable_columns`, or `0.0` when there
  are no callable columns.

## Dependencies

- `scoring` for `State` and `is_canonical`.
- `model::alignment` for `AlignmentMetrics`.
- `error` for `Error`/`Result`.

## Biological semantics

Callable columns are those where both query and reference bases are canonical;
unresolved query bases (`N`) are counted separately and excluded from identity.
This keeps ambiguous calls from inflating or deflating the reported identity.

## Tests

No dedicated unit tests; behavior is exercised through `gotoh` and `orient`.

## Status

Implemented.
