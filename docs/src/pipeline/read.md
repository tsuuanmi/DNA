# `src/pipeline/read.rs`

## Purpose

Provides the one shared orchestration path for reference-independent scientific
stages.

## Responsibilities

- Run signal-derived base calling, observational rolling signal analysis plus basecall-independent locus and trace-integrity evidence extraction, and relative quality/end trimming in order.
- Emit aggregate stage metrics including integrity warning counts and event-signal ratio without sequences or peak arrays.
- Return `ProcessedRead` with calls, signal analysis, quality control, and
  command-level warning counts.

## Non-responsibilities

No input loading, reference use, alignment, variants, report construction,
publication, or terminal error policy.

## Dependencies

`basecalling`, `signal_processing`, `quality_control`, their model/config types,
`logger`, and `error`.

## Tests

The existing scientific unit tests plus `tests/analyze.rs` and
`tests/basecall.rs` exercise both callers of this shared path.
