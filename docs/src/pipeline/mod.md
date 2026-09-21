# `src/pipeline/mod.rs`

## Purpose

Defines orchestration boundaries for production commands.

## Responsibilities

- Expose `analyze`, `basecall`, and `sample` as production command boundaries.
- Share validated input helpers, reference-independent read stages, reference-guided observation processing, multi-trace sample-read processing, and terminal operation/logging failure preservation.

## Non-responsibilities

No binary parsing, scoring loops, variant normalization, threshold fitting, or serialization-format logic.

## Key types and functions

- `analyze(args)`, `basecall(args)`, and `sample(args)`: production command entry points.
- Child modules separate input loading, read/observation science, sample-read reuse, sample metrics, and command orchestration.
- `record_failure`: shared terminal error-log and synchronization policy.

## Invariants and errors

- Production commands return success only after their result output is committed.
- Stage errors propagate as typed `Error` values.

## Dependencies

- `cli` for production command arguments.
- `error`, `logger`, and `Instant` for the shared failure boundary.

## Requirements and decisions

ADR-0001, ADR-0002, ADR-0007.

## Tests

Production integration tests cover analyze/basecall/sample.

## Status

Implemented.
