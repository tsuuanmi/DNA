# `src/error/mod.rs`

## Purpose

Defines typed failures crossing the DNA library boundary.

## Responsibilities

- Provide the shared `Result<T>` alias and errors that callers can handle without
  parsing log text.

## Key types and functions

- `Result<T> = std::result::Result<T, Error>`.
- `Error` enum with variants for path, read, config parse, config value, ABIF,
  FASTA, basecalling, signal processing, quality control, alignment, variant,
  phase measurement, report assembly,
  serialization, logging I/O, combined operation/logging failure, and output
  failures.

## Invariants and errors

Errors describe failures; they do not print, exit, or discard source context.
Errors do not contain patient trace data or full signal arrays; they may include a
safe path and stage name.

## Dependencies

- `thiserror` for error derivation.
- `toml` and `serde_json` for source error types.

## Apollo mapping

Replaces raw and inconsistently used Apollo integer error codes.

## Requirements and decisions

ADR-0002; `SRS-IN-006`, `SRS-IN-012`, `SRS-NFR-002`.

## Tests

`tests/cli.rs` verifies boundary behavior.

## Status

Implemented.
