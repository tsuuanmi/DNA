# `src/main.rs`

## Purpose

Provides the operating-system process boundary for the `signal` binary.

## Responsibilities

- Parse CLI arguments, call the library, print one concise stderr error, and map
  success or failure to a process exit status.

## Non-responsibilities

No scientific logic, file parsing, output writing, or configuration defaults.

## Key types and functions

- `main() -> ExitCode`: parses `Cli`, calls `signal::run`, and returns
  `ExitCode::SUCCESS` or `ExitCode::FAILURE`.

## Invariants and errors

Help/version are owned by Clap. The binary crate forbids unsafe Rust and denies deprecated API use; first-party suppression is rejected by the Rust source-policy gate. Library failures produce a nonzero exit and are never converted into fake success.

## Dependencies

- `clap` for `Parser`.
- `signal::cli::Cli`.

## Apollo mapping

Corresponds only to the process-boundary portion of `apollo/src/apollo.cpp`.

## Requirements and decisions

ADR-0002; `SRS-IN-011`, `SRS-IN-012`, `SRS-NFR-001`, `SRS-NFR-005`.

## Tests

`tests/cli.rs`.

## Status

Implemented.
