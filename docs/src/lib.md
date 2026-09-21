# `src/lib.rs`

## Purpose

Defines the DNA library boundary and dispatches parsed commands.

## Responsibilities

- Declare the internal subsystem graph, including the private operational logger.
- Expose the intentional CLI, configuration, error, model, and explicit local-validation surfaces.
- Route commands to `pipeline`.

## Inputs and outputs

Accepts `cli::Cli`; returns `error::Result<()>`.

## Key types and functions

- `run(cli: Cli) -> Result<()>`: dispatches `Command::Analyze`, `Command::Basecall`, and `Command::Sample` to their pipeline entry points.

## Invariants and errors

The library forbids unsafe Rust and denies deprecated API use. The separate Rust source-policy gate rejects first-party `allow`/`expect(deprecated)` suppressions. It does not print or select process exit codes.
Pipeline errors are preserved for the binary boundary.

## Dependencies

- `cli`, `config`, `error`, and `model` are public.
- `alignment`, `basecalling`, `locus`, `logger`, `pipeline`, `quality_control`, `reference`, `report`, `sample`, `signal_processing`, `trace`, and `variant_calling` are private.

## Apollo mapping

Replaces the routing role of `apollo/src/apollo.cpp`; it does not reproduce
Apollo's command set.

## Requirements and decisions

ADR-0001, ADR-0002, ADR-0006; `SRS-IN-010`, `SRS-NFR-001`, `SRS-NFR-005`, `SRS-NFR-006`.

## Tests

`tests/cli.rs`, `tests/docs_mirror.rs`.

## Status

Implemented.
