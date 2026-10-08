# ADR-0061: Use `tracing` for operational logging

- **Status:** Accepted
- **Date:** 2026-10-08
- **Related decisions:** [ADR-0007](0007-configuration-and-environment.md),
  [ADR-0059](0059-reuse-ecosystem-machinery-behind-dna-contracts.md)

## Context

SRS-OUT-008 requires append-only, run-correlated, one-line INFO/WARN/ERROR
stage records in one file per operation (`$DNA_LOG_DIR/<trace-stem>.log` or
`<sample-id>.log`), with nothing written to stdout/stderr, and SRS-API-004
forbids log files from the public Rust API.

DNA met this with a hand-written file logger: a `StageLog` trait threaded as a
`&mut` parameter through read processing and Variant Analysis, a no-op
implementation for the public API, explicit `module_path!()`/`line!()` call
sites, and a `&mut &'static str` stage label threaded through the same
functions so the terminal failure record could name the failing stage.
ADR-0059 asks DNA to reuse maintained ecosystem machinery for commodity
infrastructure, and structured instrumentation is commodity infrastructure.

## Decision

Operational logging uses the `tracing` ecosystem.

- Scientific stages emit structured `tracing` events (`event = "…"` plus
  `name = value` fields) and enter one `info` span per stage, named after the
  stage. They take no logger or stage parameter and never choose a destination.
- The CLI operation layer (`pipeline`) opens one `OperationLog` per operation
  and runs the operation inside a thread-scoped dispatcher made of a
  `tracing_subscriber` `Registry` and a DNA `RecordLayer`. The layer renders
  this crate's events at `INFO` or above in the existing record format, and
  remembers the last-entered stage span, which the terminal failure record
  reports.
- The public Rust API installs no subscriber. Library callers observe stage
  spans and events through their own subscriber, or not at all.
- Stages follow one another: each stage span is closed before the next opens,
  so subscribers see sequential, non-nested stages.
- A layer cannot return errors, so the first failed record write is captured
  together with the stage it occurred in, and surfaced as `Error::Log` when the
  operation checks or synchronizes the log: immediately after the start
  record, which keeps an unwritable log a fail-fast error, and before result
  publication. The terminal failure record names the stage of the first failed
  write.

## Consequences

- One mechanism replaces the custom trait, its no-op implementation, and the
  threaded stage label; module and line come from `tracing` metadata.
- The record format, file layout, and failure-stage attribution are unchanged.
- A write failure after the start record is reported at the next durability
  point instead of at the failing record, so the remaining stages still run.
  No result is published in either case.
- `tracing` caches callsite interest process-wide. A test that installs a
  scoped subscriber while other threads run the same instrumented code without
  one must run in its own test binary.
- Adds `tracing` (without its attribute macros) and `tracing-subscriber`
  (registry only, no formatting or environment-filter features) as
  dependencies.
