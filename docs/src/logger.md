# `src/logger.rs`

## Purpose

Provides append-only operational logging for one command operation.

## Responsibilities

- Resolve `DNA_LOG_DIR`, defaulting to `logs/`.
- Open `<operation-stem>.log` in append mode. Single-read commands use the trace stem; sample commands use `<sample-id>` so all nested trace-stage records share one sample log.
- Assign a per-open run identifier so records in an append history can be grouped.
- Format INFO, WARN, and ERROR records with a local millisecond timestamp, Rust
  module/line source location, and automatic `run_id` field.
- Escape record delimiters and control characters so each event occupies exactly
  one physical line.
- Flush and synchronize mandatory records before result publication.

## Non-responsibilities

No scientific decisions, sequence/peak serialization, JSON output, CLI parsing,
or recovery policy.

## Key types and functions

- `Logger::open(operation_stem) -> Result<Logger>`: creates the directory and opens the operation-scoped file.
- `info`, `warn`, `error`: append one Apollo-style record.
- `sync`: flushes and synchronizes the file.

## Invariants and errors

Records use `YYYY-MM-DD HH:MM:SS.mmm | LEVEL | module:line - run_id=<id> message`.
The logger writes no records to stdout or stderr. Files are append-only, and
repeated invocations for the same stem and directory share one history. An empty
`DNA_LOG_DIR` returns `Error::Path`; directory/file I/O returns `Error::Log`.
Open, write, and pre-publication synchronization failures are fatal. Concurrent
processes may append to the same file, but cross-process record order is not
guaranteed.

## Privacy

Messages may contain trace/reference/config/output paths, trace and reference
names, hashes, aggregate metrics, thresholds, timings, stage errors, and removed
variant kinds/coordinates/reasons. They omit complete sequences, alleles,
configured region contents, individual peak/call records, alignment strings, and
JSON bodies. `logs/` is ignored and follows the local AB1 privacy policy.

## Tests

Unit tests use temporary directories and cover append behavior plus one-line
escaping. CLI integration tests cover ordered success events, warning/removal
records, custom destinations, and stage-aware failures.

## Status

Implemented.
