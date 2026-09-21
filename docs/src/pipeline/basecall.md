# `src/pipeline/basecall.rs`

## Purpose

Orchestrates one reference-free AB1-to-basecalls JSON operation.

## Responsibilities

- Open the per-trace append-only logger and record command-specific start,
  input, warning, publication, and failure events.
- Load basecall inputs, delegate shared scientific stages to `pipeline::read`,
  include trace-integrity warning counts in the command summary, build and
  serialize the typed result, synchronize mandatory logs, and publish atomically
  without overwrite.

## Non-responsibilities

No reference loading, alignment, variant calling, scientific algorithms, output
format selection, or sequence payload logging.

## Dependencies

`cli`, `error`, `logger`, `pipeline::{input,read}`, and `report`.

## Tests

`tests/basecall.rs` covers deterministic success, malformed input, no-overwrite,
coexistence with analysis, and log privacy.
