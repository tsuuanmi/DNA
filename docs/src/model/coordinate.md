# `src/model/coordinate.rs`

## Purpose

Provides explicit, checked coordinate conversion at external reporting
boundaries.

## Responsibilities

- Convert a zero-based reference index to a checked one-based position.

## Non-responsibilities

No coordinate arithmetic beyond the single conversion; no I/O or algorithm logic.

## Key types and functions

- `reference_one_based(position_0based) -> Result<usize>`: adds one with overflow
  checking.

## Invariants and errors

- The conversion uses `checked_add`; overflow returns `Error::Variant`.

## Dependencies

- `error` for `Error`/`Result`.

## Biological semantics

Internal positions are zero-based; reported variant positions are one-based, as
required by the output contract. This module is the single place where that
conversion happens.

## Tests

No dedicated unit tests; behavior is exercised through `variant_calling`.

## Status

Implemented.
