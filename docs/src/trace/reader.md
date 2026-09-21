# `src/trace/reader.rs`

## Purpose

Provides bounds-checked big-endian reads over untrusted ABIF bytes.

## Responsibilities

- Wrap immutable input bytes and expose checked byte slices and big-endian
  integer reads.

## Non-responsibilities

No ABIF structure parsing or tag semantics.

## Key types and functions

- `Reader<'a>`: a view over `&'a [u8]`.
- `slice(offset, length) -> Result<&'a [u8]>`: checked byte slice.
- `u16(offset)`, `i16(offset)`, `u32(offset)`: checked big-endian reads.

## Invariants and errors

- Every read is bounds-checked; out-of-range or overflow returns `Error::Abif`.
- No read can panic on malformed input.

## Dependencies

- `error` for `Error`/`Result`.

## Biological semantics

None; this is a low-level binary cursor.

## Tests

- `reads_big_endian_signed_and_unsigned_values` covers integer decoding.
- `rejects_out_of_bounds_and_overflowing_ranges` covers both range failure modes.

## Status

Implemented.
