# `src/trace/mod.rs`

## Purpose

Decodes untrusted ABIF/AB1 binary input into a validated chromatogram.

## Responsibilities

- Re-export `load` as the module boundary.
- Verify the ABIF signature; bounds-check directory entries and payloads; extract
  A/C/G/T channels, positions, and available vendor calls and qualities.

## Non-responsibilities

No base re-calling, trimming, alignment, or reporting.

## Key types and functions

- `load(path) -> Result<Chromatogram>`: the public entry point, re-exported from
  `decode`.
- Child modules: `reader` (checked big-endian cursor), `abif` (directory parsing
  and tag lookup), `decode` (tag decoding into the model).

## Invariants and errors

All offsets and lengths are checked before slicing or allocation. Malformed,
truncated, empty, and unsupported files return typed errors, never panics. Vendor
PBAS/PCON data is retained as evidence only; P2BA.1 is not consumed.

## Dependencies

- `config` for `MAX_AB1_BYTES`.
- `checksum` for `hex_sha256`.
- `model::trace` for `Chromatogram` and `VendorEvidence`.
- `error` for `Error`/`Result`.

## Apollo mapping

`apollo/include/apollo/preprocessing/abif.h` parsing helpers and `readab`.

## Requirements and decisions

ADR-0003; `SRS-IN-001` through `SRS-IN-003`, `SRS-IN-006`,
`SRS-DATA-001` through `SRS-DATA-003`, and `SRS-VAL-001`.

## Tests

The integration tests exercise decoding through the pipeline.

## Status

Implemented.
