# `src/checksum.rs`

## Purpose

Provides the single stable SHA-256 identity helper shared by input and
configuration loading.

## Responsibilities

- Return a lowercase hexadecimal SHA-256 digest of validated input bytes.
- Be the one source of the `hex_sha256` identity used across the pipeline.

## Non-responsibilities

No file I/O, parsing, path resolution, or policy.

## Key types and functions

- `hex_sha256(bytes: &[u8]) -> String`: lowercase hex SHA-256 digest.

## Invariants and errors

- The function is total and never fails.
- The digest is deterministic: identical bytes always produce the same string.

## Dependencies

- `sha2` for `Sha256` and `Digest`.

## Biological semantics

None; this module is purely about reproducible byte identity.

## Tests

- `hashes_bytes_deterministically`: a known byte input yields the expected digest.

## Status

Implemented.
