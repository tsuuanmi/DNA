# `src/trace/decode.rs`

## Purpose

Reads and decodes one canonical analyzed ABIF/AB1 file into a validated
chromatogram.

## Responsibilities

- Read the file, enforce the size bounds, and compute the source SHA-256.
- Parse the ABIF container and decode the required tags: `FWO_.1` (channel
  order), `DATA.9-12` (channels), `PLOC.2` (base locations), and optional
  `PBAS.2` and `PCON.2` (vendor evidence). `P2BA.1` is ignored.
- Reorder channels into canonical A/C/G/T order and validate channel/PLOC structure while preserving optional vendor-series cardinality independently.

## Non-responsibilities

No base re-calling, trimming, alignment, or reporting.

## Key types and functions

- `load(path) -> Result<Chromatogram>`: the module entry point, re-exported from
  `mod.rs`.
- `decode(...) -> Result<Chromatogram>`: tag decoding into the model.
- `require_layout`, `decode_i16`, `decode_optional_string`,
  `decode_optional_bytes`, `channel_index`: helpers.

## Invariants and errors

- The file must be non-empty and at most `MAX_AB1_BYTES`; otherwise
  `Error::Abif`.
- `FWO_.1` must be a four-base A/C/G/T permutation.
- `DATA.9-12` channels must be non-empty and equally sized.
- `PLOC.2` must be non-empty, strictly increasing, and within the sample range.
- Vendor strings must use uppercase DNA/RNA IUPAC symbols. Optional PBAS/PCON
  length mismatch is retained for downstream integrity evidence and does not
  create/remove PLOC-defined loci.
- PCON.2 may use the ABIF byte or char one-byte element type.
- All offsets and lengths are checked; malformed input returns `Error::Abif`,
  never a panic.
- Vendor PBAS/PCON data is retained as evidence only; P2BA.1 is not consumed.

## Dependencies

- `config` for `MAX_AB1_BYTES`.
- `checksum` for `hex_sha256`.
- `abif` and `reader`.
- `model::trace` for `Chromatogram` and `VendorEvidence`.
- `error` for `Error`/`Result`.

## Biological semantics

The decoded chromatogram provides the four fluorescence channels and base
locations that drive DNA's own re-calling. Vendor calls are retained for
comparison but never used as the final result.

## Tests

No dedicated unit tests; behavior is exercised through the integration tests.

## Status

Implemented.
