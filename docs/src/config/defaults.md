# `src/config/defaults.rs`

## Purpose

Defines the non-scientific default path and the compiled hard resource caps.

## Responsibilities

- Provide the default strict TOML path used when `DNA_CONFIG` is unset.
- Provide the compiled safety caps: maximum config, AB1, and FASTA source bytes;
  maximum reference/indel length, representable peak height, and maximum alignment
  cells.

## Non-responsibilities

No parsing, validation, or environment handling.

## Key constants

- `DEFAULT_CONFIG_PATH`: `config/dna.toml`.
- `MAX_CONFIG_BYTES`: `1024 * 1024` (1 MiB).
- `MAX_AB1_BYTES`: `64 * 1024 * 1024` (64 MiB).
- `MAX_REFERENCE_BYTES`: `4 * 1024 * 1024` (4 MiB).
- `MAX_REFERENCE_LENGTH`: `50_000` bases.
- `MAX_INDEL_LENGTH`: `50` bases.
- `MAX_PEAK_HEIGHT`: `32767`, the largest ABIF signed-short channel value.
- `MAX_ALIGNMENT_CELLS`: `100_000_000`.

## Invariants and errors

These are compile-time constants. They are enforced by the modules that consume
them (`trace::decode`, `reference::fasta`, `config::types`, `alignment::gotoh`).
`MAX_REFERENCE_LENGTH` is validated to be non-zero during configuration
validation.

## Dependencies

None.

## Biological semantics

The caps bound resource use for direct alignment of a single read against a
reference. `MAX_REFERENCE_LENGTH` reflects the direct-alignment design for
references such as the 16.6 kb mitochondrial genome; `MAX_INDEL_LENGTH` bounds
the largest reported primary-sequence indel.

## Tests

No dedicated unit tests; the constants are exercised by the modules that consume
them.

## Status

Implemented.
