# `src/config/mod.rs`

## Purpose

Owns strict, reproducible scientific configuration: compiled safety caps, typed
validated records, and strict TOML loading with source identity.

## Responsibilities

- Re-export the intentional configuration API: `load`, the compiled caps, and the
  typed per-stage config structs. SHA-256 identity uses the shared `checksum`
  module.
- Centralize the 1 MiB config, 4 MiB FASTA source, 64 MiB AB1, 50 bp indel,
  50,000-base reference, and alignment-cell caps.

## Non-responsibilities

No `.env` parsing, per-setting environment overrides, CLI path parsing, silent
fallback from an invalid file, or hidden working-directory search beyond the
single documented default path.

## Key types and functions

- `load() -> Result<Config>`: resolves `DNA_CONFIG` (or the default path),
  parses strict TOML, and validates.
- `Config`, `ReferenceConfig`, `BasecallingConfig`, `DNAProcessingConfig`,
  `QualityControlConfig`, `AlignmentConfig`, `SampleReconciliationConfig`, `VariantCallingConfig`: typed validated records.
- Re-exported compiled caps: `MAX_AB1_BYTES`, `MAX_REFERENCE_BYTES`,
  `MAX_REFERENCE_LENGTH`, and `MAX_ALIGNMENT_CELLS`; defaults also owns internal
  config-source and indel caps.

## Invariants and errors

Configuration is strict: unknown keys and unsupported schema versions are
rejected, and out-of-range values return `Error::Config` rather than being
silently clamped. The effective config records its source path and SHA-256 for
reproducibility.

## Dependencies

- `defaults`, `load`, `types`.
- `error` for `Error`/`Result`.

## Apollo mapping

Consolidates defaults that are spread across Apollo command configuration
structs.

## Requirements and decisions

ADR-0002, ADR-0004, ADR-0007; `SRS-CFG-001` through `SRS-CFG-007`,
`SRS-ALN-003`, and `SRS-VAR-002`.

## Tests

Unit tests in `load` cover unknown-key rejection. The end-to-end integration
tests exercise configuration loading and validation.

## Status

Implemented.
