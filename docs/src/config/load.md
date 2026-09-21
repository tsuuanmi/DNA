# `src/config/load.rs`

## Purpose

Resolves, parses, and validates the one authoritative strict TOML configuration
file, and computes its identity.

## Responsibilities

- Resolve the configuration path from `DNA_CONFIG`, falling back to
  `DEFAULT_CONFIG_PATH`.
- Reject a file above the compiled 1 MiB source cap before reading it.
- Read the file, require UTF-8, and parse it as strict TOML.
- Validate the parsed `RawConfig` into a typed `Config`, attaching the source path
  and SHA-256 identity computed with the shared `checksum` module.

## Non-responsibilities

No per-setting environment overrides, `.env` parsing, or silent fallback from an
invalid file.

## Key types and functions

- `load() -> Result<Config>`: the public entry point, re-exported from `mod.rs`.
- `load_path(path) -> Result<Config>`: loads a specific path.

## Invariants and errors

- A missing or unreadable file returns `Error::Read`.
- Non-UTF-8 bytes return `Error::Config`.
- Invalid TOML returns `Error::ConfigParse`.
- Validation failures (unknown keys, unsupported schema version, missing signal/reconciliation
  settings, out-of-range values) return `Error::Config`.
- The effective config records its source path and SHA-256 for reproducibility.

## Dependencies

- `defaults` for `DEFAULT_CONFIG_PATH`.
- `types` for `Config` and `RawConfig`.
- `checksum` for `hex_sha256`.
- `error` for `Error`/`Result`.

## Biological semantics

None; this module is purely about configuration identity and validation.

## Tests

- `rejects_unknown_keys`: verifies that an unknown TOML key causes validation to
  fail.

## Status

Implemented.
