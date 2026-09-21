# ADR-0007: Use Strict TOML Configuration and a Safe Environment Policy

- **Status:** Accepted
- **Date:** 2026-08-22

## Context

Scientific defaults must be reviewable and reproducible. Apollo spreads defaults
across command structs, JSON, and environment files with conflicting values.
Committing machine-local `.env` files also creates a future secret and portability
risk.

## Options

1. Keep all settings compiled into Rust.
2. Copy Apollo's JSON and tracked `.env` behavior.
3. Use a versioned TOML contract, an environment-selected path, and a tracked
   `.env.example` while ignoring local `.env`.

## Decision

Choose option 3.

- The tracked default is `config/dna.toml` with `schema_version = 4`.
- `DNA_CONFIG` may select another TOML path.
- Environment variables do not override individual scientific values.
- The operational `DNA_LOG_DIR` may select the append-only log destination;
  it is not scientific configuration and does not affect the config checksum.
- The CLI continues to require the AB1 path and `--reference` path.
- Unknown keys, invalid values, and unsupported versions are errors.
- `.env.example` is tracked; `.env` is local and ignored.
- The binary does not automatically parse `.env`.

## Consequences

Configuration is human-readable, typed, and reproducible. Runs identify the
exact config by checksum and versioned method IDs; effective values remain in the
selected strict TOML. Local tooling remains convenient without committing machine
state. The implemented loader uses strict TOML deserialization and explicit range
validation.

## Supersession

ADR-0013 updates the strict schema version to 4 for required observational
signal settings and the minimum noisy-window run length. ADR-0011 partially supersedes the earlier "record every effective value" goal:
the compact JSON records the config checksum and versioned method IDs rather than
duplicating every effective value, which stays in the TOML. ADR-0007's strict
TOML selection, environment policy, and deterministic resolution remain in
force. A new configuration source or precedence layer still requires an ADR.
