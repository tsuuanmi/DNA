# Configuration

Owns strict loading of the configuration envelope: the schema version, the
target profile path, the source identity, and the composition of every
plugin-owned section ([ADR-0069](../../docs/decisions/adr/0069-plugin-first-modality-core-post-calling.md)).
Each section's raw record, validated record, and rules live with the stage that
owns it (`alignment::config`, `basecalling::config`, …). The Sanger and core
plugins validate their sections together as `SangerConfig` and `CoreConfig`.
Resource caps live with the code they bound.

Key children: `defaults.rs` (default path and size cap), `load.rs`, and
`types.rs` (the envelope).

Configuration path selection is an outer application concern. The current CLI
adapter resolves `DNA_CONFIG` or the default path, while pipeline operations
receive that path explicitly and never inspect process environment themselves.

The path itself is not scientific identity; reproducibility continues to use the
validated configuration content checksum recorded in result provenance.

The configuration holds method parameters only. Its root `profile` key names the
target profile, resolved against the configuration file's directory; loading
and validating that profile belongs to `profile`.

This module does not provide per-setting environment fallbacks or silently clamp
invalid values.

See [configuration requirements](../../docs/requirements/configuration.md) and
[configuration contract](../../docs/reference/configuration.md).
