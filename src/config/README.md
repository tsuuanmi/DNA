# Configuration

Owns strict configuration loading, validation, typed per-stage settings, and
compiled resource caps.

Key children: `defaults.rs`, `load.rs`, and `types.rs`.

This module does not provide per-setting environment fallbacks or silently clamp
invalid values.

See [configuration requirements](../../docs/srs/configuration.md) and
[configuration contract](../../docs/contracts/configuration.md).
