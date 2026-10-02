# Configuration

Owns strict configuration loading, validation, typed per-stage settings, and
compiled resource caps.

Key children: `defaults.rs`, `load.rs`, and `types.rs`.

Configuration path selection is an outer application concern. The current CLI
adapter resolves `DNA_CONFIG` or the default path, while pipeline operations
receive that path explicitly and never inspect process environment themselves.

This module does not provide per-setting environment fallbacks or silently clamp
invalid values.

See [configuration requirements](../../docs/requirements/configuration.md) and
[configuration contract](../../docs/reference/configuration.md).
