# Architecture

- [System overview](system.md): current module boundaries, data flow, dependency direction, resource bounds, and publication model.
- [System invariants](invariants/README.md): canonical index of cross-cutting invariant families.
- [Architecture decisions](../adr/README.md): accepted, superseded, and historical durable decisions.
- [Source navigation](../../src/README.md): implementation ownership colocated with the Rust source tree.

Architecture describes stable system structure and dependency direction. Detailed
scientific algorithms live in [methods](../methods/README.md), while implementation
ownership lives in the nearest source-directory `README.md`.

Do not maintain a second hand-written source inventory under `docs/`.
