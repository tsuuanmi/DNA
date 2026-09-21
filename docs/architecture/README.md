# Architecture

- [System overview](system.md): current module boundaries, data flow, dependency direction, resource bounds, and publication model.
- [System invariants](invariants.md): cross-cutting truths that every module must preserve.
- [Source layout](source-layout.md): current Rust module ownership and dependency direction.
- [Source mirror policy](source-mirror.md): one-to-one `src/**/*.rs` → `docs/src/**/*.md` implementation documentation.
- [Architecture decisions](../adr/README.md): accepted, superseded, and historical durable decisions.

Architecture describes where responsibilities belong. Detailed scientific algorithms live in [methods](../methods/README.md); implementation details live in `docs/src/`.
