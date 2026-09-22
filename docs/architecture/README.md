# Architecture

- [System overview](system.md): current module boundaries, data flow, dependency direction, resource bounds, and publication model.
- [System invariants](invariants/README.md): canonical index of cross-cutting invariant families.
- [Source mirror policy](source-mirror.md): one-to-one `src/**/*.rs` → `docs/src/**/*.md` implementation documentation.
- [Architecture decisions](../adr/README.md): accepted, superseded, and historical durable decisions.

The source tree and `docs/src/` mirror are the canonical implementation layout. Architecture documentation should describe stable boundaries and dependency direction rather than maintain a second hand-written file inventory.

Detailed scientific algorithms live in [methods](../methods/README.md); implementation details live in `docs/src/`.
