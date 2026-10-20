# Architecture

Architecture describes **current system structure**: stable boundaries, context, data flow, interfaces, dependency direction, and cross-cutting invariants. Diagrams are Mermaid blocks in the document that owns the fact: the [context diagram](context.md#context-diagram), the [crate dependency graph](overview.md#module-boundaries), and the [data-flow](data-flow.md) diagrams (one read, a sample, the other entries, post-calling, step contracts).

- [Overview](overview.md)
- [System context](context.md)
- [Data flow](data-flow.md)
- [Interfaces and boundaries](interfaces.md)
- [Variant lifecycle](variant-lifecycle.md)
- [System invariants](invariants/README.md)
- [Decision history](../decisions/README.md)
- [Source navigation](../../src/README.md) and [workspace crates](../../crates/README.md)

Detailed mechanisms belong in [design](../design/README.md). Exact public/configuration shapes belong in [reference](../reference/README.md). Architecture should not duplicate either.
