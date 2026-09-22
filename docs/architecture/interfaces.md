# Interfaces and Boundaries

DNA separates four interface classes:

1. **CLI/configuration boundary** — typed command arguments and strict configuration.
2. **Scientific module boundaries** — typed internal models passed between decoding, calling, signal, QC, alignment, variant, and sample stages.
3. **Public result boundary** — closed versioned JSON contracts.
4. **Filesystem/operational boundary** — atomic publication and append-only logs.

Exact syntax and serialized shapes are owned by [reference](../reference/README.md). Module responsibilities are colocated under [src](../../src/README.md). Cross-cutting rules that must hold across interfaces are indexed under [invariants](invariants/README.md).

Architecture defines where responsibilities belong; it does not duplicate algorithm details from [design](../design/README.md).
