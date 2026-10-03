# Interfaces and Boundaries

DNA separates four interface classes:

1. **Frontend/configuration boundary** — the CLI owns typed command arguments and
   resolves the current configuration path; parsed command values and the
   explicit configuration path are translated before entering operation
   orchestration.
2. **Scientific module boundaries** — typed internal models passed between
   decoding, calling, signal, QC, alignment, variant, and sample stages.
3. **Public result boundary** — closed versioned JSON contracts.
4. **Filesystem/operational boundary** — atomic publication and append-only logs.

The CLI is an outer adapter: pipeline entry points receive operation values such
as paths and sample identifiers, not `clap` argument structs. Pipeline
operations also receive the resolved configuration path explicitly rather than
reading `DNA_CONFIG` or selecting a default themselves. This prevents the
current command-line/process environment from becoming an inward dependency of
orchestration or scientific modules and leaves the operation boundary usable by
future non-CLI callers.

Single-read analysis input loading also does not own the JSON publication
destination. Output-path derivation and overwrite protection remain an outer
operation/publication concern, while decoded trace, reference, and configuration
form the scientific input boundary.

Exact syntax and serialized shapes are owned by
[reference](../reference/README.md). Module responsibilities are colocated under
[src](../../src/README.md). Cross-cutting rules that must hold across interfaces
are indexed under [invariants](invariants/README.md).

Architecture defines where responsibilities belong; it does not duplicate
algorithm details from [design](../design/README.md).
