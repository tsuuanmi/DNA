# Interfaces and Boundaries

DNA separates six interface classes:

1. **Frontend/configuration boundary** — the CLI owns typed command arguments and
   resolves the current configuration path; parsed command values and the
   explicit configuration path are translated before entering operation
   orchestration.
2. **Input adapter boundary** — source-specific loaders such as
   `input::sanger` validate external sequencing/reference/configuration inputs
   and produce validated internal models without owning CLI output publication.
3. **Public Rust capability boundary** — typed reusable operations such as
   Sanger Variant Analysis return canonical domain results without CLI
   publication side effects.
4. **Scientific module boundaries** — typed internal models passed between
   decoding, calling, signal, QC, alignment, variant, and sample stages.
5. **Public serialized result boundary** — closed versioned JSON contracts.
6. **Filesystem/operational boundary** — atomic publication and append-only logs.

The CLI is an outer adapter: pipeline entry points receive operation values such
as paths and sample identifiers, not `clap` argument structs. Pipeline
operations also receive the resolved configuration path explicitly rather than
reading `DNA_CONFIG` or selecting a default themselves. This prevents the
current command-line/process environment from becoming an inward dependency of
orchestration or scientific modules and leaves the operation boundary usable by
future non-CLI callers.

The current Sanger filesystem adapter owns validation and loading of Sanger
sequencing traces, FASTA references, and explicit configuration into validated
internal models. Its current format layer is `input::sanger::abif`; support is
determined by the ABIF container and required sequencing tags rather than a
filename suffix. `PLOC.2` is decoded at that format boundary and projected to
canonical Sanger `locus_positions`; downstream scientific modules do not depend
on the ABIF tag name. It does not derive `results/*` paths, validate overwrite targets, select log paths,
or publish outputs. Those remain outer operation/publication concerns.

Operational stage logging crosses this boundary through a minimal internal
`StageLog` capability. Shared reference-free read processing and
reference-guided Variant Analysis can emit informational and warning records
without depending on the file-backed logger; log destination selection,
terminal error logging, and synchronization stay in the outer operation layer.

The production `pipeline` composes scientific capabilities but does not own
their implementations. Reference-free read processing is crate-internal shared
science; reference-guided one-read observation is owned by
`variant_analysis` and is reused by CLI analysis, sample evidence, and the
public Rust API.

The first public Rust capability is `variant_analysis::analyze_sanger`. It reuses the same scientific stages as the CLI but returns canonical typed variant evidence directly and does not open logs or publish JSON. Its exact contract is owned by [Rust public API](../reference/rust-api.md).

Exact syntax and serialized shapes are owned by
[reference](../reference/README.md). Module responsibilities are colocated under
[src](../../src/README.md). Cross-cutting rules that must hold across interfaces
are indexed under [invariants](invariants/README.md).

Architecture defines where responsibilities belong; it does not duplicate
algorithm details from [design](../design/README.md).
