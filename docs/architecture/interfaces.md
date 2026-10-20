---
decided-by: [ADR-0007, ADR-0058, ADR-0059, ADR-0061, ADR-0069, ADR-0070, ADR-0075]
---

# Interfaces and Boundaries

DNA separates biological scope from software boundaries. Biological target is independent of sequencing modality and format: mtDNA, nuclear/genomic DNA, and targeted loci may reuse capabilities without making the platform mtDNA-specific. SNP/genotyping is an analysis use case rather than an input modality.

DNA separates six interface classes:

1. **Frontend/configuration boundary** — the CLI owns typed command arguments and
   resolves the current configuration path; parsed command values and the
   explicit configuration path are translated before entering operation
   orchestration.
2. **Input adapter boundary** — source-specific loaders such as
   `input::sanger` validate external sequencing/reference/configuration inputs
   and produce validated internal models without owning CLI output publication. Recorded contract documents (`dna.read_evidence/v1`,
   `dna.consensus/v1`, `dna.variants/v1`, `dna.signal/v1`) are input adapters too:
   a step can be replayed from the contract the previous step wrote.
3. **Public Rust capability boundary** — typed reusable operations such as
   Sanger Variant Analysis return stable typed called-variant results without CLI
   publication side effects.
4. **Scientific module boundaries** — typed internal models passed between
   decoding, calling, signal, callability, QC, alignment, variant, and sample
   stages. The modality evidence boundary is one of them: a modality plugin
   hands the core caller per-read `ReadEvidence`, and the core sees no modality
   type ([[ADR-0069]]).
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
internal models. Its current format layer is `dna_sanger::abif`; support is
determined by the ABIF container and required sequencing tags rather than a
filename suffix. `PLOC.2` is decoded at that format boundary and projected to
canonical Sanger `locus_positions`; downstream scientific modules do not depend
on the ABIF tag name. It does not derive `results/*` paths, validate overwrite
targets, select log paths, or publish outputs. Those remain outer operation/publication concerns.

Operational stage logging crosses this boundary as `tracing` instrumentation
([[ADR-0061]]). Shared
reference-free read processing and reference-guided Variant Analysis emit
structured events and one span per stage without knowing any destination; log
destination selection, record rendering, terminal error logging, and
synchronization stay in the outer operation layer.

The production `pipeline` composes scientific capabilities but does not own
their implementations. Reference-free read processing is shared Sanger science
in `dna_sanger::read_processing`; reference-guided one-read observation is owned by
`variant_analysis` and is reused by CLI analysis, sample evidence, and the
public Rust API.

The first public Rust capability is `variant_analysis::analyze_sanger`. It reuses the same scientific stages as the CLI but returns stable typed called-variant evidence directly and does not open logs or publish JSON. Optional post-calling normalization and nomenclature belong to separate capabilities; see [[architecture/variant-lifecycle|variant lifecycle]]. The exact current contracts are owned by [[reference/rust-api|Rust public API]].

`sample_analysis::analyze_sanger` is the sample-level capability: one placement
pass over a sample's traces returns the typed sample evidence and consensus. The
`sample` and `consensus` commands run the same pass. `SampleAnalysis::call`
completes the recommended `consensus` then `call` workflow in process, through
the same stages as the `call` command, without logs or publication; the
evaluation harness (`dna-eval`, ADR-0075) is its first consumer. Its reads run on worker
threads; a `ReadScope` carries their instrumentation either into the
operation log in trace order or to the caller's dispatcher, so the library still
selects no destination.

Exact syntax and serialized shapes are owned by
[reference](../reference/README.md). Module responsibilities are colocated under
[src](../../src/README.md). Cross-cutting rules that must hold across interfaces
are indexed under [invariants](invariants/README.md).

Commodity bioinformatics machinery may be supplied by maintained external libraries behind these boundaries. DNA-owned canonical contracts remain independent of dependency-specific representations; see [[engineering/dependencies|dependency policy]] and [[ADR-0059]].

Architecture defines where responsibilities belong; it does not duplicate
algorithm details from [design](../design/README.md).
