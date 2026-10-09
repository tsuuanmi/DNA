# Rust Source

This directory contains the production Rust implementation of DNA.

Use this file as the implementation router. Canonical requirements, architecture,
methods, contracts, and decisions live under [docs](../docs/README.md); detailed
API documentation lives in rustdoc and source comments.

## Module map

- [alignment](alignment/README.md) — deterministic evidence-profile alignment.
- [basecalling](basecalling/README.md) — signal-derived base re-calling.
- [callability](callability/README.md) — signal-derived per-read callability: phase state, typed mask, callable span.
- [cli](cli/README.md) — command-line syntax and typed arguments.
- [config](config/README.md) — strict configuration loading and validation.
- [error](error/README.md) — typed application failures.
- [input](input/README.md) — source-specific sequencing input adapters.
- [model](model/README.md) — validated domain vocabulary.
- [pipeline](pipeline/README.md) — end-to-end operation orchestration.
- [profile](profile/README.md) — target profiles: reference identity, regions, representation chain.
- [quality_control](quality_control/README.md) — relative quality and trimming.
- `plugin.rs` — the static plugin registry and workflow compositions, validated at compile time (ADR-0069).
- `read_evidence.rs` — the modality → core per-read evidence contract (ADR-0069).
- [read_processing](read_processing/README.md) — shared reference-free Sanger read processing and the Sanger evidence adapter.
- [reference](reference/README.md) — FASTA loading and identity.
- [report](report/README.md) — contract projection, serialization, publication.
- [sample](sample/README.md) — multi-read evidence aggregation.
- [signal_processing](signal_processing/README.md) — observation-only signal analysis.
- [variant_calling](variant_calling/README.md) — normalized primary-sequence differences.
- `variant.rs` — public canonical called-variant contracts (`dna::variant`).
- [variant_analysis](variant_analysis/README.md) — public typed raw-to-variant capability.
- [variant_nomenclature](variant_nomenclature/README.md) — optional profile-driven target nomenclature.
- `conformance.rs` — reports represented calls that break the profile's notation conventions.
- [variant_normalization](variant_normalization/README.md) — optional haplotype-preserving representation normalization.

File-only modules such as `checksum.rs`, `locus.rs`, `operation_log.rs`,
`conformance.rs`, `plugin.rs`, `read_evidence.rs`, and `variant.rs` use rustdoc/source comments. Do not create directories solely to attach README files.

## Dependency rule

Every top-level module belongs to one layer of
[ADR-0064](../docs/decisions/adr/0064-crate-ready-module-layering.md) — core,
target data, science, adapters/capabilities, delivery — and may depend only on
the same or a lower layer, without cycles;
`tools/python/scripts/validate_module_layers.py` enforces this in CI.
Dependencies point inward toward operation/scientific boundaries and shared
domain/config/error types rather than outward toward frontends. In particular,
the CLI may call the pipeline boundary, but pipeline and scientific modules must
not depend on CLI/`clap` argument types.

The same validator enforces modality neutrality
([ADR-0069](../docs/decisions/adr/0069-plugin-first-modality-core-post-calling.md)):
`plugin`, `read_evidence`, `variant`, `alignment`, `variant_calling`, and the
post-calling modules must not depend on Sanger modules or on Sanger children of
`model`. The core caller sees a read only as `ReadEvidence`, which
`read_processing` builds for Sanger.

Scientific stages report operational progress only through `tracing` events and
per-stage spans; they never select a destination. `operation_log` renders those
events into the per-operation file, and log-path selection, terminal error
records, and synchronization remain outer operation concerns.

Source-specific filesystem loading belongs to `input`; the current Sanger
adapter produces validated trace/reference/configuration models without knowing
CLI publication paths. Reference-free read processing (basecalling, signal processing, callability,
quality control) is a shared scientific core used by both basecall and Variant
Analysis. Reference-guided observation
ownership belongs to `variant_analysis`; `pipeline` consumes these
capabilities rather than owning their implementations.

Cross-cutting invariants are canonical in
[docs/architecture/invariants](../docs/architecture/invariants/README.md).

When a directory's responsibility or dependency boundary changes, update its
README in the same change.
