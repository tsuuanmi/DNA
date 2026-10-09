# Rust Source

This directory contains the production Rust implementation of DNA.

Use this file as the implementation router. Canonical requirements, architecture,
methods, contracts, and decisions live under [docs](../docs/README.md); detailed
API documentation lives in rustdoc and source comments.

## Module map

This is the `dna` facade crate. The plugin crates live in
[crates](../crates/README.md).

- [cli](cli/README.md) — command-line syntax and typed arguments.
- [config](config/README.md) — the configuration envelope composing the plugin-owned sections.
- [input](input/README.md) — loading orchestration and source adapters: Sanger traces, consensus sequences, variants documents.
- [model](model/README.md) — the Sanger read observation and the serializable result contracts.
- [pipeline](pipeline/README.md) — operation orchestration, the plugin registry, and workflow compositions.
- [report](report/README.md) — contract projection, serialization, publication.
- [variant_analysis](variant_analysis/README.md) — public typed Sanger raw-to-variant capability.
- `operation_log.rs` — append-only operational records rendered from `tracing` events.
- `lib.rs` — the public API: `cli`, `run`, `variant_analysis`, and re-exports of `error`, `variant`, `profile`, `variant_normalization`, and `variant_nomenclature` from the plugin crates.

File-only modules use rustdoc/source comments. Do not create directories solely
to attach README files.

## Dependency rule

The facade depends on every plugin crate; the plugin crates depend only on
`dna-kernel` ([ADR-0070](../docs/decisions/adr/0070-workspace-split-by-plugin-family.md)).
Each crate's module graph is acyclic; `tools/python/scripts/validate_module_layers.py`
enforces both in CI.
Dependencies point inward toward operation/scientific boundaries and shared
domain/config/error types rather than outward toward frontends. In particular,
the CLI may call the pipeline boundary, but pipeline and scientific modules must
not depend on CLI/`clap` argument types.

The core caller sees a read only as `ReadEvidence`, which `read_processing`
builds for Sanger and `input::sequence` for reviewed consensus sequences.

Scientific stages report operational progress only through `tracing` events and
per-stage spans; they never select a destination. `operation_log` renders those
events into the per-operation file, and log-path selection, terminal error
records, and synchronization remain outer operation concerns.

Source-specific filesystem loading belongs to `input`: it loads the
configuration, profile, and reference, and its Sanger adapter decodes traces,
without knowing CLI publication paths. Reference-free read processing (basecalling, signal processing, callability,
quality control) is a shared scientific core used by both basecall and Variant
Analysis. Reference-guided observation
ownership belongs to `variant_analysis`; `pipeline` consumes these
capabilities rather than owning their implementations.

Cross-cutting invariants are canonical in
[docs/architecture/invariants](../docs/architecture/invariants/README.md).

When a directory's responsibility or dependency boundary changes, update its
README in the same change.
