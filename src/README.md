# Rust Source

This directory contains the production Rust implementation of DNA.

Use this file as the implementation router. Canonical requirements, architecture,
methods, contracts, and decisions live under [docs](../docs/README.md); detailed
API documentation lives in rustdoc and source comments.

## Module map

- [alignment](alignment/README.md) — deterministic evidence-profile alignment.
- [basecalling](basecalling/README.md) — signal-derived base re-calling.
- [cli](cli/README.md) — command-line syntax and typed arguments.
- [config](config/README.md) — strict configuration loading and validation.
- [error](error/README.md) — typed application failures.
- [input](input/README.md) — source-specific sequencing input adapters.
- [model](model/README.md) — validated domain vocabulary.
- [pipeline](pipeline/README.md) — end-to-end operation orchestration.
- [quality_control](quality_control/README.md) — relative quality and trimming.
- [reference](reference/README.md) — FASTA loading and identity.
- [report](report/README.md) — contract projection, serialization, publication.
- [sample](sample/README.md) — multi-read evidence aggregation.
- [signal_processing](signal_processing/README.md) — observation-only signal analysis.
- [variant_calling](variant_calling/README.md) — normalized primary-sequence differences.
- [variant_analysis](variant_analysis/README.md) — public typed raw-to-variant capability.

File-only modules such as `checksum.rs`, `locus.rs`, `logger.rs`, and
`read_processing.rs` use rustdoc/source comments. Do not create directories solely to attach README files.

## Dependency rule

Dependencies point inward toward operation/scientific boundaries and shared
domain/config/error types rather than outward toward frontends. In particular,
the CLI may call the pipeline boundary, but pipeline and scientific modules must
not depend on CLI/`clap` argument types.

Scientific stages that emit operational progress depend on the minimal internal
`StageLog` capability rather than the concrete file-backed `Logger`. Log-path
selection, terminal error records, and synchronization remain outer operation
concerns.

Source-specific filesystem loading belongs to `input`; the current Sanger
adapter produces validated trace/reference/configuration models without knowing
CLI publication paths. Reference-free read processing is a shared scientific
core used by both basecall and Variant Analysis. Reference-guided observation
ownership belongs to `variant_analysis`; `pipeline` consumes these
capabilities rather than owning their implementations.

Cross-cutting invariants are canonical in
[docs/architecture/invariants](../docs/architecture/invariants/README.md).

When a directory's responsibility or dependency boundary changes, update its
README in the same change.
