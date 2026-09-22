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
- [model](model/README.md) — validated domain vocabulary.
- [pipeline](pipeline/README.md) — end-to-end command orchestration.
- [quality_control](quality_control/README.md) — relative quality and trimming.
- [reference](reference/README.md) — FASTA loading and identity.
- [report](report/README.md) — contract projection, serialization, publication.
- [sample](sample/README.md) — multi-read evidence aggregation.
- [signal_processing](signal_processing/README.md) — observation-only signal analysis.
- [trace](trace/README.md) — bounded ABIF/AB1 decoding.
- [variant_calling](variant_calling/README.md) — normalized primary-sequence differences.

File-only modules such as `checksum.rs`, `locus.rs`, and `logger.rs` use
rustdoc/source comments. Do not create directories solely to attach README files.

## Dependency rule

Dependencies should point toward shared domain/config/error boundaries rather than
forming cycles. Cross-cutting invariants are canonical in
[docs/architecture/invariants](../docs/architecture/invariants/README.md).

When a directory's responsibility or dependency boundary changes, update its
README in the same change.
