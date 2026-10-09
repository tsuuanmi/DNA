# dna-kernel

Owns the contracts every plugin family shares
([ADR-0069](../../../docs/decisions/adr/0069-plugin-first-modality-core-post-calling.md)).
It depends on no other DNA crate.

- [error](error/README.md) — the one typed error of every stage, re-exported as `dna::error`.
- `read_evidence.rs` — the modality → core per-read evidence contract.
- `variant.rs` — public called-variant contracts, re-exported as `dna::variant`.
- `plugin.rs` — plugin descriptor types and their compile-time validation.
- [profile](profile/README.md) — target profiles.
- [reference](reference/README.md) — FASTA loading and reference identity.
- [model](model/README.md) — shared vocabulary: nucleotides and references.
- `checksum.rs` — SHA-256 identities; `bounds.rs` — configuration range checks.

The `test-support` feature exposes evidence and profile fixtures to the tests
of the consuming crates.

See [crates](../../README.md) and the
[architecture overview](../../../docs/architecture/overview.md).
