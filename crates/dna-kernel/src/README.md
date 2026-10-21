# dna-kernel

Owns the contracts every plugin family shares
([ADR-0069](../../../docs/decisions/adr/0069-plugin-first-modality-core-post-calling.md)).
It depends on no other DNA crate.

Output contract: none; the kernel runs no step. The facade serializes its
`ReadEvidence` as [`dna.read_evidence/v1`](../../../docs/reference/schemas/read-evidence-v1.schema.json)
([step contracts](../../../docs/architecture/data-flow.md#step-contracts)).

- [error](error/README.md) — the one typed error of every stage, re-exported as `dna::error`.
- `read_evidence.rs` — the modality → core per-read evidence contract, with a
  modality's phase mixtures (`PhaseMixture`, ADR-0082), unreliable spans
  (`UnreliableSpan`, SRS-CALL-024), and tract decodes (`TractLengths`,
  SRS-CALL-025), each validated as the modality attaches it.
- `variant.rs` — public called-variant contracts, re-exported as `dna::variant`.
- `plugin.rs` — plugin descriptor types and their compile-time validation.
- [profile](profile/README.md) — target profiles.
- [reference](reference/README.md) — FASTA loading and reference identity.
- [model](model/README.md) — shared vocabulary: nucleotides, references, and
  the one-base move of 3'/right-most indel placement.
- `nnls.rs` — the one Lawson-Hanson non-negative least squares solver, shared
  by the tract decode (`dna-trace`) and the consensus's two-strand run-length
  rule (`dna-consensus`).
- `checksum.rs` — SHA-256 identities; `bounds.rs` — configuration range checks.

The `test-support` feature exposes evidence and profile fixtures to the tests
of the consuming crates.

See [crates](../../README.md) and the
[architecture overview](../../../docs/architecture/overview.md).
