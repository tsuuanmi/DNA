# dna-core

Owns the core caller
([ADR-0069](../../../docs/decisions/adr/0069-plugin-first-modality-core-post-calling.md)).
It works only on modality-neutral `ReadEvidence` and depends only on
`dna-kernel`.

- [alignment](alignment/README.md) — deterministic evidence-profile alignment.
- [variant_calling](variant_calling/README.md) — differences, mapping, anchoring, and core eligibility gates.
- [sample](sample/README.md) — multi-read evidence aggregation over `CalledRead` records.
- [consensus](consensus/README.md) — adjudication between a sample's reads into consensus segments (ADR-0073).
- `read_call.rs` — the one-read path from `ReadEvidence` to a `CalledRead`, `CoreConfig`, and the core plugin descriptor.
- [model](model/README.md) — records the core produces.

See [crates](../../README.md).
