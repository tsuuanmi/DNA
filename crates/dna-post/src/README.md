# dna-post

Owns the post-calling plugins
([ADR-0069](../../../docs/decisions/adr/0069-plugin-first-modality-core-post-calling.md)),
which consume called variants. It depends only on `dna-kernel`.

- `variant_representation.rs` — haplotype-preserving edit conversion, application, and rendering support.
- [variant_normalization](variant_normalization/README.md) — optional sequence-equivalent normalization.
- [variant_nomenclature](variant_nomenclature/README.md) — profile-driven window representation.
- `conformance.rs` — reports represented calls that break the profile's notation conventions.

See [crates](../../README.md).
