# dna-normalize

Owns the post-calling `normalization` and `nomenclature` plugins
([ADR-0069](../../../docs/decisions/adr/0069-plugin-first-modality-core-post-calling.md),
[ADR-0077](../../../docs/decisions/adr/0077-normalize-and-empop-crates.md)),
which consume called variants and write them in the lab's reviewer
convention. It depends only on `dna-kernel`.

- `variant_representation.rs` — haplotype-preserving edit conversion, application, and rendering support.
- [variant_normalization](variant_normalization/README.md) — optional sequence-equivalent normalization.
- [variant_nomenclature](variant_nomenclature/README.md) — profile-driven window representation.

See [crates](../../README.md).
