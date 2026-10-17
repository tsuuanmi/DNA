# dna-empop

Owns the post-calling `conformance` plugin
([ADR-0069](../../../docs/decisions/adr/0069-plugin-first-modality-core-post-calling.md),
[ADR-0077](../../../docs/decisions/adr/0077-normalize-and-empop-crates.md)):
forensic (EMPOP-style) notation conventions, reported and never rewritten. It
is deterministic and offline, and it depends only on `dna-kernel`.

- `conformance.rs` — reports represented calls that break the profile's notation conventions.

See [crates](../../README.md).
