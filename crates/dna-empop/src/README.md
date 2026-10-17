# dna-empop

Owns the post-calling `conformance` and `empop` plugins
([ADR-0069](../../../docs/decisions/adr/0069-plugin-first-modality-core-post-calling.md),
[ADR-0077](../../../docs/decisions/adr/0077-normalize-and-empop-crates.md)):
forensic (EMPOP-style) notation conventions, reported and never rewritten, and
the offline EMPOP query. It formats and decides nothing about the sequence; it
is deterministic and offline, and it depends only on `dna-kernel`.

- `conformance.rs` — reports represented calls that break the profile's notation conventions.
- `export.rs` — the EMPOP query of a sample (SRS-NOM-018): interpretation
  ranges over the positions exactly one record calls, the rendered calls
  inside them (an insertion by its anchor), and why a record position is cut.
  It keeps no network client and no corpus index.

See [crates](../../README.md).
