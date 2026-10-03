# ABIF Format

Owns the current Sanger source-format implementation for Applied Biosystems
Genetic Analysis Data File (ABIF) containers.

Responsibilities:

- bounds-checked big-endian reads over untrusted bytes;
- generic ABIF directory parsing and exact tag lookup;
- Sanger sequencing tag validation and decoding;
- projection into the canonical `model::trace::Chromatogram` evidence model.

The format layer does not own base re-calling, quality control, alignment,
variant calling, CLI paths, logging, or result publication.

Filename extension is not the format contract. Sanger sequencing sample files
commonly use `.ab1`, while support is determined from the ABIF signature and
required sequencing tags.

This implementation is currently retained after reuse evaluation. DNA's ABIF
boundary accepts untrusted bytes and requires checked arithmetic, bounded
directory/payload access, duplicate-tag rejection, exact Sanger tag/layout
validation, and projection into DNA-owned `Chromatogram` evidence. The evaluated
Rust `bio_files::ab1` implementation is Biopython-derived and currently
documents unresolved offset handling, so it is not an equivalent replacement
for this trust-boundary contract. This decision can be revisited when an
implementation satisfies the same requirements and passes differential,
malformed-input, and real-trace corpus validation.

See [ABIF decoding method](../../../../docs/design/abif-decoding.md) and
[input requirements](../../../../docs/requirements/input.md).
