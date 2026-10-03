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

See [ABIF decoding method](../../../../docs/design/abif-decoding.md) and
[input requirements](../../../../docs/requirements/input.md).
