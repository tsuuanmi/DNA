# Sanger Input

Owns source-specific loading for Sanger sequencing evidence.

The current source format is **ABIF**. Applied Biosystems sequencing sample files
commonly use the `.ab1` extension, but the adapter validates the ABIF container
signature and required sequencing tags rather than coupling scientific behavior
to a filename extension.

ABIF decoding lives in the Sanger crate, `dna_sanger::abif`, which owns:

- bounds-checked ABIF directory parsing;
- exact sequencing-tag lookup and validation;
- decoding analyzed DATA channels into canonical A/C/G/T order;
- PLOC.2 decoding into canonical Sanger locus positions;
- optional PBAS/PCON vendor evidence.

The output is `dna_sanger::model::sanger::Chromatogram`, the canonical Sanger evidence model.
ABIF directory entries, tag names, and format-specific structures do not cross
that boundary. In particular, `PLOC.2` is projected to canonical `locus_positions`.

This facade module (`mod.rs`) only decodes the traces a command names, after the
facade's `input` module has loaded the configuration, profile, and reference
(ADR-0070).

See [input requirements](../../../docs/requirements/input.md) and
[ABIF decoding method](../../../docs/design/abif-decoding.md).
