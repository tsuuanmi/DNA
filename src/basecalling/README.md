# Basecalling

Owns signal-derived primary and ambiguity calls at validated ABIF PLOC loci.

Entry point: `call` from `mod.rs`.

Key children: `peak.rs`, `iupac.rs`, and `call.rs`.

This module does not trim reads, align to a reference, or call variants.

See [basecalling requirements](../../docs/srs/basecalling.md) and
[basecalling method](../../docs/methods/basecalling.md).
