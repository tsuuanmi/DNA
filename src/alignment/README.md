# Alignment

Owns bounded, deterministic affine-gap alignment, strand selection, traceback,
and canonical repeat-equivalent gap placement.

Entry point: `align_best` from `mod.rs`.

Key children: `scoring.rs`, `gotoh.rs`, `traceback.rs`, `canonical.rs`,
and `orient.rs`.

This module does not extract variants or mutate upstream signal evidence.

See [alignment requirements](../../docs/srs/alignment.md),
[alignment method](../../docs/methods/alignment.md), and
[alignment invariants](../../docs/architecture/invariants/alignment.md).
