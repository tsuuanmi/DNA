# Basecalling

Owns signal-derived primary and ambiguity calls at validated canonical Sanger loci.

Entry point: `call` from `mod.rs`.

Key children: `peak.rs`, `iupac.rs`, and `call.rs`.

This module does not trim reads, align to a reference, or call variants.

`config.rs` owns the `[basecalling]` configuration section: its raw record, validated
record, and rules (ADR-0069).

See [basecalling requirements](../../docs/requirements/basecalling.md) and
[basecalling method](../../docs/design/basecalling.md).
