# Variant Normalization

Owns optional post-calling representation normalization of evidence-backed
variants.

The module consumes the public `CalledVariantSet` boundary plus an explicit
normalization policy. It reconstructs the complete alternate haplotype, preserves
the exact source variant representation, and may select another
sequence-equivalent representation only when the reconstructed haplotype remains
identical.

Key children: `edit.rs` owns source-variant → minimal-edit conversion and complete haplotype reconstruction; `mtdna.rs` owns the human-mtDNA 3'/right-most policy; `render.rs` converts selected minimal edits back to anchored public variants.

The initial implemented policy is human-mtDNA 3'/right-most indel placement with
the FASTA coordinate boundaries treated as the fixed rCRS seam.

This module does **not** call variants from raw evidence, reinterpret Sanger
signal, perform sample reconciliation, apply mtDNA nomenclature rules, or define
VCF/HGVS serialization conventions.

See [variant lifecycle](../../docs/architecture/variant-lifecycle.md),
[ADR-0060](../../docs/decisions/adr/0060-separate-variant-canonicalization-nomenclature.md),
and the forthcoming/current variant-normalization method documentation.
