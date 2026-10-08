# Variant Calling

Owns extraction, direct evidence mapping, allele anchoring, and configured
eligibility of primary-sequence SNVs and supported small indels.

Key children: `extract.rs`, `mapping.rs`, `anchor.rs`, and `filter.rs`.
`anchor.rs` constructs and validates anchored REF/ALT alleles while preserving
the alignment-selected event placement; post-calling haplotype normalization
and mtDNA nomenclature belong to `variant_normalization` and
`variant_nomenclature`.

This module does not infer genotype, heteroplasmy, phase, pathogenicity, or
clinical significance.

See [variant requirements](../../docs/requirements/variants.md),
[variant-calling method](../../docs/design/variant-calling.md), and
[scientific-state invariants](../../docs/architecture/invariants/scientific-state.md).
