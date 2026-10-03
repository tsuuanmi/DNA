# Variant Calling

Owns extraction, direct evidence mapping, allele anchoring, and configured
eligibility of primary-sequence SNVs and supported small indels.

Key children: `extract.rs`, `mapping.rs`, `normalize.rs`, and `filter.rs`.
`normalize.rs` currently constructs/validates anchored REF/ALT alleles while
preserving the alignment-selected event placement; it is not the future
post-calling haplotype canonicalization or mtDNA nomenclature layer.

This module does not infer genotype, heteroplasmy, phase, pathogenicity, or
clinical significance.

See [variant requirements](../../docs/requirements/variants.md),
[variant-calling method](../../docs/design/variant-calling.md), and
[scientific-state invariants](../../docs/architecture/invariants/scientific-state.md).
