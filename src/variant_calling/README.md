# Variant Calling

Owns extraction, mapping, normalization, and configured eligibility of
primary-sequence SNVs and supported small indels.

Key children: `extract.rs`, `mapping.rs`, `normalize.rs`, and `filter.rs`.

This module does not infer genotype, heteroplasmy, phase, pathogenicity, or
clinical significance.

See [variant requirements](../../docs/srs/variants.md),
[variant-calling method](../../docs/methods/variant-calling.md), and
[scientific-state invariants](../../docs/architecture/invariants/scientific-state.md).
