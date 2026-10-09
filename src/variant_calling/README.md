# Variant Calling

Owns extraction, direct evidence mapping, allele anchoring, and configured
eligibility of primary-sequence SNVs and supported small indels. Reportable
regions come from the target profile and are passed in by the caller.

Key children: `extract.rs`, `mapping.rs`, `anchor.rs`, `eligibility.rs`, and
`filter.rs`. `eligibility.rs` decides which calls of a read can support a
variant: calls near the ends of the trim interval (`read_end`, ADR-0062) and
masked calls, which carry the reason of their phase segment from the read's
[callability](../callability/README.md) (ADR-0067). A masked call is no SNV
candidate in `extract.rs`.
`anchor.rs` constructs and validates anchored REF/ALT alleles while preserving
the alignment-selected event placement; post-calling haplotype normalization
and target nomenclature belong to `variant_normalization` and
`variant_nomenclature`.

This module does not infer genotype, heteroplasmy, phase, pathogenicity, or
clinical significance.

See [variant requirements](../../docs/requirements/variants.md),
[variant-calling method](../../docs/design/variant-calling.md), and
[scientific-state invariants](../../docs/architecture/invariants/scientific-state.md).
