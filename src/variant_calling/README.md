# Variant Calling

Owns extraction, direct evidence mapping, allele anchoring, and configured
eligibility of primary-sequence SNVs and supported small indels. Reportable
regions come from the target profile and are passed in by the caller.

Key children: `extract.rs`, `mapping.rs`, `anchor.rs`, `callability.rs`, and
`filter.rs`. `callability.rs` derives the calls of a read that cannot support a
variant (read ends and post-homopolymer windows, ADR-0062) from the read's own
calls in sequencing order.
`anchor.rs` constructs and validates anchored REF/ALT alleles while preserving
the alignment-selected event placement; post-calling haplotype normalization
and target nomenclature belong to `variant_normalization` and
`variant_nomenclature`.

This module does not infer genotype, heteroplasmy, phase, pathogenicity, or
clinical significance.

See [variant requirements](../../docs/requirements/variants.md),
[variant-calling method](../../docs/design/variant-calling.md), and
[scientific-state invariants](../../docs/architecture/invariants/scientific-state.md).
