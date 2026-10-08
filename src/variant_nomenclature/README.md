# Variant Nomenclature

This module owns the optional target-specific representation boundary after
variant calling and, when selected by a workflow, normalization.

The public `NomenclatureInput` view over `VariantNormalizationResult` carries
reference identity, exact source variants, reconstructed alternate haplotype,
and normalized variants without copying or rewriting them.

Current production behavior implements the human-mtDNA control-region policy,
`mtdna::apply_control_region`: declarative rCRS windows (HVS-II 303-315 and
HVS-I 16181-16193 poly-C, HVS-III 513-524 AC repeat), each with ordered
representation rules, while proving that the complete alternate haplotype is
unchanged.

Target-independent edit/application/render mechanics live in the crate-internal
`variant_representation` module; this module owns only nomenclature policy.
Its failure vocabulary, `NomenclatureError`, lives in `error`.

It does **not** implement Sanger-specific repeat artifact interpretation,
primer callable ranges, sample reconciliation, or NGS behavior. Decimal
notation rendering lives in `report::notation`, and the `sample` workflow
composes this policy per read.

Nomenclature must never reinterpret sequencing signal, change caller
eligibility, manufacture phase, or change the represented biological haplotype.

See [variant lifecycle](../../docs/architecture/variant-lifecycle.md),
[ADR-0060](../../docs/decisions/adr/0060-separate-variant-canonicalization-nomenclature.md),
and the [variant-nomenclature design](../../docs/design/variant-nomenclature.md).
