# Variant Nomenclature

This module owns the optional target-specific representation boundary after
variant calling and, when selected by a workflow, normalization.

The public `NomenclatureInput` view over `VariantNormalizationResult` carries
reference identity, exact source variants, reconstructed alternate haplotype,
and normalized variants without copying or rewriting them.

`apply` runs the nomenclature windows of an explicit target profile
(`crate::profile`): within each window the first rule whose candidate
reproduces the window haplotype wins, and the complete alternate haplotype is
proven unchanged. `windows.rs` owns the generic engine and the rule
implementations; no target windows, bases, or references are hard-coded here.
The shipped human-mtDNA profile declares the HVS-II 303-315 and HVS-I
16181-16193 poly-C windows and the HVS-III 513-524 AC repeat.

Target-independent edit/application/render mechanics live in the crate-internal
`variant_representation` module; this module owns only nomenclature policy.
Its failure vocabulary, `NomenclatureError`, lives in `error`.

It does **not** implement Sanger-specific repeat artifact interpretation,
primer callable ranges, sample reconciliation, or NGS behavior. Decimal
notation rendering lives in `report::notation`, and the `sample` workflow
composes it per read when the profile declares notation.

Nomenclature must never reinterpret sequencing signal, change caller
eligibility, manufacture phase, or change the represented biological haplotype.

See [variant lifecycle](../../docs/architecture/variant-lifecycle.md),
[ADR-0060](../../docs/decisions/adr/0060-separate-variant-canonicalization-nomenclature.md),
and the [variant-nomenclature design](../../docs/design/variant-nomenclature.md).
