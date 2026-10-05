# Variant Nomenclature

This module owns the optional target-specific representation boundary after
variant calling and, when selected by a workflow, normalization.

The public `NomenclatureInput` view over `VariantNormalizationResult` carries
reference identity, exact source variants, reconstructed alternate haplotype,
and normalized variants without copying or rewriting them.

Current production behavior implements one deliberately narrow human-mtDNA rule:
`mtdna::apply_hv2_polyc` canonicalizes the rCRS 303-315 HVS-II poly-C window
around T310 into run-length changes at the 309 and 315 boundaries while proving
that the complete alternate haplotype is unchanged.

Target-independent edit/application/render mechanics live in the crate-internal
`variant_representation` module; this module owns only nomenclature policy and
its typed error boundary.

It does **not** currently implement the 513-524 HVS-III AC repeat, HVS-I
16189/16193 policy, Sanger-specific repeat artifact interpretation, decimal
notation rendering, sample reconciliation, or NGS behavior.

Nomenclature must never reinterpret sequencing signal, change caller
eligibility, manufacture phase, or change the represented biological haplotype.
