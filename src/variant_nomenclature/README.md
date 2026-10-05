# Variant Nomenclature

This module owns the optional target-specific representation boundary after
variant calling and, when selected by a workflow, normalization.

Current production scope is deliberately small: it exposes a borrowed
`NomenclatureInput` view over `VariantNormalizationResult` containing the
reference identity, exact source variants, reconstructed alternate haplotype,
and normalized variants.

It does **not** currently implement human-mtDNA special-region rules such as
309/315, 513-524, 16189, or 16193. Those rules require their own requirements,
tests, and scientific validation before they become production behavior.

Nomenclature must never reinterpret sequencing signal, change caller
eligibility, manufacture phase, or change the represented biological haplotype.
