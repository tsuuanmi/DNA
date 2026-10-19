---
decided-by: [ADR-0060]
---

# Variant Normalization Requirements

**Requirement namespace:** `SRS-VN-*`

These requirements are part of the canonical [[requirements/SRS|DNA SRS]].

- **SRS-VN-001:** Post-calling variant normalization MUST be an optional capability consuming an evidence-backed `CalledVariantSet`, the matching reference, and an explicit normalization policy. It MUST NOT reinterpret source sequencing signal or change caller eligibility.
- **SRS-VN-002:** Normalization MUST verify that the supplied reference name and sequence identity match the reference identity carried by the called-variant set before transforming any representation.
- **SRS-VN-003:** The source called variants MUST be converted to sequence edits and applied jointly to reconstruct one alternate haplotype before any positional movement is accepted.
- **SRS-VN-004:** A normalized edit set MUST reconstruct exactly the same alternate haplotype as the source edit set. A candidate movement that changes the complete reconstructed sequence MUST be rejected.
- **SRS-VN-005:** The normalization result MUST preserve the exact source variant representation, the reconstructed alternate sequence, and the selected normalized representation; normalization MUST NOT destructively overwrite the source description.
- **SRS-VN-006:** The right-alignment policy (`RightAligned`, profile `indel_placement = "right"`) MUST shift only sequence-equivalent pure insertion/deletion edits toward the 3'/right-most coordinate and MUST NOT rotate an event across the fixed FASTA coordinate seam.
- **SRS-VN-007:** Multiple nearby edits MUST be normalized against the complete reconstructed haplotype rather than shifted independently when an independent move would change phase or sequence.
- **SRS-VN-008:** SNVs MUST remain positionally unchanged by the right-alignment policy, except in a deletion trade (variant-normalization design): a deletion that cannot shift further, including one whose shift would overlap a substitution, MAY move 3' by at most the policy's `max_trade` bases (the profile's `normalization.max_trade`) across substitutions only, when the same haplotype needs no more edits and, where the caller gives trade bounds (the sample, call, and notation workflows do), every base of the stretch it crosses is read and the stretch crosses no boundary of a nomenclature window; the public `normalize` gives none; the 3'-most such placement is chosen and the substitutions it passes are rewritten. Invalid alleles, reference mismatches, unsupported overlapping source edits, and unsupported replacement edits MUST return typed failures rather than be guessed or silently rewritten.
- **SRS-VN-009:** Normalized variants MUST be emitted deterministically in reference/position/allele order.
- **SRS-VN-010:** Variant normalization MUST NOT apply nomenclature-window policy, HGVS formatting, VCF interchange normalization, sample reconciliation, genotype interpretation, or clinical interpretation.
- **SRS-VN-011:** A normalization result intended for downstream representation policy MUST retain the exact reference identity alongside source variants, reconstructed alternate sequence, and normalized variants.

## Implemented in

<!-- docs:begin backlinks -->
- **Implemented in:** [design/variant-normalization](../design/variant-normalization.md), [reference/profiles](../reference/profiles.md), [reference/rust-api](../reference/rust-api.md)
<!-- docs:end backlinks -->
