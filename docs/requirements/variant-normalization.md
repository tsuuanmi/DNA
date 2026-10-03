# Variant Normalization Requirements

**Requirement namespace:** `SRS-VN-*`

These requirements are part of the canonical [DNA SRS](SRS.md).

- **SRS-VN-001:** Post-calling variant normalization MUST be an optional capability consuming an evidence-backed `CalledVariantSet`, the matching reference, and an explicit normalization policy. It MUST NOT reinterpret source sequencing signal or change caller eligibility.
- **SRS-VN-002:** Normalization MUST verify that the supplied reference name and sequence identity match the reference identity carried by the called-variant set before transforming any representation.
- **SRS-VN-003:** The source called variants MUST be converted to sequence edits and applied jointly to reconstruct one alternate haplotype before any positional movement is accepted.
- **SRS-VN-004:** A normalized edit set MUST reconstruct exactly the same alternate haplotype as the source edit set. A candidate movement that changes the complete reconstructed sequence MUST be rejected.
- **SRS-VN-005:** The normalization result MUST preserve the exact source variant representation, the reconstructed alternate sequence, and the selected normalized representation; normalization MUST NOT destructively overwrite the source description.
- **SRS-VN-006:** The human-mtDNA right-alignment policy MUST shift only sequence-equivalent pure insertion/deletion edits toward the 3'/right-most coordinate and MUST NOT rotate an event across the fixed FASTA/rCRS coordinate seam.
- **SRS-VN-007:** Multiple nearby edits MUST be normalized against the complete reconstructed haplotype rather than shifted independently when an independent move would change phase or sequence.
- **SRS-VN-008:** SNVs MUST remain positionally unchanged by the mtDNA right-alignment policy. Invalid alleles, reference mismatches, unsupported overlapping source edits, and unsupported replacement edits MUST return typed failures rather than be guessed or silently rewritten.
- **SRS-VN-009:** Normalized variants MUST be emitted deterministically in reference/position/allele order.
- **SRS-VN-010:** Variant normalization MUST NOT apply mtDNA nomenclature-window policy, HGVS formatting, VCF interchange normalization, sample reconciliation, genotype interpretation, or clinical interpretation.
