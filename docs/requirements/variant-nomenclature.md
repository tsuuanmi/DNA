# Variant Nomenclature Requirements

**Requirement namespace:** `SRS-NOM-*`

These requirements are part of the canonical [DNA SRS](SRS.md).

- **SRS-NOM-001:** Target-specific variant nomenclature MUST be an optional post-calling representation capability. It MUST NOT reinterpret sequencing signal, change caller eligibility, manufacture phase, or change the represented biological haplotype.
- **SRS-NOM-002:** A nomenclature operation consuming a normalization result MUST verify that the supplied reference name and sequence identity match the reference identity carried by that result.
- **SRS-NOM-003:** The represented variant set MUST reconstruct exactly the same complete alternate sequence as the normalized variant set. Source variants, normalized variants, reference identity, and the reconstructed alternate sequence MUST remain available in the result.
- **SRS-NOM-004:** The initial human-mtDNA HVS-II policy applies only to the validated rCRS 303-315 poly-C window `CCCCCCCTCCCCC`, whose single `T` anchor is rCRS position 310.
- **SRS-NOM-005:** When the HVS-II alternate window contains only `C` bases around exactly one `T` anchor, sequence-equivalent movement of that anchor MUST be represented as C-run length change at the 309 and 315 boundaries rather than as movement of the anchor base itself.
- **SRS-NOM-006:** Multiple inserted C bases at the left HVS-II run MUST remain anchored after position 309 and multiple inserted C bases at the right run MUST remain anchored after position 315 so an outer notation layer can render ordered forms such as `309.1C`, `309.2C`, and `315.1C`.
- **SRS-NOM-007:** Variants outside the HVS-II nomenclature window MUST be preserved unchanged. An edit crossing the validated window boundary MUST fail explicitly rather than be partially rewritten.
- **SRS-NOM-008:** If the HVS-II window does not match the rule's validated anchored-homopolymer pattern, the rule MUST preserve the normalized representation rather than invent another decomposition, except when the supplied reference itself does not contain the validated rCRS window, which MUST fail explicitly.
- **SRS-NOM-009:** The HVS-II policy MUST NOT implement HVS-III 513-524, HVS-I 16189/16193, Sanger repeat-artifact interpretation, decimal notation rendering, sample reconciliation, or NGS-specific behavior.
