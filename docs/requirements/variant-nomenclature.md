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
- **SRS-NOM-010:** When the reference sequence is the rCRS (normalized-sequence SHA-256 `f156ff3f65bbcc80c7ebb9936dceb96b1477b4f8f535c4e1dbe7baea225cbc66`), the `sample` workflow MUST represent each read's eligible variants independently: human-mtDNA right alignment followed by the HVS-II rule. A read with an edit crossing the HVS-II window MUST keep its right-aligned representation. Against any other reference the workflow MUST NOT produce notation.
- **SRS-NOM-011:** Notation MUST render represented variants one call per changed base: a substitution as `<position><base>`, each deleted base as `<position>DEL`, and each inserted base as `<anchor>.<ordinal><base>` after the preceding reference base (anchor `0` before position 1), ordered by position and then ordinal.
- **SRS-NOM-012:** Sample notation MUST list each distinct rendered call with every read whose represented variants contain it, in read-registry order. It MUST NOT include ineligible observations, infer reference support from a read's silence, or adjudicate disagreement between reads.
