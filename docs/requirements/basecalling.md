# DNA-derived Base Re-calling Requirements

**Requirement namespace:** `SRS-BC-*`

These requirements are part of the canonical [DNA SRS](SRS.md).

- **SRS-BC-001:** Calls MUST derive from analyzed A/C/G/T signals at validated PLOC loci. PBAS MUST NOT become final algorithm output.
- **SRS-BC-002:** Windows MUST use neighboring PLOC midpoints with bounded first/last extrapolation.
- **SRS-BC-003:** Each channel MUST use its strongest positive local maximum under the documented plateau rule, or its PLOC sample with an explicit fallback marker.
- **SRS-BC-004:** No positive strongest signal or exact strongest tie MUST yield primary N.
- **SRS-BC-005:** For a uniquely strongest positive selected peak, a channel MUST qualify only when both its positive selected peak divided by the primary peak height and its positive signal at the primary peak sample divided by the primary peak height are greater than or equal to `secondary_peak_ratio`. One qualifying channel MUST call its base canonically; two qualifying channels MUST use the standard two-base IUPAC code on the strongest primary; three qualifying channels MUST retain the strongest primary with unresolved `N` ambiguity; four qualifying channels MUST be unresolved `N` for both primary and ambiguity.
- **SRS-BC-006:** Every call MUST internally retain its original index, PLOC, call-window sample bounds, four peak heights/positions/sources, primary, ambiguity, qualifying channels, and vendor agreement. Compact analysis JSON MUST expose each variant-associated call as reference-oriented base plus co-located A/C/G/T primary-event channel heights and uncalibrated quality; selected peak positions/sources and ABIF implementation coordinates remain internal.
- **SRS-BC-007:** No sample-specific poly-C correction or reference-aware basecall rescue is permitted; configured biological regions apply only after variant normalization.
- **SRS-BC-008:** The MVP is a signal-derived **re-caller at vendor-defined loci**, not an independent de novo locus detector. `PLOC.2` defines the call loci for the current method; `PBAS.2` and `PCON.2` remain optional vendor evidence and MUST NOT determine DNA's final base identity.
