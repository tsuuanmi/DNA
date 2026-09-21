# Evidence Invariants

These invariants are part of the canonical [system invariant set](../invariants.md).

- **INV-EVID-001:** Decoded analyzed A/C/G/T channel arrays are immutable source evidence after validation.
- **INV-EVID-002:** Derived signal features, corrected waveforms, quality values, alignments, and variants never overwrite the source evidence from which they were derived.
- **INV-EVID-003:** Vendor PBAS/PCON are vendor evidence, not authoritative DNA output.
- **INV-EVID-004:** `LocusEvidence` and `EvidenceProfile` are derived directly from immutable analyzed A/C/G/T channel values and PLOC-defined geometry; primary/ambiguity calls, selected basecall peaks, and qualifying-channel thresholds MUST NOT determine profile membership or weights.
- **INV-EVID-005:** A zero-positive-signal locus has no evidence profile; the system MUST NOT synthesize a uniform or reference-guided profile as a fallback.
- **INV-EVID-006:** Reference placement MAY consume `EvidenceProfile`, but reference context MUST NOT mutate or rewrite upstream locus evidence or base calls.
- **INV-EVID-007:** PLOC validity and optional vendor-series cardinality are distinct evidence dimensions. Valid PLOC loci remain the authoritative current-method event anchors; PBAS/PCON length mismatch cannot add/remove DNA loci and must remain explicit non-fatal integrity evidence.
- **INV-EVID-008:** Exact signed-16-bit clipping and whole-trace event-signal scale observations are evidence only. They cannot mutate channels, calls, trim bounds, alignment, or variants, and an amplitude ratio cannot become an artifact/dye-blob label without a separately specified method.
- **INV-EVID-009:** A locus evidence event MUST remain anchored to the PLOC-local chromatogram event: among positive local maxima of total corrected A/C/G/T signal inside the locus window, the nearest to PLOC wins before amplitude. A stronger but more distant neighboring event cannot steal the locus profile. If no positive local event exists, PLOC itself is used.
