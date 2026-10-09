# Validation Acceptance Criteria

**Requirement namespace:** `SRS-VAL-*`

These normative criteria define the minimum validation evidence required by the canonical [DNA SRS](../requirements/SRS.md).

- **SRS-VAL-001:** Parser, calling, signal processing, callability, QC, alignment, normalization, JSON, and atomic publication MUST have focused boundary/adversarial tests.
- **SRS-VAL-002:** A synthetic canonical ABIF MUST exercise end-to-end forward/reverse and variant behavior without identifying data.
- **SRS-VAL-003:** Real-trace release evidence MUST follow [data policy](../governance/data.md); ignored local data MUST never be a build/test prerequisite.
- **SRS-VAL-004:** Rust source-policy validation, documentation-structure validation, format, check, Clippy warnings-denied, all tests, rustdoc, schema/example validation, TOML validation, source-directory README coverage, and rCRS identity gates MUST pass.
- **SRS-VAL-005:** Synthetic ABIF fixtures MUST exercise one- and two-sided slippage shadow ladders behind a long homopolymer, double peaks that the shadow model does not explain, and an amplitude-collapsed tail, and MUST assert the published phase-state segments, repeat attribution, callable span, and masked-call count without identifying data. They MUST also show the mask acting: a dephased tail trimmed away, an internal mixed stretch aligned as unresolved with split callable reference segments and no variant inside it, and an SNV after a clean long homopolymer reported.
- **SRS-VAL-006:** Real-trace evidence quoted for callability thresholds MUST name the DNA revision, be produced outside the repository (ADR-0065), be quoted in aggregate without sample identifiers, and MUST NOT be a build, test, or CI prerequisite.
