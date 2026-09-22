# Validation Acceptance Criteria

**Requirement namespace:** `SRS-VAL-*`

These normative criteria define the minimum validation evidence required by the canonical [DNA SRS](../requirements/SRS.md).

- **SRS-VAL-001:** Parser, calling, signal processing, QC, alignment, normalization, JSON, and atomic publication MUST have focused boundary/adversarial tests.
- **SRS-VAL-002:** A synthetic canonical ABIF MUST exercise end-to-end forward/reverse and variant behavior without identifying data.
- **SRS-VAL-003:** Real-trace release evidence MUST follow [data policy](../governance/data.md); ignored local data MUST never be a build/test prerequisite.
- **SRS-VAL-004:** Rust source-policy validation, documentation-structure validation, format, check, Clippy warnings-denied, all tests, rustdoc, schema/example validation, TOML validation, source-directory README coverage, and rCRS identity gates MUST pass.
