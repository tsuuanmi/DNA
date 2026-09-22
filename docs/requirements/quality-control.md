# Quality Control Requirements

**Requirement namespace:** `SRS-QC-*`

These requirements are part of the canonical [DNA SRS](SRS.md).

- **SRS-QC-001:** Each call MUST receive a documented ambiguity/spacing penalty and bounded relative score.
- **SRS-QC-002:** The relative score MUST be explicitly uncalibrated. PCON MUST remain separate internally and apply only when PBAS agrees with the re-called primary; compact JSON MUST omit both the calibration flag and vendor evidence.
- **SRS-QC-003:** A zero maximum penalty MUST produce maximum relative scores without division by zero.
- **SRS-QC-004:** Trimming MUST remove only left/right tails and internally retain the best section, half-open trim interval, retained sequence, penalty, score, and vendor quality applicability. Compact JSON MUST expose only the global trim interval and each variant-associated call's uncalibrated relative score under the public field `quality`; sequences, penalties, calibration flags, and vendor quality MUST be omitted.
- **SRS-QC-005:** Fewer than configured minimum retained bases MUST fail with a typed QC error.
