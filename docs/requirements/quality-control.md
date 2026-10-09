# Quality Control Requirements

**Requirement namespace:** `SRS-QC-*`

These requirements are part of the canonical [DNA SRS](SRS.md).

- **SRS-QC-001:** Each call MUST receive a documented ambiguity/spacing penalty over `quality_control.penalty_window_size` calls and a bounded relative score.
- **SRS-QC-002:** The relative score MUST be explicitly uncalibrated. PCON MUST remain separate internally and apply only when PBAS agrees with the re-called primary; compact JSON MUST omit both the calibration flag and vendor evidence.
- **SRS-QC-003:** A zero maximum penalty MUST produce maximum relative scores without division by zero.
- **SRS-QC-004:** The trim interval MUST be the read's callable span (SRS-CALL-005), widened on each side by at most `variant_calling.read_end_margin` calls of an adjacent `dephased` segment and by no call of any other masked state; it MUST be one half-open interval, and masked segments inside it MUST stay in the retained sequence. The trim interval, retained sequence, penalty, score, and vendor quality applicability MUST be retained internally. Compact JSON MUST expose only the global trim interval and each variant-associated call's uncalibrated relative score under the public field `quality`; sequences, penalties, calibration flags, and vendor quality MUST be omitted.
- **SRS-QC-005:** A callable span that does not fit the read's calls MUST fail with a typed QC error; a read with too few callable calls fails earlier with a typed callability error (SRS-CALL-009).
