# Configuration Requirements

**Requirement namespace:** `SRS-CFG-*`

These requirements are part of the canonical [DNA SRS](SRS.md).

- **SRS-CFG-001:** DNA MUST load `DNA_CONFIG` or `config/dna.toml`; it MUST NOT parse `.env`.
- **SRS-CFG-002:** Schema version MUST be integer `5`; all documented sections/keys are required.
- **SRS-CFG-003:** Unknown/duplicate keys, missing fields, non-finite values, unsupported enums/versions, and invalid ranges MUST fail.
- **SRS-CFG-004:** Environment MUST NOT override individual scientific settings. `DNA_LOG_DIR` MAY select only the operational log directory.
- **SRS-CFG-005:** JSON MUST include the configuration checksum but MUST omit method constants, the local config path, and the effective-value expansion. Effective values and configuration schema version 5 remain in the selected strict TOML.
- **SRS-CFG-006:** Hard caps MUST bound AB1 bytes, reference length, alignment cells, indel length, and reachable peak thresholds before unsafe allocation.
- **SRS-CFG-007:** Configuration MUST require a positive `sample_reconciliation.minimum_comparable_bases` and finite `sample_reconciliation.minimum_overlap_agreement` in `(0,1]`; these settings apply only after independently placed reads enter sample reconciliation.
