# Configuration Requirements

**Requirement namespace:** `SRS-CFG-*`

These requirements are part of the canonical [DNA SRS](SRS.md).

- **SRS-CFG-001:** DNA MUST load `DNA_CONFIG` or `config/dna.toml`; it MUST NOT parse `.env`.
- **SRS-CFG-002:** Schema version MUST be integer `7`; all documented sections/keys are required, including a non-empty root `profile` naming the target profile (SRS-PRF-001), resolved against the configuration file's directory when relative.
- **SRS-CFG-003:** Unknown/duplicate keys, missing fields, non-finite values, unsupported enums/versions, and invalid ranges MUST fail.
- **SRS-CFG-004:** Environment MUST NOT override individual scientific settings. `DNA_LOG_DIR` MAY select only the operational log directory.
- **SRS-CFG-005:** JSON MUST include the configuration checksum but MUST omit method constants, the local config path, and the effective-value expansion. Effective values and configuration schema version 7 remain in the selected strict TOML; target knowledge is identified separately by the profile identity (SRS-PRF-007).
- **SRS-CFG-006:** Hard caps MUST bound AB1 bytes, reference length, alignment cells, indel length, and reachable peak thresholds before unsafe allocation.
- **SRS-CFG-007:** Configuration MUST require a positive `sample_reconciliation.minimum_comparable_bases` and finite `sample_reconciliation.minimum_overlap_agreement` in `(0,1]`; these settings apply only after independently placed reads enter sample reconciliation.
- **SRS-CFG-008:** Configuration MUST require a `[callability]` section with `window_calls` in `8..=64`, finite `onset_defect_fraction` in `(0,1]`, finite `exit_defect_fraction` in `[0,1]` below the onset fraction, finite `minimum_main_share` in `(0,1]`, finite `maximum_far_share` in `[0,1]`, finite `minimum_shadow_share` in `(0,1]`, finite `weak_amplitude_fraction` in `[0,0.5]`, `repeat_min_length` of at least `2`, and a positive `minimum_callable_calls`; these settings control read callability, which the trim interval, alignment, eligibility, and sample aggregation consume (SRS-CALL-009..012).
