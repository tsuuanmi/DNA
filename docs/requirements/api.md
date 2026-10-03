# Public API Requirements

**Requirement namespace:** `SRS-API-*`

These requirements are part of the canonical [DNA SRS](SRS.md).

- **SRS-API-001:** DNA MUST expose a Rust library capability for reference-guided Sanger Variant Analysis without requiring callers to construct CLI/`clap` types.
- **SRS-API-002:** The Sanger Variant Analysis API MUST accept explicit AB1, reference FASTA, and configuration paths and MUST reuse the same scientific read-processing, alignment, and variant-calling behavior as the production CLI.
- **SRS-API-003:** The API MUST return typed canonical result data containing input/configuration identities, reference identity, covered reference segments, and normalized reportable variants. The public result MUST NOT expose the versioned JSON report DTO as its domain contract.
- **SRS-API-004:** Calling the Variant Analysis API MUST NOT create operational log files, derive or validate CLI output paths, serialize JSON, or publish result files.
- **SRS-API-005:** Public reference segments MUST use explicit zero-based half-open coordinates and public variant positions MUST use explicit one-based coordinates.
