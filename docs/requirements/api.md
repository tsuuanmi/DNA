# Public API Requirements

**Requirement namespace:** `SRS-API-*`

These requirements are part of the canonical [DNA SRS](SRS.md).

- **SRS-API-001:** DNA MUST expose a Rust library capability for reference-guided Sanger Variant Analysis without requiring callers to construct CLI/`clap` types.
- **SRS-API-002:** The Sanger Variant Analysis API MUST accept explicit AB1, reference FASTA, and configuration paths and MUST reuse the same scientific read-processing, alignment, and variant-calling behavior as the production CLI.
- **SRS-API-003:** The Sanger Variant Analysis API MUST return typed result data containing input/configuration identities, reference identity, covered reference segments, and evidence-backed reportable variants. The public result MUST NOT expose the versioned JSON report DTO as its domain contract.
- **SRS-API-004:** Calling the Variant Analysis API MUST NOT create operational log files, derive or validate CLI output paths, serialize JSON, or publish result files.
- **SRS-API-005:** Public reference segments MUST use explicit zero-based half-open coordinates and public variant positions MUST use explicit one-based coordinates.

- **SRS-API-006:** `VariantAnalysisResult` MUST expose a typed `CalledVariantSet` projection that carries one reference identity plus the reportable called variants without implying a right/left alignment, nomenclature, or interchange-format normalization policy.
- **SRS-API-007:** DNA MUST expose optional post-calling variant normalization as a separate Rust capability accepting an explicit normalization policy; using Variant Analysis MUST NOT implicitly invoke that capability.

- **SRS-API-008:** DNA MUST expose a reusable Sanger analysis context that validates and materializes one reference plus scientific configuration once, can analyze multiple independent ABIF traces without rereading those shared inputs, and is safe to share across threads. The reusable context MUST produce the same typed scientific result as the one-shot Sanger Variant Analysis API for the same trace/reference/configuration.
- **SRS-API-009:** The reusable Sanger analysis context MUST NOT own a global thread pool, create operational logs, publish files, or impose batch scheduling policy. Parallel work scheduling belongs to outer orchestration.
