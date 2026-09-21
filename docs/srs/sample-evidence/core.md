# Core Sample Evidence Requirements

**Requirement namespace:** `SRS-SAMPLE-001..SRS-SAMPLE-006`

These requirements are part of the canonical [sample evidence SRS](../sample-evidence.md) and cover independent read processing, identity, sparse evidence, and preserved observations.

- **SRS-SAMPLE-001:** Multi-read sample analysis MUST process every trace independently through the same authoritative read-observation path used by single-read analysis before cross-read reconciliation.
- **SRS-SAMPLE-002:** Cross-read overlap and support MUST be discovered from mapped reference coordinates/variants rather than inferred from canonical F/R pair names.
- **SRS-SAMPLE-003:** A usable read MUST NOT be rejected solely because a nominal F/R partner is missing.
- **SRS-SAMPLE-004:** Sample evidence MUST contain one deterministic read registry sorted by input SHA-256. Reviewer-facing filename-stem provenance, stable SHA-256 identity, trace-integrity evidence, derived orientation, and mapped post-trim coverage MUST be defined once per read; locus and variant evidence MUST reference that registry by unique reviewer-facing read name. Filename, amplicon, primer, or replicate grouping MUST NOT be treated as a scientific placement or exclusive merge key. Duplicate trace content MUST NOT contribute twice.
- **SRS-SAMPLE-005:** Sample evidence MUST serialize locus observations only at positions where at least one covering read is alternate, unresolved, or deleted. Every covering read at such a differential locus MUST remain explicit, including canonical reference support with its observed base and quality. For a read, a position inside its mapped reference segments but absent from `locus_differences[]` is a canonical reference match; a position outside its mapped segments is uncovered. Missing coverage MUST NOT become reference support.
- **SRS-SAMPLE-006:** Normalized variant support MUST retain the contributing read name, configured eligibility, exact exclusion reasons, and reviewer-facing reference-oriented base/peak/quality evidence; read-level filtering MUST NOT erase an observed normalized variant. Future sample variants and discordance states MUST derive from this evidence; a consensus sequence MAY be emitted only as a downstream projection.
