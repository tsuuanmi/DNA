---
implements: [SRS-SAMPLE]
decided-by: [ADR-0023, ADR-0030, ADR-0053, ADR-0068]
verified-by: [crates/dna-sample/src]
---

# Core Sample Evidence Requirements

**Requirement namespace:** `SRS-SAMPLE-001..SRS-SAMPLE-006`

These requirements are part of the canonical [sample evidence SRS](README.md) and cover independent read processing, identity, sparse evidence, and preserved observations.

- **SRS-SAMPLE-001:** Multi-read sample analysis MUST process every trace independently through the same authoritative read-observation path used by single-read analysis before cross-read reconciliation. Traces MAY be processed concurrently; the result MUST NOT depend on processing order or scheduling, and a read failure that is not a rejection (SRS-SAMPLE-027, SRS-SAMPLE-029, SRS-SAMPLE-030) MUST be reported for the first such trace in input order.
- **SRS-SAMPLE-002:** Cross-read overlap and support MUST be discovered from mapped reference coordinates/variants rather than inferred from canonical F/R pair names.
- **SRS-SAMPLE-003:** A usable read MUST NOT be rejected solely because a nominal F/R partner is missing.
- **SRS-SAMPLE-004:** Sample evidence MUST contain one deterministic read registry sorted by input SHA-256. Reviewer-facing filename-stem provenance, stable SHA-256 identity, trace-integrity evidence, derived orientation, and mapped post-trim coverage MUST be defined once per read; locus and variant evidence MUST reference that registry by unique reviewer-facing read name. Filename, amplicon, primer, or replicate grouping MUST NOT be treated as a scientific placement or exclusive merge key. Duplicate trace content MUST NOT contribute twice.
- **SRS-SAMPLE-005:** Sample evidence MUST serialize locus observations only at positions where at least one covering read is alternate, unresolved, or deleted by an unmasked call (SRS-SAMPLE-024). Every covering read at such a differential locus MUST remain explicit, including canonical reference support with its observed base and quality and masked observations. For a read, a position inside its callable reference segments but absent from `locus_differences[]` is a canonical reference match; a position inside its mapped segments but outside its callable segments is covered by a masked or an unresolved call and is neither reference nor difference; a position outside its mapped segments is uncovered. Missing or masked coverage MUST NOT become reference support.
- **SRS-SAMPLE-006:** Normalized variant support MUST retain the contributing read name, configured eligibility, exact exclusion reasons, and reviewer-facing reference-oriented base/peak/quality evidence; read-level filtering MUST NOT erase an observed normalized variant. Future sample variants and discordance states MUST derive from this evidence; a consensus sequence MAY be emitted only as a downstream projection.

## Implemented in

<!-- docs:begin backlinks -->
- **Implemented in:** [design/pipeline](../../design/pipeline.md), [design/sample-evidence](../../design/sample-evidence/README.md), [reference/rust-api](../../reference/rust-api.md), [reference/rust-api-sample](../../reference/rust-api-sample.md), [reference/sample-evidence/overlaps-loci](../../reference/sample-evidence/overlaps-loci.md), [reference/sample-evidence/reads-coverage](../../reference/sample-evidence/reads-coverage.md), [requirements/api](../api.md), [requirements/callability](../callability.md), [requirements/consensus](../consensus.md), [requirements/output](../output.md), [requirements/sample-evidence/callability](callability.md), [requirements/signal-processing](../signal-processing.md), [validation/evaluation-harness](../../validation/evaluation-harness.md), [validation/evaluation-measures](../../validation/evaluation-measures.md), [validation/evaluation-triage](../../validation/evaluation-triage.md)
<!-- docs:end backlinks -->
