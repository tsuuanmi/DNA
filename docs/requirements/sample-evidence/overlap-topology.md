---
decided-by: [ADR-0023, ADR-0030, ADR-0032, ADR-0033, ADR-0035, ADR-0053, ADR-0068]
verified-by: [crates/dna-sample/src]
---

# Overlap and Topology Requirements

**Requirement namespace:** `SRS-SAMPLE-007..SRS-SAMPLE-012`

These requirements are part of the canonical [sample evidence SRS](README.md) and cover pairwise overlap admission, coverage topology, and normalized-variant support topology.

- **SRS-SAMPLE-007:** After deterministic SHA-256 read ordering, sample reconciliation MUST discover an unordered read-pair overlap only from shared mapped reference coordinates. A pair with no shared coordinate MUST produce no overlap edge, and lack of an overlapping or canonical F/R partner MUST NOT reject an otherwise valid read.
- **SRS-SAMPLE-008:** Pairwise overlap evidence MUST distinguish all shared reference positions from comparable nucleotide positions. Only coordinates where both reads have canonical A/C/G/T query bases that the core trusts for their base (scored, clear of the signal edges, and with no modality reason against the base, as variant calling and the consensus judge them) MAY enter the nucleotide agreement denominator; unresolved symbols and deletions MUST remain outside that denominator, and indel/gap evidence MUST remain in its dedicated locus/variant evidence rather than receive fabricated nucleotide agreement.
- **SRS-SAMPLE-009:** A materialized overlap pair MUST be eligible for later consensus reconciliation exactly when its comparable canonical-base count is at least `minimum_comparable_bases` and its canonical-base agreement fraction is at least `minimum_overlap_agreement`. Failed rules MUST remain explicit as `comparable_bases_below_minimum` and/or `agreement_below_minimum` without mutating either read.
- **SRS-SAMPLE-010:** Compact sample evidence MUST expose pairwise overlap/admission evidence in addition to the existing read registry, sparse differential loci, and normalized variant support. This evidence MUST NOT itself produce a consensus sequence, sample-level variant verdict, majority vote, genotype, heteroplasmy estimate, or reference vote.
- **SRS-SAMPLE-011:** Sample evidence MUST expose a run-length encoded reference-coordinate coverage topology derived only from selected post-trim read segments. Each maximal interval MUST report total read depth plus forward/reverse orientation depth, with total depth equal to the orientation-depth sum. Coverage MUST include every independently placed read regardless of pairwise overlap eligibility and MUST NOT imply canonical-base agreement, consensus admission, or biological strand independence.
- **SRS-SAMPLE-012:** Every normalized sample variant MUST retain a factorized support topology derived from its authoritative per-read support records and selected read orientations: total observing reads, eligible observing reads, forward/reverse observing reads, and eligible forward/reverse observing reads. Each read MUST contribute at most once to one normalized variant identity. The topology MUST include ineligible observations, MUST NOT count covering reference/unresolved/competing-event reads as support for the variant, and MUST NOT be interpreted as confidence, vote weight, sample-level verdict, or biological strand independence.

## Implemented in

<!-- docs:begin backlinks -->
- **Implemented in:** [design/pipeline](../../design/pipeline.md), [design/sample-evidence](../../design/sample-evidence/README.md), [reference/rust-api](../../reference/rust-api.md), [reference/rust-api-sample](../../reference/rust-api-sample.md), [reference/sample-evidence/overlaps-loci](../../reference/sample-evidence/overlaps-loci.md), [reference/sample-evidence/reads-coverage](../../reference/sample-evidence/reads-coverage.md), [requirements/callability](../callability.md), [requirements/consensus](../consensus.md), [requirements/output](../output.md), [requirements/sample-evidence/callability](callability.md), [requirements/sample-evidence/core](core.md), [requirements/signal-processing](../signal-processing.md), [validation/evaluation-harness](../../validation/evaluation-harness.md), [validation/evaluation-measures](../../validation/evaluation-measures.md), [validation/evaluation-triage](../../validation/evaluation-triage.md)
<!-- docs:end backlinks -->
