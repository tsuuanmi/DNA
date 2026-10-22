---
decided-by: [ADR-0065, ADR-0075]
verified-by: [crates/dna-eval/src]
unverified:
  - SRS-COMPAT-001 process rule (reference-validation policy)
  - SRS-COMPAT-002 process rule (approved differential evidence)
---

# System Constraints

**Requirement namespace:** `SRS-COMPAT-*`

These cross-cutting constraints govern external-reference comparison, compatibility interpretation, and where comparison is performed. They are part of the canonical [[requirements/SRS|DNA SRS]].

- **SRS-COMPAT-001:** Apollo comparisons MUST follow [[governance/reference-validation|reference validation policy]]; known defects are intentional divergences, not parity failures or backward-compatibility obligations.
- **SRS-COMPAT-002:** Approved differential evidence MUST compare exact decoded arrays and unaffected deterministic results; normalized variants compare by full tuple without ignoring extras/missing calls.
- **SRS-COMPAT-003:** Logic that compares DNA's results with external call sets, reviewer truth, or other tools (for example Sequencher), and parsers for their formats, MUST live only in the evaluation harness crate `dna-eval` or in downstream pipelines on DNA's published results. The `dna` facade and the plugin crates MUST NOT contain them, and other crates MUST NOT depend on `dna-eval` ([[ADR-0065]], [[ADR-0075]]).

## Implemented in

<!-- docs:begin backlinks -->
- **Implemented in:** [architecture/overview](../architecture/overview.md), [governance/data](../governance/data.md), [validation/evaluation-harness](../validation/evaluation-harness.md), [validation/evaluation-measures](../validation/evaluation-measures.md), [validation/evaluation-workload](../validation/evaluation-workload.md)
<!-- docs:end backlinks -->
