# System Constraints

**Requirement namespace:** `SRS-COMPAT-*`

These cross-cutting constraints govern external-reference comparison, compatibility interpretation, and where comparison is performed. They are part of the canonical [DNA SRS](SRS.md).

- **SRS-COMPAT-001:** Apollo comparisons MUST follow [reference validation policy](../governance/reference-validation.md); known defects are intentional divergences, not parity failures or backward-compatibility obligations.
- **SRS-COMPAT-002:** Approved differential evidence MUST compare exact decoded arrays and unaffected deterministic results; normalized variants compare by full tuple without ignoring extras/missing calls.
- **SRS-COMPAT-003:** DNA MUST NOT contain logic that compares its results with external call sets, reviewer truth, or other tools (for example Sequencher), nor parsers for their formats. Comparison runs in downstream pipelines on DNA's published results ([ADR-0065](../decisions/adr/0065-result-comparison-downstream.md)).
