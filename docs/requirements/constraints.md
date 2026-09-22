# System Constraints

**Requirement namespace:** `SRS-COMPAT-*`

These cross-cutting constraints govern external-reference comparison and compatibility interpretation. They are part of the canonical [DNA SRS](SRS.md).

- **SRS-COMPAT-001:** Apollo comparisons MUST follow [reference validation policy](../governance/reference-validation.md); known defects are intentional divergences, not parity failures or backward-compatibility obligations.
- **SRS-COMPAT-002:** Approved differential evidence MUST compare exact decoded arrays and unaffected deterministic results; normalized variants compare by full tuple without ignoring extras/missing calls.
