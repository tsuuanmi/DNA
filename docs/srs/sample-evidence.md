# Sample Evidence Requirements

**Requirement namespace:** `SRS-SAMPLE-*`

These requirements are part of the canonical [DNA SRS](README.md). The family is split by durable sample-evidence responsibility so contributors can read the relevant subset without duplicating requirements.

| Range | Scope | Document |
|---|---|---|
| `SRS-SAMPLE-001..006` | independent reads, identity, sparse observations, preserved variant support | [Core](sample-evidence/core.md) |
| `SRS-SAMPLE-007..012` | pairwise overlap admission, coverage, support topology | [Overlap and topology](sample-evidence/overlap-topology.md) |
| `SRS-SAMPLE-013..018` | locus profiles, noisy context, signal evidence, contribution eligibility | [Signal evidence](sample-evidence/signal-evidence.md) |
| `SRS-SAMPLE-019..023` | support accumulation, mean profiles, heterogeneity, directional distance | [Profile aggregation](sample-evidence/profile-aggregation.md) |

Requirement IDs are stable. Cross-family documents should reference these IDs rather than restating their normative language.
