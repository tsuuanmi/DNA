# Sample Evidence Aggregation Method

This is the canonical method entry point for sample aggregation after every trace has independently completed the one-read observation path.

Normative requirements live in [SRS-SAMPLE-*](../srs/sample-evidence.md), public serialization semantics live in the [sample evidence contract](../contracts/sample-evidence.md), and rationale lives in the relevant ADR decision family.

## Methods

- [Read registry](sample-evidence/read-registry.md): validation, duplicate-content rejection, deterministic SHA ordering, and reviewer-facing read identity.
- [Coverage and overlap](sample-evidence/coverage-overlap.md): run-length mapped coverage plus Tracy-derived pairwise overlap/admission.
- [Call signal projection](sample-evidence/call-signal.md): one authoritative reference-oriented projection of per-call signal evidence.
- [Differential loci and profile geometry](sample-evidence/locus-profiles.md): sparse locus selection, support topology, contribution eligibility, profile support/means, heterogeneity, and directional distance.
- [Normalized variants](sample-evidence/variants.md): normalized-variant grouping, eligibility preservation, reviewer evidence, and support topology.

## Interpretation boundary

These methods preserve and summarize evidence for future reconciliation. They do not emit a consensus sequence, adjudicated sample-level variant verdict, genotype, quantitative heteroplasmy, or biological-independence claim.
