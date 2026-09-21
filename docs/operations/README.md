# Development and Release Operations

- [Development readiness](development-readiness.md): whether the active scientific core is sufficiently defined for focused implementation.
- [Local batch orchestration](batch.md): local corpus selection, clean reruns, and `scripts/analyze_samples.py` operational behavior.
- [CI lanes](ci.md): required verification gates.
- [Release operations](release.md): release procedure and evidence.
- [Security and trust boundaries](security.md): operational security boundaries.
- [Delivery record](delivery-record.md): historical delivery acceptance record.
- [Production readiness ADR](../adr/0018-production-readiness-release-contract.md): release-quality decision contract.

Data privacy, provenance, and approval policy live under [governance/data.md](../governance/data.md). Repository commands that agents and developers must run are defined by root `AGENTS.md` and CI. Operational documentation must not redefine scientific semantics.
