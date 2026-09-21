# Development and Release Operations

- [Local batch orchestration](batch.md): local corpus selection, clean reruns, and `scripts/analyze_samples.py` operational behavior.
- [CI lanes](ci.md): required verification gates.
- [Release operations](release.md): release procedure and evidence.
- [Security and trust boundaries](security.md): operational security boundaries.
- [Delivery record](delivery-record.md): historical implementation/delivery acceptance record.
- [Production readiness ADR](../adr/0018-production-readiness-release-contract.md): release-quality decision contract.

Current and future product priorities live in the [roadmap](../roadmap.md). Scientific stage contracts belong to the [SRS](../srs/README.md) and [methods](../methods/README.md), not to a separate readiness checklist.

Data privacy, provenance, and approval policy live under [governance/data.md](../governance/data.md). Repository commands that agents and developers must run are defined by root `AGENTS.md` and CI. Operational documentation must not redefine scientific semantics.
