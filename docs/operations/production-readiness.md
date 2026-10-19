---
decided-by: [ADR-0018, ADR-0021]
---

# Production Readiness

Production readiness is an evidence-backed release state, not a label inferred from a successful build.

## Required evidence

- [ ] normative requirements are current;
- [ ] architecture/design/reference documentation matches intended behavior;
- [ ] required CI/CD gates pass;
- [ ] validation acceptance criteria pass;
- [ ] approved real-trace evidence supports the documented scientific scope;
- [ ] performance evidence is documented for the release environment;
- [ ] dependency/advisory review is complete;
- [ ] security/trust-boundary review is current;
- [ ] operational logs/failure handling are documented;
- [ ] release artifact identity and provenance are recorded.

Evidence locations:

- [requirements](../requirements/README.md)
- [[engineering/ci-cd|CI/CD]]
- [[engineering/release|release engineering]]
- [[engineering/platform-support|platform support]]
- [[operations/release-evidence-template|release evidence template]]
- [[governance/repository|repository governance]]
- [validation](../validation/README.md)
- [security](../security/README.md)
- [[ADR-0018|ADR-0018 rationale]]

The checklist is living current-state readiness documentation. ADR-0018 preserves why this release contract exists.
