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
- [CI/CD](../engineering/ci-cd.md)
- [release engineering](../engineering/release.md)
- [validation](../validation/README.md)
- [security](../security/README.md)
- [ADR-0018 rationale](../decisions/adr/0018-production-readiness-release-contract.md)

The checklist is living current-state readiness documentation. ADR-0018 preserves why this release contract exists.
