---
decided-by: [ADR-0065, ADR-0075]
---

# Validation Strategy

DNA selects validation by the failure mode being protected.

## Engineering confidence

- formatting, compilation, lint/static analysis, and rustdoc;
- unit and integration tests for deterministic behavior;
- schema/example validation for public contracts;
- focused regression tests for fixed defects;
- property/fuzz/adversarial testing where malformed input or invariants justify it.

## Scientific confidence

- synthetic fixtures establish deterministic boundary behavior without sensitive data;
- approved real AB1 evidence is required for scientific release confidence;
- external implementations are differential evidence, not automatic ground truth;
- disagreements are analyzed rather than hidden by compatibility fallbacks;
- comparison with reviewer truth or other tools runs in the [[validation/evaluation-harness|evaluation harness]] (ADR-0075) or in downstream pipelines on DNA's published results (ADR-0065); its outputs are evidence about a DNA revision, not DNA contracts;
- local real traces are exploratory evidence: threshold choices are recorded in decisions as aggregates measured against a named revision, and release claims still require the approval record of [[governance/data|data governance]].

## Release confidence

A production-ready release must satisfy [[validation/acceptance-criteria|acceptance criteria]], required CI gates, approved evidence, dependency/security review, and the [[operations/production-readiness|production-readiness record]].

Coverage alone is not treated as proof of domain correctness.
