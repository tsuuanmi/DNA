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
- disagreements are analyzed rather than hidden by compatibility fallbacks.

## Release confidence

A production-ready release must satisfy [acceptance criteria](acceptance-criteria.md), required CI gates, approved evidence, dependency/security review, and the [production-readiness record](../operations/production-readiness.md).

Coverage alone is not treated as proof of domain correctness.
