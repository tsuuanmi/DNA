# Testing

Testing protects encoded properties; it does not by itself establish scientific correctness.

Use the test type that matches the risk:

- unit tests for focused deterministic behavior;
- integration tests for command and cross-module behavior;
- regression tests for fixed defects;
- schema/contract tests for public shapes;
- property/fuzz/adversarial tests for invariant and untrusted-input classes;
- approved real-data validation for scientific claims.

Normative minimums are in [validation acceptance criteria](../validation/acceptance-criteria.md). Broader evidence strategy is in [validation](../validation/README.md).
