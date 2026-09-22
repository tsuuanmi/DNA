# Supply-Chain and Vulnerability Management

Dependency and release integrity are part of the production trust boundary.

## Dependencies

- use locked Rust/Python dependency graphs;
- add dependencies only for a concrete required capability;
- review source, maintenance, license, advisory, and transitive cost;
- do not suppress advisory findings merely to keep a release green.

## Release provenance

A production release records source revision, toolchain, lockfile identity, artifact SHA-256, schema/config versions, CI status, validation status, and dependency-review status.

## Vulnerabilities

A newly discovered vulnerability is triaged by:

1. affected dependency/component and reachable code path;
2. impact on confidentiality, integrity, availability, or scientific correctness;
3. fixed-version or mitigation availability;
4. required source/config/contract changes;
5. regression/security validation before release.

If a mitigation changes architecture, public behavior, or trust boundaries, update the corresponding decision and canonical docs rather than maintaining a hidden workaround.

See [dependency policy](../engineering/dependencies.md) and [release engineering](../engineering/release.md).
