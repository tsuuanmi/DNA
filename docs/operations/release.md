# Release Operations

This is the operational companion to ADR-0018.

A release record should identify the exact:

- DNA semantic version;
- source revision;
- Rust/Cargo toolchain;
- Cargo.lock identity;
- supported target triple(s);
- artifact SHA-256;
- schema/config versions;
- dependency audit status;
- CI gate status;
- fuzz/adversarial status required for the release;
- synthetic regression status;
- approved real-AB1 validation status;
- documented runtime and peak-memory measurement.

## Release states

Use explicit language:

- **build verified** — engineering build/check gates passed;
- **scientifically validated for the documented corpus/domain** — approved real-trace evidence passed;
- **production-ready release** — ADR-0018 release contract is satisfied.

Do not collapse these into one status.

## Artifact behavior

Existing versioned scientific JSON contracts must not gain unversioned build metadata. Software/build provenance belongs in a separately designed versioned contract or a new result-schema version when required.


## Automated delivery

Tags matching `v*` trigger `.github/workflows/release.yml`. The workflow verifies that the tag matches the crate version, reruns the required Rust engineering gates, builds with the pinned release toolchain and locked dependency graph, records source/toolchain/lockfile provenance, produces SHA-256 checksums, and publishes a Linux x86_64 GitHub Release artifact.

This automation establishes a **build verified** delivery artifact. It does not by itself satisfy the scientific-validation, dependency-audit, fuzz/adversarial, performance, and real-corpus evidence required by ADR-0018 for a **production-ready release**.
