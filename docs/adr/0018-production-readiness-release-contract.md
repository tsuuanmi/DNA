# ADR-0018: Production readiness is an explicit release contract

## Status

Accepted

## Context

DNA already has a strong deterministic MVP baseline: checked ABIF parsing, typed failures, strict configuration, warnings-denied Clippy, focused unit/integration tests, versioned JSON schemas, atomic publication, bounded alignment, and explicit biological limitations.

Those properties are necessary but not sufficient for a production-ready scientific tool. A production claim also requires a stable toolchain policy, dependency governance, adversarial parser validation, release provenance, explicit platform support, and repeatable release evidence.

Treating "production-ready" as a feature milestone would make the claim ambiguous. The project instead needs a release contract that can be checked independently of any particular feature roadmap.

## Decision

DNA will treat production readiness as an evidence-backed release profile, not as a synonym for "the MVP is implemented".

A release may be described as production-ready only when all mandatory gates below are satisfied and recorded for the exact release revision.

### 1. Toolchain and compiler policy

- `Cargo.toml` MUST declare the supported `rust-version`.
- CI MUST verify the declared minimum Rust version and the project-selected stable release toolchain.
- The release toolchain MUST be explicit and reproducible from repository metadata; silently following an unpinned moving toolchain is insufficient for a release build.
- Formatting and lint behavior used as release gates MUST therefore be tied to the selected toolchain.
- Supported target platforms MUST be documented. DNA MUST NOT imply support for a platform that is not built and tested by the release process.

### 2. Source-quality gates

The normal pull-request lane MUST remain warnings-clean and include, at minimum:

```text
cargo fmt --all --check
cargo check --all-targets
cargo clippy --all-targets -- -D warnings
cargo test --all-targets
cargo doc --no-deps
```

Repository-specific Python/schema/reference checks remain part of the same quality contract.

Production Rust continues to:

- forbid `unsafe`;
- avoid `unwrap`/`expect` in production paths;
- use typed errors for externally triggered failures;
- keep scientific functions separate from filesystem and logging side effects.

### 3. Adversarial and property validation

ABIF is untrusted binary input. In addition to example-based unit tests, DNA MUST maintain adversarial testing for parsing and coordinate-sensitive logic.

The production validation program will include:

- fuzz targets for ABIF directory/tag decoding and other bounded binary entry points;
- regression fixtures for every minimized crashing or invariant-breaking fuzz input;
- property/invariant tests for coordinate transforms, circular projection, reverse-strand call mapping, normalization idempotence, and serialization invariants.

Fuzzing does not need to run exhaustively on every pull request. A bounded smoke target may run in CI, while longer campaigns can run on a schedule or before releases.

### 4. Dependency and supply-chain policy

Every release dependency graph MUST be checked for:

- known RustSec advisories;
- disallowed or unexpected dependency sources;
- license policy violations;
- yanked or explicitly banned dependencies according to project policy.

Exceptions MUST be documented with a reason and review date. A vulnerability or policy exception must never be silently ignored.

### 5. Release identity and provenance

A production binary MUST have a stable, inspectable identity that can be connected to the source revision and dependency graph used to build it.

At minimum the release record MUST retain:

- DNA semantic version;
- source revision;
- Rust/Cargo toolchain identity;
- `Cargo.lock` identity;
- release artifact SHA-256;
- supported target triple.

Existing JSON contracts MUST NOT be mutated retroactively. If software/build provenance is added to scientific output, it requires a new versioned output contract or another explicitly versioned provenance record.

### 6. Determinism and failure behavior

For identical input bytes, reference bytes, scientific configuration, algorithm versions, and supported execution environment, scientific output MUST remain deterministic.

Externally supplied malformed input MUST result in a bounded typed failure, not a panic, uncontrolled allocation, partial scientific result, or silent fallback.

Atomic no-overwrite result publication remains part of the production contract.

### 7. Release evidence

A release checklist MUST identify the exact revision and record the outcome of:

- source-quality gates;
- schema/configuration/reference validation;
- dependency audit;
- adversarial/fuzz validation status;
- synthetic regression suite;
- approved real-AB1 regression suite;
- documented runtime and peak-memory measurement;
- artifact checksum and toolchain identity.

Passing source checks without approved real-trace scientific evidence is not sufficient for a production scientific release.

## Consequences

### Positive

- "Production-ready" becomes a falsifiable claim with an auditable checklist.
- Toolchain drift and dependency risk become explicit instead of implicit.
- Binary parsing receives validation proportional to its attack/corruption surface.
- Release artifacts can be connected to source, dependencies, and validation evidence.
- Scientific validation remains distinct from generic software correctness.

### Cost

- CI and release maintenance become more involved.
- Fuzzing and real-corpus validation require retained corpora and periodic execution.
- Adding build provenance may require a future output-schema version.
- Some checks will run outside the fast pull-request lane.

## Non-goals

This ADR does not require:

- claiming bit-for-bit reproducible builds across all hosts;
- supporting every Rust target or operating system;
- running long fuzz campaigns on every commit;
- adding a network service, daemon, installer, or package registry publication;
- changing current scientific algorithms.

## References

- Rust Cargo Book: Continuous Integration and `rust-version`.
- RustSec Advisory Database and `cargo-audit`/`cargo-deny`.
- rust-fuzz `cargo-fuzz` documentation.
