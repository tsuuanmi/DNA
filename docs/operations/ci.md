# CI and Verification Lanes

CI exists to protect documented invariants and security boundaries, not to maximize the number of badges.

The repository root is a Rust project. First-party production source under `src/` is Rust-only. Python is isolated under `tools/python/` and is permitted only for research, validation, test, and repository tooling.

All third-party GitHub Actions are pinned to immutable full commit SHAs. Dependabot owns routine updates to those pins.

## Pull-request lane

### GitHub Actions policy

Workflow changes are checked three ways:

- the repository validator rejects mutable third-party action pins and privileged untrusted triggers;
- `actionlint` validates GitHub Actions syntax and expressions;
- `zizmor` performs GitHub Actions security analysis.

### Rust quality


The release toolchain is pinned by `rust-toolchain.toml`. Every pull request runs:

```bash
cargo fmt --all --check
cargo shear --deny-warnings
cargo check --locked --all-targets --all-features
cargo clippy --locked --all-targets --all-features -- -D warnings
cargo test --locked --all-targets --all-features
RUSTDOCFLAGS="-D warnings" cargo doc --locked --no-deps --all-features
cargo build --locked --release
```

`cargo shear --deny-warnings` rejects unused/misplaced dependencies and unlinked Rust source files. `--locked` prevents CI from silently changing dependency resolution. Rustdoc warnings are release-blocking alongside compiler and Clippy warnings.

### Minimum supported Rust version

`Cargo.toml` declares the MSRV. CI independently installs that exact compiler and verifies:

```bash
cargo +1.88.0 check --locked --all-targets
```

The MSRV and release toolchain are intentionally separate: the former is a compatibility promise; the latter is the reproducible toolchain used for release-quality checks.

### Dependency policy

Dependency verification has three layers:

1. `cargo-deny` checks advisories, yanked crates, licenses, trusted sources, wildcard requirements, banned/replacement crates, and duplicate-version policy.
2. pinned `cargo-audit 0.22.2` independently checks the committed `Cargo.lock` against RustSec.
3. GitHub dependency review rejects pull requests that introduce dependencies with moderate-or-higher known vulnerabilities.

`deny.toml` is the authoritative source/license/bans policy. Exceptions must include a concrete reason and review date rather than silently weakening the global policy.

### Static security analysis

CodeQL analyzes Rust on pull requests, `main`, and a weekly schedule with the `security-extended` query suite. Results are published to GitHub code scanning.

### Adversarial parser validation

The `ABIF fuzz smoke` job exercises the bounds-checked ABIF directory parser with `cargo-fuzz`:

- 30-second campaigns on pull requests and `main`;
- longer scheduled campaigns;
- pinned nightly toolchain and cargo-fuzz version;
- retained minimized regressions when a defect is found.

Fuzzing complements deterministic malformed-input tests; it does not replace them.

### Repository policy and Python companion tooling

Python dependencies are locked under `tools/python/`; they are not runtime dependencies of the DNA binary.

The repository-policy job:

- rejects any non-`.rs` file under `src/`;
- validates the explicit Rust source policy;
- runs Ruff formatting and lint checks;
- runs basedpyright;
- runs Python tooling tests;
- validates result schemas and examples;
- parses the strict TOML configuration template;
- verifies the documented environment template;
- verifies the rCRS checksum and reference length.

The Rust source-policy gate complements compiler/Clippy checks by rejecting explicit production compatibility scaffolding that could otherwise be intentionally suppressed.

### Aggregate required check

The `CI success` job waits for every mandatory job in `.github/workflows/ci.yml` and fails unless all applicable gates succeeded. Branch rules should require this aggregate check instead of duplicating every internal job name, reducing protection drift as CI evolves.

## Scheduled security posture

OpenSSF Scorecard runs on `main` and weekly. Its SARIF output is retained briefly as an Actions artifact and uploaded to GitHub code scanning.

Dependabot monitors:

- Cargo dependencies;
- the isolated uv Python tooling project;
- full-SHA GitHub Actions pins;
- the pinned Rust release toolchain.

## Release / delivery lane

Pushing a version tag matching `v*` triggers the release workflow only after the tagged commit is verified to be reachable from `main`.

The workflow:

1. verifies the tag exactly matches the crate version and belongs to `main`;
2. reruns formatting, compilation, Clippy, tests, and Rustdoc with the locked dependency graph;
3. reruns `cargo-deny` and RustSec audit;
4. builds the release binary with pinned `cargo-auditable` and audits the produced binary;
5. records Rust/Cargo identity, source revision, and `Cargo.lock` checksum;
6. generates an SPDX JSON SBOM from the auditable binary with a pinned Syft version;
7. packages the Linux `x86_64-unknown-linux-gnu` artifact;
8. produces SHA-256 checksums;
9. creates GitHub/Sigstore build-provenance and SBOM attestations;
10. publishes the archive, SBOM, and checksums to the GitHub Release.

Release builds deliberately do not restore shared CI caches; artifact-producing workflows build from source and the locked dependency graph to avoid cache-poisoning risk.

The current automated binary support claim is therefore Linux x86_64 only. Other platforms are not implied to be release-supported until they are built, tested, and published by the release process.

## Extended scientific lane

Checks that cannot be reduced to normal public CI remain release evidence:

- property/invariant expansion;
- long fuzz campaigns;
- performance and peak-memory measurements;
- approved real-AB1 regression corpus;
- ground-truth biological comparison and disagreement analysis.

A check is added only when its protected failure mode is documented.

## Enforcement

CI is not itself an enforcement mechanism. The protected-`main` and protected-`v*` rules described in [repository governance](repository-governance.md) make the required checks non-bypassable.

## Failure ownership

- workflow-policy/actionlint/zizmor failure: CI supply-chain or workflow-security defect;
- source-policy/formatter/lint/compiler/Rustdoc failure: engineering defect;
- MSRV failure: declared compatibility or dependency-resolution defect;
- cargo-deny/audit/dependency-review failure: supply-chain or licensing blocker;
- CodeQL failure/alert: security review blocker until triaged;
- fuzz failure: parser/adversarial correctness blocker;
- schema/example mismatch: contract defect;
- synthetic test failure: algorithm/implementation regression;
- real-trace disagreement: scientific validation issue requiring analysis, not automatic suppression;
- provenance/SBOM/attestation failure: delivery blocker.

See [release operations](release.md), [repository governance](repository-governance.md), the [production-readiness ADR](../adr/0018-production-readiness-release-contract.md), and the [data policy](../data.md).
