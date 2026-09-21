# CI and Verification Lanes

CI exists to protect known invariants, not to maximize the number of badges.

The repository root is a Rust project. First-party production source under `src/` is Rust-only. Python is isolated under `tools/python/` and is permitted only for research, validation, test, and repository tooling.

## Pull-request lane

Every pull request runs four independent required jobs.

### Rust quality

The release toolchain is pinned by `rust-toolchain.toml`. The Rust quality job runs:

```bash
cargo fmt --all --check
cargo check --locked --all-targets
cargo clippy --locked --all-targets -- -D warnings
cargo test --locked --all-targets
RUSTDOCFLAGS="-D warnings" cargo doc --locked --no-deps
cargo build --locked --release
```

`--locked` prevents CI from silently changing dependency resolution. Rustdoc warnings are release-blocking alongside compiler and Clippy warnings.

### Minimum supported Rust version

`Cargo.toml` declares the MSRV. CI independently installs that exact compiler and verifies:

```bash
cargo +1.85.0 check --locked --all-targets
```

The MSRV and release toolchain are intentionally separate: the former is a compatibility promise; the latter is the reproducible toolchain used for release-quality checks.

### Dependency audit

CI installs the pinned `cargo-audit 0.22.2` tool and audits the committed `Cargo.lock` against the RustSec advisory database.

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

The Rust source-policy gate complements compiler/Clippy checks by rejecting explicit production compatibility scaffolding that could otherwise be intentionally suppressed: `#[deprecated]` APIs, legacy/backward-compatibility feature gates or declarations, and `allow`/`expect` escape hatches for deprecated/dead/unreachable/unused code. It is deliberately narrow: it does not claim to prove that all conceptual legacy code has been detected.

## Release / delivery lane

Pushing a version tag matching `v*` triggers the release workflow.

The workflow:

1. verifies the tag exactly matches the crate version;
2. reruns formatting, compilation, Clippy, tests, and Rustdoc with the locked dependency graph;
3. builds the release binary;
4. records Rust/Cargo identity, source revision, and `Cargo.lock` checksum;
5. packages the Linux `x86_64-unknown-linux-gnu` artifact;
6. produces SHA-256 checksums;
7. publishes the artifact and checksums to the GitHub Release for that tag.

The current automated binary support claim is therefore Linux x86_64 only. Other platforms are not implied to be release-supported until they are built, tested, and published by the release process.

## Extended lane

Checks with higher runtime or specialized toolchains may run on a schedule, release candidate, or targeted change:

- fuzz campaigns;
- mutation testing;
- property-test expansion;
- dependency/license/advisory audit;
- performance regression measurements;
- approved real-AB1 regression corpus.

A check should be added only when its protected failure mode is documented.

## Failure ownership

- Rust source-policy failure: obsolete/compatibility scaffolding or a diagnostic suppression that must be removed or explicitly redesigned;
- formatter/lint/compiler/Rustdoc failure: engineering defect;
- MSRV failure: declared compatibility or dependency-resolution defect;
- schema/example mismatch: contract defect;
- synthetic test failure: algorithm/implementation regression;
- real-trace disagreement: scientific validation issue requiring analysis, not automatic suppression;
- dependency audit failure: supply-chain/release blocker unless explicitly reviewed;
- release provenance/checksum failure: delivery blocker.

See [release operations](release.md), the [production-readiness ADR](../adr/0018-production-readiness-release-contract.md), and the [data policy](../data.md).
