# CI and Verification Lanes

CI exists to protect known invariants, not to maximize the number of badges.

## Pull-request lane

The fast required lane should cover deterministic checks that are expected on every code change:

```bash
uv run ruff format --check scripts/
uv run ruff check scripts/
uv run basedpyright scripts/
uv run python scripts/validate_result_schemas.py
uv run python scripts/validate_rust_source_policy.py
cargo fmt --all --check
cargo check --all-targets
cargo clippy --all-targets -- -D warnings
cargo test --all-targets
cargo doc --no-deps
```

The Rust source-policy gate complements compiler/Clippy checks by rejecting explicit production compatibility scaffolding that could otherwise be intentionally suppressed: `#[deprecated]` APIs, legacy/backward-compatibility feature gates or declarations, and `allow`/`expect` escape hatches for deprecated/dead/unreachable/unused code. It is deliberately narrow: it does not claim to prove that all conceptual legacy code has been detected.

Repository-specific reference/config/docs-mirror checks remain required when present.

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
- formatter/lint/compiler failure: engineering defect;
- schema/example mismatch: contract defect;
- synthetic test failure: algorithm/implementation regression;
- real-trace disagreement: scientific validation issue requiring analysis, not automatic suppression;
- dependency audit failure: supply-chain/release blocker unless explicitly reviewed.

See [release operations](release.md) and the [data policy](../data.md).
