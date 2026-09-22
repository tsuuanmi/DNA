# Contributing

DNA is scientific software. Changes are expected to preserve both software
invariants and the documented interpretation boundary.

## Before changing code

Read the repository `README.md`, then use `docs/README.md` to locate the relevant
requirements, architecture/invariants, decisions, design, reference contracts, and
validation policy. Read the nearest source-directory `README.md` for implementation
ownership. Research and proposals are non-normative until promoted through the
documented lifecycle.

## Source boundaries

- Executable production source under `src/` is Rust; source-local `README.md` files are the only non-Rust files allowed there.
- Python under `tools/python/` is for research, validation, orchestration, and
  tests only and must not become a runtime dependency of the `dna` binary.
- Scientific behavior changes require matching tests and documentation.
- Breaking JSON output changes require a new contract/schema version.

## Required verification

```bash
cargo fmt --all --check
cargo shear --deny-warnings
cargo check --locked --all-targets --all-features
cargo clippy --locked --all-targets --all-features -- -D warnings
cargo test --locked --all-targets --all-features
RUSTDOCFLAGS="-D warnings" cargo doc --locked --no-deps --all-features
cargo build --locked --release
```

For Python companion tooling:

```bash
uv sync --project tools/python --frozen
cd tools/python
uv run ruff format --check scripts tests
uv run ruff check scripts tests
uv run basedpyright scripts tests
uv run python -m unittest discover -s tests -p 'test_*.py'
```

CI additionally runs workflow security analysis, dependency policy/review, MSRV checks,
CodeQL, fuzz smoke, release-package smoke testing, schema/reference validation, Python validation-tool tests, and scheduled security checks. The
aggregate `CI success` job is the required branch-protection signal for the main CI workflow.

## Pull requests

Keep changes focused and explain:

- the failure mode or requirement being addressed;
- scientific behavior changes, if any;
- tests/evidence added;
- contract/schema impact;
- operational or migration impact.

Do not suppress warnings or preserve legacy compatibility paths unless an active
requirement explicitly calls for them.
