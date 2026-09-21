# Python Companion Tooling

DNA's production core and runtime are Rust. This directory contains Python only for repository support work:

- research and exploratory analysis;
- validation and contract checks;
- batch/corpus orchestration;
- tests for the Python tooling itself.

Python code here must not become a runtime dependency of the `dna` binary or be placed under the Rust production `src/` tree.

## Environment

```bash
uv sync --project tools/python --frozen
```

When working inside this directory:

```bash
uv run ruff format --check scripts tests
uv run ruff check scripts tests
uv run basedpyright scripts tests
uv run python -m unittest discover -s tests -p 'test_*.py'
uv run python scripts/validate_result_schemas.py
uv run python scripts/validate_rust_source_policy.py
```

The lockfile in this directory belongs only to these Python development tools. Rust dependencies remain authoritative in the root `Cargo.toml` and `Cargo.lock`.
