# Python Companion Tooling

DNA's production core and runtime are Rust. This directory contains Python only for repository support work:

- repository checks and CI validators (source policy, module layering,
  workflows, documentation structure and the documentation knowledge graph,
  result schemas, the shipped configuration against its reference table);
- measurement and validation helpers that produce evidence about a DNA build;
- research and exploratory analysis that no production path depends on;
- tests for the Python tooling itself.

Orchestrating DNA runs, converting results, and comparing them with other
sources belong to downstream pipelines or to the Rust evaluation harness
`dna-eval`, not here (ADR-0065, ADR-0066, ADR-0075).

Research code lives under `research/<topic>/`, one self-contained `uv` script (PEP 723)
per step, with its own README. It is outside the CI checks and the lockfile below
(for example [the trace simulator](research/trace_simulator/README.md)).

Python code here must not become a runtime dependency of the `dna` binary or be placed under the Rust production `src/` or `crates/*/src/` trees.

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
uv run python scripts/validate_docs_structure.py
uv run python scripts/docs_graph.py check
uv run python scripts/docs_graph.py sync --check
uv run python scripts/validate_module_layers.py
uv run python scripts/validate_workflow_policy.py
uv run python scripts/validate_configuration_reference.py
```

The lockfile in this directory belongs only to these Python development tools. Rust dependencies remain authoritative in the workspace `Cargo.toml` files and the root `Cargo.lock`.


## Documentation graph

`scripts/docs_graph.py` implements the checks, generators, and commands of the
documentation standard (`docs/governance/documentation-architecture.md`) for the
profile in `docs/governance/documentation.md`: `check`, `sync [--check]`,
`trace <path|ID>`, `impact <base>`, `candidates`, and `retire`. Run `sync` after
changing front matter edges or statuses; never edit a generated block by hand.

## Performance evidence

For Linux release-validation runs, measure one exact command and write machine-readable
runtime/CPU/peak-RSS evidence:

```bash
uv run python scripts/measure_command.py \
  --output ../../results/validation/performance.json \
  -- ../../target/release/dna analyze TRACE.ab1 --reference ../../references/rCRS.fasta
```

The wrapper executes the command directly without a shell and propagates its exit status.
