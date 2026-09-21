# Local Batch Orchestration

`scripts/analyze_samples.py` is an external development/validation wrapper around the core DNA CLI. It selects local traces, performs a guarded clean rerun, preserves per-trace results, and publishes one sample aggregate only after every selected trace for that sample succeeds.

This workflow does not change the core CLI contract: `analyze` still accepts exactly one AB1 trace, while `sample` accepts an explicit set of traces for one sample.

## Default local layout

The wrapper defaults to a local ignored corpus shaped like:

```text
data/
├── MS_010426_001.txt
└── raw/
    └── MS_010426_001/
        └── *.ab1
```

The exact corpus, identifiers, and file count are local state and are not repository contracts. Data handling is governed by [data and fixture governance](../governance/data.md).

## Default run

```bash
uv run python scripts/analyze_samples.py
```

By default the wrapper reads the configured manifest prefix, discovers matching traces, uses the bundled reference/configuration, and groups outputs by sample:

```text
results/<sample-id>/
├── <trace-stem>.json
├── ...
└── <sample-id>.json
```

The persistent operational log is `logs/<sample-id>.log`. Helper per-trace analyses run in temporary workspaces so their temporary logs are not promoted as authoritative sample logs.

Use `uv run python scripts/analyze_samples.py --help` to select a different manifest, trace directory, reference, configuration, output/log root, binary, or sample limit.

## Preflight and cleanup

Before deleting any selected artifact, the wrapper:

1. reads and validates the selected manifest prefix;
2. discovers every selected trace;
3. rejects missing matches, ambiguous ownership, duplicate selected IDs, trace-stem collisions, unsafe target types, path escapes, and symlinks;
4. validates that output/log roots do not overlap protected inputs or each other;
5. validates all selected cleanup targets;
6. builds the release binary unless `--no-build` is selected, then requires a regular binary.

Only after successful preflight/build does cleanup remove selected `results/<sample-id>/` directories and selected `logs/<sample-id>.log` files. Unselected results and unrelated/legacy logs are preserved.

## Execution and publication

Each selected trace is run through an isolated one-file `signal analyze` invocation. Its generated JSON is then placed into the sample result directory without overwrite.

After every trace for one sample succeeds, the wrapper runs:

```text
signal sample <sample-id> <trace.ab1>... --reference <reference.fasta>
```

over the complete selected trace set and atomically publishes the aggregate as `results/<sample-id>/<sample-id>.json`.

The wrapper must not weaken the core CLI's atomic no-overwrite behavior.

## Failure semantics

The clean rerun is intentionally not one transaction across the whole selected workload.

- If preflight or build fails, prior selected artifacts are preserved.
- If a trace fails after cleanup, that sample aggregate is not published.
- Earlier successful new results may remain when a later trace/sample fails.
- Removed prior selected artifacts are not restored automatically.
- Final status reports trace and sample failures explicitly.

These semantics are operational behavior of the wrapper, not scientific interpretation.
