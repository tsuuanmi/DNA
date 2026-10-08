# ADR-0066: Python in DNA is limited to repository tooling

- **Status:** Accepted
- **Date:** 2026-10-08
- **Related decisions:** [ADR-0021](0021-scientific-core-confidence-floor.md),
  [ADR-0065](0065-result-comparison-downstream.md)

## Context

DNA's production runtime is Rust. Python under `tools/python/` grew two kinds of
code:

- repository tooling: CI validators for source policy, module layering,
  workflows, documentation structure, and result schemas; a command
  measurement helper for release evidence; and their tests;
- `analyze_samples.py`, a batch runner with its own requirements (`SRS-BAT-*`)
  and runbook. It selected manifest samples, deleted earlier results, ran the
  CLI per trace and per sample, and arranged per-sample result folders.

The batch runner is application orchestration, not a check of the repository.
Selecting samples, running DNA over a corpus, and organising its outputs is the
job of downstream pipelines such as `mtdna_raw`, which also convert and compare
results (ADR-0065).

## Decision

Python in DNA is limited to:

- repository checks and CI validators;
- tests of that tooling;
- measurement and validation helpers that produce evidence about a DNA build;
- research and exploratory analysis that no production path depends on.

Orchestration of DNA runs, conversion of results, and comparison belong to
downstream pipelines. The batch runner, `SRS-BAT-*`, and the batch runbook are
removed. This revises ADR-0021's list of capabilities to keep: batch
orchestration leaves DNA for a boundary reason, not to simplify the MVP.

## Consequences

- DNA's operational surface is the `analyze`, `basecall`, and `sample` CLI and
  the Rust library; downstream pipelines drive them.
- Per-sample result folders (`results/<sample-id>/`) are no longer produced by
  anything in DNA.
- No Python file may become a runtime dependency of `dna` or a required step of
  producing DNA results.
