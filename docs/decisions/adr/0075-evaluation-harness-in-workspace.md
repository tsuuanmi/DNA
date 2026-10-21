---
id: ADR-0075
status: accepted
date: 2026-10-14
family: [crates, comparison]
supersedes: [ADR-0065]
related: [ADR-0018, ADR-0065, ADR-0066, ADR-0070, PROP-0005]
---

# ADR-0075: Evaluation harness as a workspace crate

<!-- docs:begin lifecycle -->
- **Status:** Accepted
- **Date:** 2026-10-14
- **Supersedes:** [ADR-0065](0065-result-comparison-downstream.md)
- **Related:** [ADR-0018](0018-production-readiness-release-contract.md), [ADR-0066](0066-python-limited-to-repository-tooling.md), [ADR-0070](0070-workspace-split-by-plugin-family.md), [PROP-0005](../../proposals/0005-evaluation-feedback-loop.md)
<!-- docs:end lifecycle -->

<!-- docs:begin backlinks -->
- **Governs:** [architecture/context](../../architecture/context.md), [architecture/interfaces](../../architecture/interfaces.md), [architecture/overview](../../architecture/overview.md), [engineering/ci-cd](../../engineering/ci-cd.md), [engineering/dependencies](../../engineering/dependencies.md), [engineering/development](../../engineering/development.md), [governance/data](../../governance/data.md), [requirements/constraints](../../requirements/constraints.md), [validation/acceptance-criteria](../../validation/acceptance-criteria.md), [validation/artefact-bench](../../validation/artefact-bench.md), [validation/evaluation-harness](../../validation/evaluation-harness.md), [validation/evaluation-measures](../../validation/evaluation-measures.md), [validation/evaluation-triage](../../validation/evaluation-triage.md), [validation/strategy](../../validation/strategy.md)
<!-- docs:end backlinks -->

## Context

[[ADR-0065]] keeps every comparison with
reviewer truth out of DNA. Downstream pipelines read DNA's published JSON,
convert it, and compare it there.

[[PROP-0005]] needs a loop
that scores a change on a fixed test panel in seconds to minutes. Through the
CLI, every sample costs two processes and two JSON documents, and a rule change
repeats the whole run. A Rust harness that calls DNA in process removes that
cost. The user chose to keep that harness in the DNA workspace (PROP-0005,
option B), so that a change to DNA and its evaluation can land together.

## Decision

1. The evaluation harness is the workspace crate `dna-eval` in `crates/`. It
   builds the `dna-eval` binary, which is not part of the `dna` binary or the
   release package.
2. Truth adapters (parsers of reviewer truth and panel manifests) and
   comparison with reviewer truth or other tools may live in `dna-eval` and
   nowhere else in the repository.
3. `dna-eval` may depend on the `dna` facade and on the plugin crates. No crate
   may depend on `dna-eval`. The `dna` facade and the plugin crates
   (`dna-kernel`, `dna-core`, `dna-sanger`, `dna-post`) stay free of truth
   formats and comparison. `validate_module_layers.py` enforces this graph.
4. The harness runs DNA through the public library, the same path as the
   CLI. It does not reach into DNA internals to compare. Its comparison
   adapts the lab's rules (`mtdna_raw`), so that a number means what the lab
   means by it.
5. Harness outputs are validation evidence about a DNA revision. They are not
   DNA result contracts and are kept outside the repository under the
   [[governance/data|data policy]]. Tests of the harness use synthetic
   fixtures only.

## Alternatives

- **A separate Rust repository** (PROP-0005, option A). It keeps ADR-0065
  whole, but DNA and the harness would be versioned apart and local work
  needs a path dependency on a DNA checkout.
- **Downstream Python driving the CLI** (option C). Nothing new to build, but
  a process per sample and JSON on every stage cannot meet the loop budgets.
- **Comparison inside the `dna` binary.** Rejected: lab truth formats would
  become part of the shipped product and its contracts.

## Consequences

- One change can update DNA and its evaluation together, and the harness
  builds and tests with the workspace.
- DNA exposes the `consensus` then `call` workflow as a library call
  ([[reference/rust-api-sample#consensus-call|`SampleAnalysis::call`]]), so the
  harness needs no CLI process or JSON parse per sample.
- The release still builds and ships only `dna` (`cargo build -p dna`).
- Downstream pipelines may keep comparing published results; ADR-0065 still
  holds for them.
- [[ADR-0066]] is unchanged:
  Python in DNA still runs no DNA, converts no results, and compares nothing.
  The harness is Rust.

## Supersession

[[ADR-0065]] is superseded in part:

- its rule that DNA contains no comparison logic and no parsers of external
  call sets or truth gives way to decision 2 above;
- its other parts still hold: DNA's public contracts carry what a downstream
  comparison needs, quoted metrics name the DNA revision they were measured
  against, and DNA's own tests assert DNA's behaviour against synthetic
  expectations.
