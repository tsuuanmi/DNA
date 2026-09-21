# AGENTS.md — Agent Development Guide

This file defines the default workflow for coding agents working in this repository.

The guide is intentionally repository-agnostic. Do not assume a particular documentation layout, language, package manager, test framework, or file naming convention. Discover the repository's own structure and authoritative instructions before changing anything.

## 1. Main workflow

Follow this workflow for every non-trivial task:

```text
discover
  ↓
understand
  ↓
plan
  ↓
implement
  ↓
verify
  ↓
reconcile
  ↓
review
  ↓
report
```

### Discover

Before editing:

- inspect repository status and current branch/worktree state;
- find repository-level instructions and contributor guidance;
- identify the build/package system, CI configuration, test layout, and documentation entry points;
- locate the files and modules that appear to own the requested behavior;
- notice unrelated or in-progress changes and avoid touching them.

Do not infer the repository structure from conventions used in another project.

### Understand

Read enough context to understand both **intent** and **implementation** before changing code.

For the area being changed, discover and read the repository's equivalents of:

- product or software requirements;
- architecture and cross-cutting invariants;
- accepted design/decision records;
- domain or algorithm documentation;
- public API, schema, configuration, CLI, or protocol contracts;
- implementation/module documentation;
- relevant source files;
- relevant tests and fixtures;
- validation, release, or operational requirements;
- roadmap/research material when useful for context.

Read affected source files in full for architectural, scientific, security-sensitive, or broad refactoring work. Search snippets alone are not sufficient for understanding a module.

### Plan

Before implementation, identify:

- the intended behavior;
- what must remain unchanged;
- the smallest coherent change;
- affected source, tests, contracts, docs, generated artifacts, configuration, and dependencies;
- the invariants and public behavior at risk;
- the verification required to demonstrate correctness.

Do not expand scope merely because nearby cleanup is possible.

### Implement

Make the smallest production-quality change that satisfies the task.

Prefer:

- one authoritative implementation;
- clear ownership and dependency boundaries;
- explicit states and failures;
- existing project patterns and utilities;
- deterministic behavior where the domain requires it;
- removal of obsolete logic when the current specification no longer requires compatibility.

Avoid:

- speculative abstractions;
- duplicate implementations;
- fallback paths added only to preserve old behavior;
- unrelated refactors;
- silent semantic changes;
- unnecessary dependencies;
- broad formatting churn.

### Verify

Choose verification based on the failure modes introduced by the change.

Discover the repository's required commands from CI, build configuration, contributor docs, or project scripts rather than inventing them.

Possible verification layers include:

- formatting;
- compilation/build;
- lint/static analysis;
- type checking;
- unit tests;
- integration tests;
- property/invariant tests;
- schema/contract validation;
- fuzz/adversarial tests;
- performance/resource checks;
- dependency/security audits;
- domain-specific or real-data validation.

Run the narrowest checks that provide sufficient confidence during iteration, then run all required gates before finalizing when the repository policy requires them.

A green test suite proves only the properties encoded by those tests. It does not automatically prove domain correctness.

### Reconcile

After the implementation is correct, synchronize every affected representation of the behavior.

Depending on the repository, this may include:

- requirements;
- architecture/invariants;
- design decisions;
- public contracts or schemas;
- module/implementation documentation;
- examples;
- tests and fixtures;
- generated artifacts;
- changelog/release notes;
- traceability or validation records.

Do not update documentation mechanically. Update the authoritative layer that actually changed.

### Review

Before finalizing:

- inspect the complete diff;
- confirm every changed line belongs to the task;
- verify no obsolete or parallel implementation remains unintentionally;
- check imports, exports, references, types, tests, configuration, and docs for stale dependencies;
- confirm no accidental lockfile, generated-file, or formatting churn;
- ensure public behavior and documented intent agree;
- ensure the change fits the repository's architecture rather than merely passing tests.

### Report

The final response should state:

- what changed;
- important files/modules affected;
- verification performed and its result;
- tests or checks intentionally not run and why;
- documentation/changelog impact;
- remaining risks, assumptions, or follow-up work.

Keep the report concise and factual.

## 2. How to read repository documentation

Do not rely on fixed filenames. Infer each document's **role**.

Typical roles include:

### Normative requirements

Describe what the current system is intended to do.

Examples of forms this may take:

- SRS;
- product requirements;
- specification;
- protocol requirements;
- acceptance criteria.

Treat normative language such as MUST/SHOULD/MAY according to the repository's documented convention.

### Public contracts

Define interfaces that users, tools, or other systems depend on.

Examples:

- API definitions;
- schemas;
- CLI behavior;
- configuration formats;
- file formats;
- database contracts;
- protocol definitions;
- serialization/versioning rules.

Machine-readable contracts may be authoritative for exact syntax or shape while prose explains semantics.

### Architecture and invariants

Describe:

- component boundaries;
- dependency direction;
- ownership;
- data flow;
- lifecycle;
- trust boundaries;
- coordinate/unit conventions;
- states that must always remain true.

Use these to decide **where** a change belongs, not just what code currently happens to do.

### Decision records

Explain why a non-obvious design choice exists, alternatives considered, tradeoffs, and consequences.

The filename or format may be ADR, RFC, design note, proposal, decision log, or something else.

Check status and supersession before relying on a decision.

### Domain and method documentation

Explains domain semantics, algorithms, assumptions, limitations, and interpretation.

This is especially important in scientific, financial, security, distributed, data-processing, and protocol-heavy systems where technically valid code can still model the domain incorrectly.

### Implementation documentation

Explains module responsibility, inputs/outputs, dependencies, failure modes, internal invariants, and ownership.

Implementation docs are descriptive. They should not silently override a normative requirement or public contract.

### Validation and operations

Describe how correctness, performance, security, release quality, migration, deployment, or real-world behavior is demonstrated.

Do not confuse "the code builds" with "the system satisfies its domain or operational requirements."

### Roadmap and research

Describe possible future work, experiments, hypotheses, or proposed architecture.

Treat exploratory material as non-normative unless the repository explicitly promotes it into the current specification.

## 3. Resolve authority by role, not by filename

When artifacts disagree, do not silently choose whichever one is convenient.

Determine:

1. which artifact is intended to be normative for that concern;
2. whether a newer decision supersedes an older one;
3. whether source code is ahead of docs or the implementation is incomplete;
4. whether the mismatch is a bug, stale documentation, or an intentional transition.

Useful distinction:

```text
requirements/contracts = intended behavior
source                 = current executable behavior
decision records       = rationale
architecture           = responsibility and invariants
implementation docs    = current ownership/details
validation             = evidence
roadmap/research       = possible future behavior
```

A mismatch between intended behavior and executable behavior must be surfaced and resolved deliberately.

## 4. Invariants first

Before changing a mature system, identify the invariants affected by the task.

Examples:

- identity must be preserved across transformations;
- coordinate/unit domains must not be mixed;
- source evidence must remain immutable;
- operations must be atomic;
- public schemas are versioned;
- invalid or unresolved states must remain explicit;
- resource use must be bounded;
- deterministic inputs must yield deterministic outputs;
- security/trust boundaries must remain intact.

When practical, move important invariants from comments and developer memory into types, constructors, APIs, validators, schemas, static analysis, or tests.

Do not create type-level or abstraction complexity unless it prevents a concrete failure mode.

## 5. Domain-sensitive changes

For systems where domain correctness matters, review a change through three questions:

### Evidence / inputs

- What was actually measured, received, or observed?
- What transformations are applied?
- Which information is original and which is derived?
- Are assumptions and units explicit?

### Domain meaning

- What claim does the result actually support?
- Which states should remain unknown, unresolved, or ambiguous?
- What alternative explanations remain possible?
- Is the implementation accidentally making a stronger claim than the evidence supports?

### Engineering

- Which invariant or contract changes?
- Which component owns the behavior?
- How are failures represented?
- What tests or validation demonstrate correctness?
- Does the change require a contract/version/migration decision?

A technically elegant implementation is not sufficient if it models the wrong domain concept.

## 6. Research and proposal promotion

Exploratory work should not silently become production behavior.

A common promotion path is:

```text
research / proposal
        ↓
accepted decision
        ↓
normative requirement / public contract
        ↓
implementation
        ↓
tests
        ↓
validation
        ↓
release
```

Repositories may use different names, but preserve the distinction between **exploration** and **current production truth**.

Agreement with another implementation is comparison evidence, not automatically ground truth.

## 7. Tests and validation

Choose tests according to the risk being protected.

Examples:

- parser/input boundary -> malformed/adversarial cases;
- state/coordinate transform -> invariant/property tests;
- API/schema -> contract tests;
- bug fix -> focused regression test;
- concurrency -> race/interleaving tests where supported;
- performance-sensitive logic -> benchmark/resource measurements;
- security-sensitive boundary -> threat-specific tests/audits;
- scientific/domain method -> independent or real-world validation where required.

Coverage is not a substitute for meaningful assertions.

Mutation, fuzzing, property testing, static analysis, formal methods, or sanitizers are useful only when they address an identified failure class.

## 8. Source and documentation synchronization

A behavior change may require updates across multiple layers.

Use the repository's own structure, but think in terms of impact:

| Change | Consider updating |
|---|---|
| domain/algorithm behavior | requirement + method/domain docs + tests + validation |
| architectural boundary | architecture/invariants + decision record + implementation docs |
| public API/schema/config/CLI | contract + examples + tests + migration/changelog |
| internal ownership | module/implementation docs |
| dependency/toolchain | lockfiles + build/release/security docs |
| research only | research/proposal area only until promoted |

Do not duplicate the same specification across many documents. Link to the authoritative source instead.

## 9. Dependencies and generated artifacts

Treat dependency and lockfile changes as code changes.

Before adding a dependency:

- confirm the capability is not already available;
- justify the new dependency and its maintenance/security cost;
- follow the repository's package-management policy;
- review the resulting lockfile diff.

Do not edit generated artifacts directly when an authoritative generator or source definition exists. Modify the source, regenerate, and review the generated diff.

## 10. Worktree and version-control safety

Assume other work may exist in the same repository.

- inspect status before editing;
- avoid destructive repository-wide reset/clean operations;
- do not overwrite unrelated work;
- stage or commit only intentional files;
- do not hide unrelated failures;
- do not create temporary backups inside the repository unless repository policy explicitly requires it;
- never bypass required hooks or verification merely to make a change mergeable.

If a conflict involves unrelated work you do not understand, stop rather than overwriting it.

## 11. Completion standard

A task is complete when:

- the requested behavior is implemented;
- the implementation fits existing architectural boundaries;
- affected invariants and public contracts are preserved or intentionally revised;
- relevant tests and required checks pass;
- documentation reflects the intended behavior;
- obsolete paths introduced or superseded by the task are removed;
- the final diff contains no unrelated changes;
- remaining uncertainty is stated explicitly.

Passing tests alone is not the completion criterion.
