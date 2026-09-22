# ADR-0022: Govern documentation as an executable knowledge system

## Status

Accepted

**Supersedes:** ADR-0006

## Context

DNA contains multiple kinds of knowledge: requirements, architecture, decisions,
scientific methods, public contracts, implementation ownership, validation,
operations, research, roadmap material, and agent instructions.

These artifacts do not have equal authority. A parallel source-documentation tree
also creates a second hierarchy that can drift from the code it describes.

The goal is therefore not to maximize documentation volume. The goal is to make
current truth easy to locate, hard to duplicate, and synchronized with the code.

## Decision

DNA documentation is governed as a layered knowledge system.

### 1. Current production knowledge is authoritative by role

- **SRS** defines what the current system MUST/SHOULD/MAY do.
- **Accepted ADRs** record why durable architectural or scientific choices were made.
- **Architecture** defines stable boundaries, dependency direction, data flow, and
  cross-cutting invariants.
- **Methods** define current scientific and algorithmic behavior.
- **Contracts** define user-visible and machine-visible interfaces.
- **Source-directory README files** define implementation ownership and route
  readers to relevant code and canonical documentation.
- **Validation, operations, and governance** define verification and operating policy.
- **Roadmap and research** are non-normative.

A disagreement between source and normative production documentation is a defect
to reconcile, not permission to choose whichever artifact is convenient.

### 2. README files are routers

README files orient readers; they do not become duplicate specifications.

- The repository `README.md` routes users and contributors into the project.
- `docs/README.md` routes readers to canonical knowledge by role.
- Each documentation-folder `README.md` indexes that knowledge area.
- Each source-directory `README.md` describes module responsibility and links to
  deeper SRS, architecture, methods, contracts, ADRs, and tests.

### 3. Implementation documentation is colocated with code

Every directory under `src/` has an up-to-date `README.md`.

The README describes:

- responsibility and non-responsibilities;
- boundary entry points;
- important child modules/files;
- dependencies and local invariants;
- links to canonical requirements, methods, decisions, contracts, and tests.

File-only modules continue to use rustdoc/module comments. They are not converted
into directories solely to create documentation.

DNA does not maintain a `docs/src/` shadow tree. This replaces ADR-0006's
one-to-one source/manual mirroring model with a lower-duplication, colocated
ownership model.

### 4. Cross-cutting invariants have one explicit home

Stable invariants spanning modules belong under
`docs/architecture/invariants/`. Source README files link to those invariants
rather than restating them inconsistently.

### 5. Traceability connects intent to implementation and evidence

`docs/validation/traceability.md` maps requirement families to current architecture/methods,
owning source modules, tests/evidence, and public contracts. It is a navigation
aid, not a second specification.

### 6. Research is explicitly non-normative

Exploratory work belongs under `docs/research/<topic>/`.

Promotion is explicit:

```text
research evidence
  ↓
accepted decision when needed
  ↓
current requirement / architecture / method / contract
  ↓
implementation + source-local README + tests
  ↓
validation
```

Research does not change production behavior merely by existing.

### 7. AGENTS.md is routing plus invariants

Root `AGENTS.md` remains concise and repository-agnostic. It tells an agent:

- how to find repository and documentation entry points;
- how to resolve authority by role;
- that the nearest source README defines implementation ownership;
- that code and affected documentation change together;
- which repository-wide implementation invariants must be preserved;
- how to discover verification from CI and operations docs.

Repository-specific scientific detail, commands, and contracts remain in their
canonical homes instead of being copied into `AGENTS.md`.

### 8. ADRs are deduplicated by decision boundary

One accepted ADR is the canonical rationale for one durable decision scope.

Implementation progress, schema revisions, validation observations, and ordinary
documentation cleanup do not require new ADRs. A material replacement creates a
successor that explicitly supersedes the prior decision in whole or in part.

## Consequences

### Positive

- implementation documentation is visible where developers browse the code;
- there is no second one-to-one source hierarchy to keep synchronized;
- README files remain compact navigation surfaces;
- agents can move from requirement to implementation without treating research as truth;
- source changes have an explicit documentation-impact rule;
- CI can enforce both documentation-folder indexes and source-directory README coverage.

### Cost

- module owners must keep colocated README files current when responsibilities change;
- documentation review remains necessary because structural checks cannot prove prose accuracy;
- migration removes the existing `docs/src/` manuals and their dedicated mirror test.

## Revision

2026-09-22: refined the accepted documentation-system decision by replacing the
ADR-0006 `docs/src/` one-to-one mirror with colocated source-directory README
routers and narrowing `AGENTS.md` to routing plus invariants.


## 2026-09-22 taxonomy amendment

The layered-knowledge decision remains unchanged, but the current role names were
refined to reduce ambiguity and make lifecycle explicit:

- `srs/` -> `requirements/` with `SRS.md`, quality attributes, constraints,
  and stable requirement-family documents;
- `methods/` -> `design/` for current mechanisms and algorithms;
- `contracts/` -> `reference/` for public/configuration/schema semantics;
- `adr/` -> `decisions/adr/`, with `decisions/` as the historical rationale layer;
- `proposals/` is the evolutionary layer between research and accepted current truth;
- `validation/` owns acceptance criteria, datasets/evidence, limitations, benchmarks,
  and requirement-to-implementation traceability;
- `engineering/` owns development/testing/CI/CD/release/dependencies/code quality;
- `operations/` owns readiness, observability, runbooks, and playbooks;
- `security/` owns trust boundaries, threat modeling, and supply-chain/vulnerability policy;
- `governance/` owns documentation lifecycle, ownership, versioning, and data policy.

The original wording above is retained as historical decision context. Current
navigation and authority are defined by [documentation governance](../../governance/documentation.md)
and [docs/README.md](../../README.md).
