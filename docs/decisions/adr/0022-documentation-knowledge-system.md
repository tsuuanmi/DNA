# ADR-0022: Govern documentation as an executable knowledge system

- **Status:** Accepted
- **Date:** 2026-09-22
- **Supersedes:** ADR-0006

## Context

DNA has multiple kinds of knowledge: normative requirements, architecture,
current mechanisms, public contracts, durable decisions, proposed changes,
research evidence, validation, engineering/release process, operations, security,
governance, implementation ownership, and agent instructions.

These artifacts do not have equal authority or lifecycle. Without an explicit
model, humans and coding agents can treat research as production truth, use an ADR
as a current design manual, duplicate specifications across folders, or let
implementation documentation drift away from source.

The goal is not to maximize documentation volume. The goal is one canonical home
per fact, explicit lifecycle, clear ownership, traceability, and executable
repository checks.

## Decision

DNA documentation is governed as a layered knowledge system.

### 1. Knowledge is separated by role

| Area | Responsibility | Lifecycle class |
|---|---|---|
| `requirements/` | normative intended behavior and quality attributes | canonical / living |
| `architecture/` | system structure and cross-cutting invariants | canonical / living |
| `design/` | current mechanisms and algorithms | canonical / living |
| `reference/` | public/configuration/schema/coordinate semantics | canonical / living |
| `validation/` | acceptance, evidence, limitations, traceability | canonical / living |
| `engineering/` | development, testing, CI/CD, release, dependencies | canonical / living |
| `operations/` | readiness, observability, runbooks, playbooks | canonical / living |
| `security/` | trust boundaries, threats, supply-chain policy | canonical / living |
| `governance/` | documentation/data/versioning/ownership/lifecycle policy | canonical / living |
| `decisions/` | durable decision rationale | historical / durable |
| `proposals/` | reviewed changes under consideration | evolutionary |
| `research/` | evidence, experiments, candidate ideas | exploratory |

Machine-readable schemas remain authoritative for the exact shape of their named
version. Source code is executable reality. A disagreement between source and
normative current-state documentation is a defect to reconcile, not a reason to
silently choose whichever artifact is convenient.

### 2. README files are routers

- root `README.md` routes users and contributors;
- `docs/README.md` routes knowledge by authority/lifecycle;
- every documentation folder has one canonical `README.md` index;
- every real source directory has a colocated `README.md` describing ownership;
- `AGENTS.md` contains routing plus repository-wide invariants.

Routers link to deeper canonical material rather than duplicating it.

### 3. Implementation documentation is colocated with source

Every directory under `src/` has an up-to-date `README.md` describing
responsibility, non-responsibilities, entry points, dependency direction, local
invariants, and links to relevant canonical docs/tests.

File-only Rust modules use rustdoc/source comments. DNA does not maintain a
parallel `docs/src/` shadow tree.

### 4. Current truth, history, change, and evidence are distinct

- requirements/architecture/design/reference describe current intended truth;
- ADRs explain why durable choices were made;
- proposals describe changes not yet current;
- research supplies evidence and hypotheses;
- validation demonstrates whether claims/requirements are supported.

Historical ADRs are not rewritten into current design manuals. Current-state docs
are not retained as legacy copies after replacement; Git preserves their history.

### 5. Change promotion is explicit

```text
research
  ↓
proposal / review
  ↓
decision when architecturally significant
  ↓
requirements + architecture + design + reference
  ↓
source + source README + tests
  ↓
validation
  ↓
release / operations
```

Accepted research/proposals do not become production behavior until the relevant
current-state authorities, implementation, tests, and validation are updated.

### 6. Traceability connects intent to evidence

`docs/validation/traceability.md` maps requirement families to architecture/design,
owning source, tests/evidence, and public reference contracts. It is a navigation
aid, not a duplicated specification.

### 7. ADRs are deduplicated by decision boundary

One accepted ADR is the canonical rationale for one durable decision scope.
Implementation progress, schema revisions, validation observations, or ordinary
documentation cleanup do not require new ADRs. A materially changed choice
creates a successor that explicitly supersedes the prior ADR in whole or in part.

### 8. Documentation structure is executable policy

CI validates:

- one README index per documentation folder;
- one README per source directory;
- no `docs/src/` shadow mirror;
- no duplicate `topic.md` + `topic/README.md` entry points;
- no legacy/archive/history/temp current-state paths;
- repository-local Markdown links resolve.

## Consequences

### Positive

- humans and agents can determine which artifact is authoritative;
- research and proposals can be rich without contaminating production semantics;
- decision history is preserved without becoming current-state duplication;
- source ownership is visible where developers browse code;
- navigation drift is caught by CI;
- requirement -> design -> source -> test/validation paths are explicit.

### Cost

- behavior changes may require coordinated updates across several authorities;
- module owners must maintain source-local README files when responsibilities move;
- lifecycle and link validation add repository-maintenance discipline.

## Non-goal

The taxonomy does not require empty enterprise-style folders. Documentation areas
are created when the system has a real artifact or governance boundary to place
there.

The repository-agnostic rationale and adoption guidance live in the reusable
[Documentation Architecture Standard](../../governance/documentation-architecture.md).
DNA-specific enforcement lives in the [DNA documentation policy](../../governance/documentation.md).
