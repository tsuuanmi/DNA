# Documentation Architecture Standard

This document defines a reusable documentation architecture for production-grade
software repositories. It is intentionally technology- and product-agnostic so
other repositories can adopt the model without copying DNA-specific requirements,
algorithms, commands, or operational assumptions.

The goal is not a large `docs/` tree. The goal is a **knowledge system** in which
a reader or coding agent can answer:

1. What must be true?
2. What evidence or alternatives have been explored?
3. What change is being proposed?
4. Why was a durable choice made?
5. How is the current system structured?
6. How does the current mechanism work?
7. What interface or contract is exact?
8. Where is the behavior implemented?
9. How is the behavior verified?
10. How is the software built, released, operated, secured, and evolved?

## Core principles

### One fact, one canonical home

A fact should have one authoritative location. Other documents link to it instead
of restating it.

Duplication makes documentation drift independently from the code and from itself.
A small set of current, owned documents is more valuable than a larger set of
partially stale documents.

### Documentation changes with the system

When behavior, contracts, architecture, ownership, operations, or validation
expectations change, the affected documentation changes in the same change set.

Documentation is therefore part of the correctness system, not an after-the-fact
description.

### Current truth, history, change, and exploration are different things

Do not mix these knowledge classes:

- **canonical / living** — what is true now;
- **historical / durable** — why a decision was made or what happened;
- **evolutionary** — a change being proposed or implemented;
- **exploratory** — evidence, experiments, or ideas not yet accepted;
- **executable reality** — source and tests.

A proposal is not a specification. An ADR is not a current design guide. Research
is not production authority.

## Reference architecture

A mature repository may use this shape:

```text
docs/
├── README.md
├── requirements/
├── architecture/
├── design/
├── decisions/
│   └── adr/
├── proposals/
├── research/
├── validation/
├── engineering/
├── operations/
│   ├── runbooks/
│   ├── playbooks/
│   └── incidents/        # only when incidents are a real lifecycle artifact
├── security/
├── reference/
└── governance/
```

This is a **semantic taxonomy**, not a mandatory folder checklist. Create a
category only when the repository has real knowledge that belongs there.

For example, a local CLI should not invent service SLO, deployment-topology, or
disaster-recovery documents merely to resemble a distributed service repository.

## Knowledge areas

### `requirements/` — what must be true

Owns normative functional requirements, non-functional quality attributes,
constraints, and acceptance-level intent.

Good contents:

- SRS or requirement-family documents;
- quality attributes;
- platform/domain constraints.

Does not own:

- implementation mechanism;
- decision history;
- test results.

### `architecture/` — where responsibilities belong

Owns current system structure, context, boundaries, dependency direction, data
flow, interfaces between major components, and cross-cutting invariants.

Architecture answers **where and under what constraints**, not line-by-line
implementation details.

### `design/` — how current mechanisms work

Owns current algorithms, component mechanisms, data-model mechanics, pipeline
behavior, and other implementation-facing design that remains true after the
proposal phase.

Design is living current-state knowledge.

### `decisions/` — why durable choices were made

Owns append-only decision history, commonly ADRs.

An accepted ADR records context, alternatives, outcome, tradeoffs, and
consequences. If the decision changes, create a successor that supersedes the
prior record. Do not rewrite history to look current.

### `proposals/` — what change is being considered

Owns reviewable changes that are not yet current truth.

A proposal should make the problem, goals, non-goals, evidence, proposed design,
alternatives, validation, rollout/rollback, operational impact, security impact,
and implementation status explicit.

Typical lifecycle:

```text
draft -> proposed -> accepted -> implementing -> implemented
                  \-> rejected
implemented/accepted -> superseded
```

### `research/` — what has been learned

Owns non-normative experiments, external-system analysis, benchmarks under
investigation, scientific exploration, and candidate approaches.

Research can support a proposal or decision but does not change production
behavior on its own.

### `validation/` — how claims are demonstrated

Owns validation strategy, acceptance criteria, datasets/evidence, benchmark
evidence, known limitations, and requirement-to-implementation traceability.

Tests are evidence, but validation is broader than test execution when domain,
safety, performance, scientific, or operational claims require independent
evidence.

### `engineering/` — how software is built and delivered

Owns developer setup, testing strategy, CI/CD, code-quality policy, dependency
policy, release engineering, and tooling expectations.

This area describes the software-development system rather than runtime support.

### `operations/` — how the released system is supported

Owns production/readiness checks, observability, routine procedures, investigation
guides, and incident learning where those concerns exist.

Distinguish:

- **runbook** — a known procedure to achieve a known outcome;
- **playbook** — a guided investigation when the symptom is known but the root
  cause is not;
- **postmortem** — durable learning from a real incident.

Do not create incident/SLO/DR documentation before the repository actually has
those operational boundaries.

### `security/` — how trust boundaries are protected

Owns threat modeling, security architecture, sensitive-data rules,
dependency/supply-chain controls, vulnerability handling, and secrets policy when
secrets exist.

Security documentation should follow real trust boundaries, not a generic
checklist.

### `reference/` — what is exact

Owns precise stable lookup material: APIs, schemas, file formats, configuration
semantics, coordinate conventions, CLI/reference contracts, examples, and
glossaries.

Machine-readable schemas can be authoritative for exact shape while prose explains
semantics that the schema cannot express.

### `governance/` — how the knowledge system evolves

Owns documentation policy, lifecycle, ownership, versioning, naming, data policy,
and other repository-wide rules.

Governance should be small and stable. Feature-specific behavior does not belong
here.

## Implementation documentation belongs with implementation

Do not create a hand-maintained `docs/src/` shadow tree.

For a real package/module/component directory, colocate a small `README.md` that
answers:

- what this boundary owns;
- what it intentionally does not own;
- primary entry points;
- important child components;
- dependency direction and local invariants;
- where canonical requirements/design/reference/decisions/tests live.

File-local API/item detail belongs in the language-native documentation system
(for example rustdoc, Javadoc, docstrings, or generated API docs).

The exact granularity is repository-specific. Do not turn every source file into a
directory merely to attach a README.

## README files are routers

A README should orient, not duplicate.

### Repository README

Answers:

- what is this project;
- simplest use case;
- how to build/run/test at a high level;
- where should a user, contributor, or agent read next.

### `docs/README.md`

Routes by knowledge role and authority.

### Folder README

Explains:

- what belongs in this area;
- what does not belong here;
- which documents are canonical entry points;
- lifecycle/authority where ambiguity is possible.

### Source-directory README

Routes from implementation ownership to canonical project knowledge.

## AGENTS.md is routing plus invariants

Agent instructions should not become a second project encyclopedia.

Keep:

- repository/documentation navigation;
- authority-resolution rules;
- repository-wide invariants;
- the rule that code changes update affected docs;
- how verification is discovered;
- final-review expectations.

Keep detailed product knowledge, algorithms, commands, schemas, and operational
procedures in their canonical homes.

## Knowledge lifecycle

The architecture supports a closed learning loop:

```text
need / problem
      ↓
requirements and constraints
      ↓
research evidence
      ↓
proposal / review
      ↓
decision when durable rationale is needed
      ↓
architecture + design + reference
      ↓
source + tests
      ↓
validation
      ↓
engineering / release
      ↓
operations / incidents / feedback
      ↓
new requirement or research
```

Not every change needs every stage. A small bug fix may move directly from current
requirements/design to source/tests. A new architectural feature may use the full
flow.

The important rule is that a lower-authority artifact cannot silently become
production truth.

## Traceability

For important behavior, a reader should be able to navigate forward:

```text
requirement
  -> proposal/decision when applicable
  -> architecture/design/reference
  -> owning source
  -> tests/validation evidence
```

and backward:

```text
source behavior
  -> design/reference
  -> decision rationale
  -> requirement
```

Traceability is a navigation aid, not a duplicated specification.

## Lifecycle and durability

### Canonical / living

Examples: requirements, architecture, design, reference, validation policy,
engineering, operations, security, governance.

Update them when current truth changes.

### Historical / durable

Examples: ADRs, postmortems, implemented/rejected proposals with durable rationale,
important retained validation reports.

Preserve history; do not edit it into current-state documentation.

### Evolutionary

Proposals carry explicit status and implementation links.

### Exploratory

Research remains non-normative until promotion.

### Temporary

Scratch notes, migration checklists, one-off planning, and intermediate
investigation output should be promoted into a durable canonical/historical home
or deleted.

Git history is the archive; do not create a generic `legacy/` graveyard.

## Lightweight metadata

Use structured metadata when lifecycle or tooling benefits from it, for example:

```yaml
---
id: PROP-0017
type: proposal
status: implementing
owners: [core]
created: 2026-09-22
related-requirements: [REQ-042]
related-decisions: [ADR-0054]
implementation: [PR-48]
---
```

Do not add front matter mechanically to documents that have no useful lifecycle
state.

## Repository checks

Where practical, make documentation structure executable policy.

Useful CI checks include:

- every documentation directory has one README index;
- every selected source/package directory has a colocated README;
- forbidden shadow documentation trees do not reappear;
- duplicate entry points such as `topic.md` and `topic/README.md` are rejected;
- repository-local Markdown links resolve;
- generated schemas/examples validate;
- stale/legacy/temporary current-state paths are rejected where the repository
  deliberately forbids them.

CI cannot prove prose correctness, so review ownership remains necessary.

## Adoption in another repository

Start small.

### Minimum

```text
README.md
AGENTS.md                    # if coding agents are used
docs/
├── README.md
├── requirements/
├── architecture/
├── decisions/
└── engineering/
```

Add `design/`, `reference/`, `validation/`, `proposals/`, `research/`,
`operations/`, `security/`, or `governance/` when real artifacts need those
homes.

### Migration sequence

1. inventory existing docs and identify duplicated facts;
2. classify each document by role and lifecycle;
3. choose one canonical home per fact;
4. move current truth before deleting duplicate/legacy paths;
5. preserve durable decision/incident history;
6. add README routers;
7. colocate package/module ownership docs with code;
8. update links and traceability;
9. add structural/link checks to CI;
10. delete temporary compatibility pointers once migration is complete.

Do not preserve old paths merely for documentation compatibility unless the
repository explicitly requires published stable documentation URLs.

## Why this model

This architecture combines several mature practices:

- Google documentation guidance emphasizes small fresh docs, updating docs with
  code, README-based orientation, deleting dead documentation, and avoiding
  duplicated guides:
  <https://google.github.io/styleguide/docguide/best_practices.html>
- Microsoft Azure Well-Architected guidance treats ADRs as an append-only decision
  log and recommends new records that supersede old decisions instead of editing
  accepted history:
  <https://learn.microsoft.com/en-us/azure/well-architected/architect-role/architecture-decision-record>
- Swift Evolution separates proposal review from implementation and maintains
  explicit proposal states:
  <https://github.com/swiftlang/swift-evolution/blob/main/process.md>
- Kubernetes KEPs use explicit lifecycle metadata and production-readiness review
  for changes intended for production:
  <https://github.com/kubernetes/enhancements/blob/master/keps/sig-architecture/0000-kep-process/README.md>
- AWS Well-Architected separates operational readiness, runbooks for known
  procedures, and playbooks for investigation:
  <https://docs.aws.amazon.com/wellarchitected/latest/operational-excellence-pillar/operational-readiness.html>

These are influences, not a claim that any one organization uses this exact
directory tree.

## Mental model

> Requirements define intent.  
> Research provides evidence.  
> Proposals explore change.  
> Decisions preserve rationale.  
> Architecture defines structure.  
> Design defines mechanisms.  
> Reference defines exact interfaces.  
> Code realizes the design.  
> Tests and validation demonstrate behavior.  
> Engineering delivers it.  
> Operations supports it.  
> Security protects its trust boundaries.  
> Governance keeps the knowledge system coherent.

The directory tree is only the storage layout. The lifecycle and authority model
are the architecture.
