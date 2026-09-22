# DNA Documentation

This is the canonical knowledge router for humans and coding agents.

DNA documentation implements the reusable [Documentation Architecture Standard](governance/documentation-architecture.md). It is organized by **role, authority, lifecycle, and ownership**. One fact has one canonical home; other documents link to it rather than duplicating it.

## Read order for a change

1. [Requirements](requirements/README.md) — what must be true.
2. [Architecture](architecture/README.md) — where responsibilities and invariants belong.
3. [Design](design/README.md) — how current mechanisms work.
4. [Decisions](decisions/README.md) — why durable choices were made.
5. [Proposals](proposals/README.md) and [research](research/README.md) — changes/evidence that are not current truth.
6. [Reference](reference/README.md) — exact public/configuration/schema semantics.
7. The nearest [source README](../src/README.md) and affected source — implementation ownership and executable behavior.
8. [Validation](validation/README.md) — evidence that protects the behavior.
9. [Engineering](engineering/README.md), [operations](operations/README.md), and [security](security/README.md) — delivery and production support.
10. [Governance](governance/README.md) — lifecycle, ownership, versioning, and documentation policy.

## Authority and durability

| Area | Role | Class |
|---|---|---|
| [Requirements](requirements/README.md) | normative intended behavior | canonical / living |
| [Architecture](architecture/README.md) | stable structure and invariants | canonical / living |
| [Design](design/README.md) | current mechanisms | canonical / living |
| [Reference](reference/README.md) | exact interfaces/contracts | canonical / living |
| source + tests | executable behavior/evidence | executable reality |
| [Validation](validation/README.md) | acceptance/evidence policy | canonical / living |
| [Engineering](engineering/README.md) | development/delivery process | canonical / living |
| [Operations](operations/README.md) | runtime/readiness/support procedures | canonical / living |
| [Security](security/README.md) | trust boundaries and controls | canonical / living |
| [Governance](governance/README.md) | lifecycle/ownership/versioning policy | canonical / living |
| [Decisions](decisions/README.md) | historical rationale | historical / durable |
| [Proposals](proposals/README.md) | reviewed change under evolution | evolutionary |
| [Research](research/README.md) | evidence/experiments | exploratory |

If source and normative current-state documentation disagree, surface and reconcile the mismatch. Do not silently choose whichever artifact is convenient.

## Directory layout

```text
docs/
├── README.md
├── requirements/
├── architecture/
│   └── invariants/
├── design/
├── decisions/
│   └── adr/
├── proposals/
├── research/
├── validation/
├── engineering/
├── operations/
│   ├── runbooks/
│   └── playbooks/
├── security/
├── reference/
└── governance/
```

Categories are created when they have real artifacts or governance value. DNA intentionally does not create empty service-only documentation such as deployment topology, service SLOs, or disaster recovery while those system boundaries do not exist.

## Core principle

> Specs define truth. Research provides evidence. Proposals explore change. ADRs preserve decisions. Architecture describes structure. Design describes mechanisms. Code realizes design. Tests and validation prove behavior. Operations keep it supportable. Learning feeds new requirements and research.

See the reusable [documentation architecture standard](governance/documentation-architecture.md), the [DNA documentation policy](governance/documentation.md), and [lifecycle](governance/lifecycle.md).
