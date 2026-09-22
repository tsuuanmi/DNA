# Documentation and Code Ownership

Ownership follows responsibility, not document age.

## Canonical knowledge ownership

| Area | Owns |
|---|---|
| `requirements/` | normative intended behavior and quality attributes |
| `architecture/` | stable system structure and cross-cutting invariants |
| `design/` | current mechanisms and algorithms |
| `reference/` | public/configuration/schema/coordinate semantics |
| `decisions/` | durable historical rationale |
| `proposals/` | reviewed changes not yet current truth |
| `validation/` | evidence strategy, acceptance, datasets, limitations, traceability |
| `engineering/` | development, testing, CI/CD, release engineering, dependencies |
| `operations/` | runtime procedures, readiness, observability, troubleshooting |
| `security/` | trust boundaries, threats, supply-chain/vulnerability policy |
| `research/` | non-normative evidence and experiments |
| `governance/` | lifecycle, ownership, versioning, data/documentation policy |

## Source ownership

Every real source directory contains a colocated `README.md` describing responsibility, boundary, and navigation. File-local API detail belongs in rustdoc/source comments.

Do not create a second source documentation tree.

## Change ownership

The change that modifies behavior owns the corresponding documentation/test update. Do not defer required documentation synchronization to a follow-up without explicitly recording the inconsistency.
