# DNA Documentation Policy

DNA adopts the reusable [Documentation Architecture Standard](documentation-architecture.md).
This document records the DNA-specific application of that standard and the
repository rules enforced by CI.

## Current DNA knowledge areas

| Area | DNA-specific responsibility |
|---|---|
| [requirements](../requirements/README.md) | normative DNA SRS, quality attributes, and constraints |
| [architecture](../architecture/README.md) | DNA system context, module boundaries, flow, interfaces, and invariants |
| [design](../design/README.md) | current scientific and algorithmic mechanisms |
| [decisions](../decisions/README.md) | durable DNA ADR history |
| [proposals](../proposals/README.md) | DNA changes under review plus the non-normative roadmap |
| [research](../research/README.md) | exploratory scientific/engineering evidence |
| [validation](../validation/README.md) | acceptance, datasets/evidence, limitations, benchmarks, and traceability |
| [engineering](../engineering/README.md) | development, testing, CI/CD, release, dependencies, and code quality |
| [operations](../operations/README.md) | readiness, observability, batch runbook, and investigation playbook |
| [security](../security/README.md) | local-CLI trust boundaries, threat model, and supply-chain policy |
| [reference](../reference/README.md) | JSON schemas, result/configuration semantics, coordinates, examples, glossary |
| [governance](README.md) | documentation/data/versioning/ownership/lifecycle policy |

DNA does not create empty service-only documentation such as deployment topology,
service SLOs, incident/postmortem collections, or disaster recovery while those
boundaries do not exist.

## Navigation

- root `README.md` is the product/contributor router;
- `docs/README.md` is the knowledge router;
- every documentation folder has one `README.md` index;
- every directory under `src/` has a colocated `README.md`;
- `AGENTS.md` contains agent routing plus repository-wide invariants.

README files do not duplicate full specifications.

## Authority

For DNA:

- requirements define intended behavior;
- architecture/design/reference define current structure, mechanisms, and exact
  interfaces;
- source is executable reality;
- tests and validation provide evidence;
- accepted ADRs preserve rationale but do not override newer current-state docs;
- proposals and research are non-authoritative until explicitly promoted.

If source and normative current-state documentation disagree, treat that as a
defect to reconcile in the owning change.

## Source-local documentation

Every `src/**/` directory has an up-to-date `README.md` describing its
responsibility boundary and navigation.

File-only Rust modules use rustdoc/source comments. DNA does not maintain a
`docs/src/` mirror.

## DNA promotion flow

```text
research
  ↓
proposal / review
  ↓
ADR when durable rationale is needed
  ↓
requirements + architecture + design + reference
  ↓
source + source README + tests
  ↓
validation
  ↓
engineering / release
  ↓
operations feedback
```

Not every change needs every stage. The affected canonical layers must still be
updated atomically with the implementation.

## Change impact

- internal refactor without boundary/behavior change -> source/tests only as needed;
- module responsibility/dependency change -> source README + architecture;
- scientific/algorithm behavior change -> requirements + design + tests + validation;
- public schema/config/CLI change -> requirements + reference + examples/tests;
- durable architectural/scientific/security choice -> ADR or successor ADR;
- unaccepted change -> proposal;
- exploratory work -> research;
- build/test/release process -> engineering;
- runtime/support procedure -> operations;
- trust-boundary or vulnerability policy -> security.

## ADR policy

DNA ADRs are append-only decision history.

- search existing current-state docs and ADRs before creating another record;
- do not create ADRs for ordinary implementation progress;
- create a successor for a materially changed decision;
- mark supersession in both directions where applicable;
- never reuse an ADR identifier.

[ADR-0022](../decisions/adr/0022-documentation-knowledge-system.md) is the owning
decision for this documentation system.

## CI enforcement

`scripts/validate_docs_structure.py` enforces:

- one README per documentation folder;
- one README per source directory;
- no `docs/src/` shadow tree;
- no duplicate `topic.md` + `topic/README.md` entry points;
- no repository-defined legacy/archive/history/temp current-state paths;
- repository-local Markdown links resolve.

Other CI gates validate source policy, formatting, static analysis, tests, schemas,
configuration, and reference identity.

Structural validation cannot prove prose accuracy. Reviewers still own semantic
correctness.

## Cleanup

Do not keep renamed compatibility pointers, old/current parallel docs, or a generic
legacy archive after migration. Git preserves history.

When research or temporary planning is fully promoted, delete duplicated text
unless it contains unique evidence worth retaining.

See also [lifecycle](lifecycle.md), [ownership](ownership.md), and
[versioning](versioning.md).
