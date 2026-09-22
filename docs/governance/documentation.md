# Documentation Governance

DNA documentation is organized as a **knowledge system** with explicit authority,
lifecycle, ownership, and traceability. Folder names reflect the role of the
knowledge they contain; they are not a second copy of source layout.

## Navigation layers

| Layer | Purpose |
|---|---|
| repository `README.md` | product/use/contributor router |
| `docs/README.md` | canonical knowledge router |
| folder `README.md` files | area-specific router and authority boundary |
| source-directory `README.md` files | implementation ownership colocated with code |
| `AGENTS.md` | agent routing plus repository-wide invariants |

README files orient and link. They must not become duplicate specifications.

## Authority model

| Area | Owns | Class |
|---|---|---|
| `requirements/` | normative intended behavior and quality attributes | canonical / living |
| `architecture/` | stable system structure and cross-cutting invariants | canonical / living |
| `design/` | current mechanisms and algorithms | canonical / living |
| `reference/` | exact public/configuration/schema/coordinate semantics | canonical / living |
| source + tests | executable behavior and encoded evidence | executable reality |
| `validation/` | evidence strategy, acceptance, datasets, limits, traceability | canonical / living |
| `engineering/` | development/testing/CI/CD/release/dependency process | canonical / living |
| `operations/` | readiness, observability, procedures, investigation | canonical / living |
| `security/` | trust boundaries and security controls | canonical / living |
| `governance/` | lifecycle, ownership, versioning, data/document policy | canonical / living |
| `decisions/` | durable historical rationale | historical / durable |
| `proposals/` | reviewed changes under evolution | evolutionary |
| `research/` | evidence and experiments | exploratory |

A machine-readable schema is authoritative for the exact serialized shape of its
named version. Requirements govern intended behavior; design/reference explain
current mechanisms/interfaces. Accepted ADRs preserve rationale rather than
overriding newer current-state documents.

If source and normative current-state documentation disagree, surface the mismatch
and reconcile it in the change that owns the behavior. Do not silently choose the
artifact that is easiest to implement.

## Canonical-home rule

One fact has one authoritative home.

- what must be true -> requirements;
- stable system structure/invariants -> architecture;
- how current mechanisms work -> design;
- exact interfaces/shapes -> reference;
- why a durable choice was made -> decisions/ADR;
- what change is being considered -> proposals;
- what was learned but not accepted -> research;
- how correctness is demonstrated -> validation;
- how software is developed/released -> engineering;
- how it is run/investigated -> operations;
- how trust boundaries are protected -> security;
- how knowledge/process evolves -> governance.

Other documents link to the canonical home instead of copying its specification.

## Promotion path

```text
research evidence
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
engineering/release
      ↓
operations / production feedback
```

Research and proposals do not become production truth merely by existing or being
accepted. Current-state authorities, implementation, tests, and validation must
be updated explicitly.

## Source-local README rule

Every directory under `src/` must contain an up-to-date `README.md` describing
the responsibility boundary and routing to canonical documentation.

A source README may describe responsibility, non-responsibilities, entry points,
child modules, dependency direction, local invariants, and links. It must not
translate source line-by-line or duplicate full algorithm specifications.

File-only Rust modules use rustdoc/source comments. DNA does not maintain a
parallel `docs/src/` mirror.

## Change impact

Update only the authorities affected by the change:

- internal refactor with unchanged responsibility -> source README usually unchanged;
- ownership/dependency boundary -> source README + architecture as needed;
- scientific/algorithm behavior -> requirements + design + tests + validation;
- public schema/config/CLI -> requirements + reference + examples/tests;
- durable architectural/scientific/security choice -> ADR or ADR successor;
- proposed but not accepted change -> proposal only, plus research evidence;
- research-only work -> research only until promotion;
- build/test/release process -> engineering;
- runtime/support procedure -> operations;
- trust boundary/dependency vulnerability policy -> security.

Code and affected documentation should land in the same change.

## ADR lifecycle

Treat ADRs as an append-only decision log.

- search current requirements/architecture/design/reference and existing ADRs first;
- do not create an ADR for ordinary implementation progress or documentation cleanup;
- a materially changed choice creates a successor that explicitly supersedes the prior ADR;
- historical ADR identifiers are never reused;
- current truth belongs in living documentation, not in a rewritten historical ADR.

## Folder indexes and cleanup

- every documentation folder under `docs/` has one canonical `README.md`;
- every source directory under `src/` has one current `README.md`;
- do not keep sibling `topic.md` + `topic/README.md` entry points;
- do not keep compatibility pointers, legacy copies, or old/new parallel docs after migration;
- promoted temporary research/planning is deleted once no unique evidence remains;
- Git history is the archive for removed stale current-state documentation.

## Naming

Prefer stable responsibility names such as:

- `requirements/SRS.md`;
- `architecture/overview.md`;
- `design/pipeline.md`;
- `decisions/adr/`;
- `validation/strategy.md`;
- `engineering/ci-cd.md`;
- `operations/runbooks/`;
- `security/threat-model.md`;
- `reference/schemas/`.

Avoid catch-all or lifecycle-noise names such as `NOTES.md`, `NEW.md`,
`FINAL.md`, `legacy/`, or `old/`.

See [lifecycle](lifecycle.md), [ownership](ownership.md), and [versioning](versioning.md).
