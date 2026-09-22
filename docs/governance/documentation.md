# Documentation Governance

DNA documentation is organized by authority and responsibility, not by file age,
document length, or a shadow copy of the source tree.

## Documentation architecture

DNA uses four navigation layers with different responsibilities:

| Layer | Purpose |
|---|---|
| repository `README.md` | product and contributor router |
| `docs/README.md` and folder `README.md` files | canonical knowledge routers |
| source-directory `README.md` files | implementation ownership routers colocated with code |
| `AGENTS.md` | agent routing and repository-wide change invariants |

README files orient the reader and link to deeper canonical material. They must not
become duplicate specifications.

## Authority model

| Layer | Purpose | Normative? |
|---|---|---|
| Requirements (`requirements/`) | intended system requirements | Yes |
| machine-readable schemas / configuration contract | exact machine-visible contract for the named version | Yes |
| accepted ADR | durable decision and rationale | Yes for the decision it governs |
| methods | detailed current scientific/algorithmic semantics | Yes |
| architecture / invariants | stable boundaries, dependency direction, and cross-cutting truths | Yes |
| source | behavior executed by the current revision | Executable reality |
| source-directory README | implementation ownership, boundary, and navigation | Descriptive; must track source |
| validation / operations / governance | verification and operational policy | Yes for the policy it governs |
| roadmap | future direction | No |
| research | exploration and candidate designs | No |

If source and normative documentation disagree, do not silently choose one.
Surface the mismatch as a defect, incomplete implementation, or stale
documentation and resolve it in the change that owns the behavior.

## Canonical-home rule

One fact should have one authoritative home.

- Requirements belong in the SRS.
- Durable rationale belongs in ADRs.
- Current algorithms belong in methods.
- Public/configuration/serialization semantics belong in contracts.
- Stable system boundaries and cross-cutting invariants belong in architecture.
- Implementation ownership belongs beside code in source-directory README files.
- Operational procedures belong in operations.
- Exploratory evidence belongs in research.

Routers and related documents link to the authoritative source instead of copying
its content.

## Promotion path

Research becomes production behavior only through explicit promotion:

```text
research evidence
        ↓
proposal/decision when a durable choice is required
        ↓
current requirement / architecture / method / contract
        ↓
source + source-local README + tests
        ↓
validation and release evidence
```

A research note is never production authority merely because it exists. Once
exploratory content becomes current production truth, move that truth into its
persistent authoritative home and delete the temporary duplicate.

## Source-local README rule

Every directory under `src/` must contain an up-to-date `README.md`.

A source-directory README is a compact router for a real implementation boundary.
It should describe or link to:

- responsibility and non-responsibilities;
- public or crate-level entry points;
- important child files/modules;
- dependency direction and local invariants;
- relevant requirements, architecture, methods, contracts, ADRs, and tests.

It must not translate source line by line or duplicate the full algorithm
specification.

File-only Rust modules do not need to be converted into directories merely to gain
a README. Use rustdoc/module comments for file-local API documentation.

DNA does not maintain a parallel `docs/src/` tree.

## Change impact

A code change updates only the documentation layers it actually affects.

- Internal refactor with unchanged responsibility: source README usually unchanged.
- Module responsibility/dependency change: update the nearest source README.
- Scientific behavior change: update requirements, design docs, tests, validation
  implications, and the source README when its boundary changes.
- Schema/config/CLI change: update machine and human contracts, SRS, tests, and examples.
- New architectural dependency or boundary: update architecture and the owning ADR
  when the durable decision changes.
- Research-only work: keep it under `docs/research/<topic>/` until promotion.

Code and affected documentation should land in the same change.

## ADR lifecycle and deduplication

Treat ADRs as a decision log, not a changelog.

- Search the SRS, architecture/invariants, methods/contracts, and ADR index before
  creating a new ADR.
- If the decision boundary and core choice are unchanged, update the owning
  current-state documentation and amend the existing ADR only when clarification
  of the rationale is useful.
- Create a successor ADR only for a material new choice, then mark the previous ADR
  `Superseded` or `Superseded in part` and link the relationship both ways.
- Never keep two accepted ADRs claiming authority over the same decision scope.
- Superseded ADRs remain as historical provenance; numbers are never reused.

## Folder indexes and cleanup

- Every documentation folder under `docs/` must contain one canonical
  `README.md` index.
- Every source directory under `src/` must contain one current `README.md`.
- A folder README is a router; do not keep a sibling `<folder>.md` as a second index.
- Do not keep compatibility pointers, legacy copies, renamed duplicates, or
  old/new parallel documentation paths after migration.
- When temporary research is fully promoted, delete the promoted temporary copy.
- Git history is the archive for deleted stale documentation; do not create a
  generic legacy-document graveyard.

## Naming

Use stable responsibility-based names.

Prefer:

- `srs/README.md`, `architecture/system.md`, `contracts/README.md`,
  `operations/ci.md`;
- `src/<module>/README.md` for implementation ownership;
- `docs/research/<topic>/` for explorations.

Avoid catch-all names such as `NOTES.md`, `NEW.md`, or `FINAL.md` when the
content has an existing authoritative home.
