# Documentation Governance

DNA documentation is organized by authority, not by file age or document length.

## Authority model

| Layer | Purpose | Normative? |
|---|---|---|
| SRS (`srs/`) | Intended system requirements | Yes |
| JSON Schema / configuration contract | Exact machine-visible contract for the named version | Yes |
| Accepted ADR | Decision and rationale | Yes for the decision it governs |
| Method docs | Detailed current scientific/algorithmic semantics | Yes |
| Architecture / invariants | Module boundaries and cross-cutting truths | Yes |
| Source | Actual behavior of the current revision | Executable reality |
| docs/src mirror | Module ownership and implementation manual | Descriptive; must track source |
| Validation docs | Required evidence and acceptance method | Yes for validation policy |
| Research | Exploration and candidate designs | No |
| Roadmap | Future direction and priorities | No |

If source and normative documentation disagree, do not silently choose one. Surface the mismatch as a defect, incomplete implementation, or stale documentation and resolve it in the change that owns the behavior.

## Promotion path

Research becomes production behavior only through an explicit promotion path:

```text
docs/research/<topic>/
        ↓
promote durable decisions into docs/adr/ when needed
        ↓
promote current truth into the appropriate persistent authority
        ↓
source + docs/src + tests
        ↓
validation and release evidence
```

A research note is not production authority. Once exploratory content becomes current production truth, move that truth into its persistent authoritative document and remove the temporary duplicate.

## Change impact

A code change should update only the documentation layers it actually affects.

- Internal refactor with unchanged behavior: update `docs/src` only when ownership/responsibility changes.
- Scientific behavior change: update SRS, method docs, relevant ADR if needed, tests, validation implications, and `docs/src`.
- Schema/config/CLI change: update the machine contract, human contract, SRS, tests and examples.
- New architectural dependency or boundary: update architecture and usually an ADR.
- Research-only work: keep it under `docs/research/<topic>/`; do not edit root production contracts until promotion.

## ADR lifecycle and deduplication

Treat ADRs as a decision log, not a changelog.

- Search the SRS, architecture/invariants, methods/contracts, and ADR index before creating a new ADR.
- Reuse the existing ADR when the decision boundary and core choice are unchanged; implementation progress, schema revisions, validation results, and clarifications normally belong in their owning documentation.
- Create a replacement ADR only for a material new choice, then mark the previous ADR `Superseded` or `Superseded in part` and link the relationship both ways.
- Do not keep two `Accepted` ADRs that claim authority over the same decision scope.
- Keep superseded ADRs for provenance and never reuse their numbers.

## Staleness rules

- A `docs/src` file without a matching source file is stale.
- A source module without the required mirror is undocumented.
- An ADR marked Superseded must point to the replacing decision.
- Two Accepted ADRs must not claim authority over the same decision scope.
- Examples must validate against their named schema.
- Roadmap or research text must not be used to justify current production behavior.

## Folder indexes and promotion cleanup

- Every documentation folder under `docs/` MUST contain a canonical `README.md` index, except `docs/src/` and its descendants because that tree must remain an exact one-to-one mirror of Rust source files.
- A folder index defines the purpose, authority, and navigation for that folder. Do not keep a sibling `<folder>.md` as a second index.
- When exploratory/temporary documentation is fully promoted into persistent SRS, architecture, ADR, method, contract, governance, or operations documentation, delete the promoted temporary document or folder.
- Do not keep compatibility pointers, legacy copies, renamed duplicates, or “old/new” parallel documentation paths after migration.
- Retain a research document only while it still contains active non-normative research that has not been promoted.

## Naming

Use stable role-based names over temporary project names.

Prefer:

- `srs/README.md`, `architecture/system.md`, `contracts/README.md`, `operations/ci.md`, `roadmap.md`;
- `docs/research/<topic>/` for explorations;
- `docs/src/<same-relative-path>.md` for implementation manuals.

Avoid new catch-all files such as `NOTES.md`, `NEW.md`, or `FINAL.md` when the content has an existing authoritative home.
