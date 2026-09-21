# Documentation Governance

DNA documentation is organized by authority, not by file age or document length.

## Authority model

| Layer | Purpose | Normative? |
|---|---|---|
| SRS | Intended system requirements | Yes |
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
docs/research/<topic>
        ↓
root ADR if a decision is architectural/scientific
        ↓
root SRS and contract change when behavior changes
        ↓
source + docs/src + tests
        ↓
validation evidence
        ↓
changelog/release evidence when user-visible
```

A research ADR is not a production ADR.

## Change impact

A code change should update only the documentation layers it actually affects.

- Internal refactor with unchanged behavior: update `docs/src` only when ownership/responsibility changes.
- Scientific behavior change: update SRS, method docs, relevant ADR if needed, tests, validation implications, and `docs/src`.
- Schema/config/CLI change: update the machine contract, human contract, SRS, tests, examples, and changelog.
- New architectural dependency or boundary: update architecture and usually an ADR.
- Research-only work: keep it under `docs/research/<topic>/`; do not edit root production contracts until promotion.

## Staleness rules

- A `docs/src` file without a matching source file is stale.
- A source module without the required mirror is undocumented.
- An ADR marked Superseded must point to the replacing decision.
- Examples must validate against their named schema.
- Roadmap or research text must not be used to justify current production behavior.

## Naming

Use stable role-based names over temporary project names.

Prefer:

- `requirements.md`, `architecture.md`, `operations/ci.md`, `roadmap.md`;
- `docs/research/<topic>/` for explorations;
- `docs/src/<same-relative-path>.md` for implementation manuals.

Avoid new catch-all files such as `NOTES.md`, `NEW.md`, or `FINAL.md` when the content has an existing authoritative home.
