---
decided-by: [ADR-0022, ADR-0081]
level: [3]
areas:
  - requirements=requirement
  - architecture=architecture
  - design=design
  - reference=reference
  - decisions=decision
  - proposals=proposal
  - research=research
  - validation=validation
  - engineering=engineering
  - operations=operations
  - security=security
  - governance=governance
id-prefixes: [decision=ADR, proposal=PROP, research=RES]
source-roots: [src, crates/*/src]
forbidden-names: [archive, archives, deprecated, history, legacy, old, temp, temporary]
ledger: [governance/retired.md]
reachability-depth: [4]
review-days: [365]
max-lines: [500]
---

# DNA Documentation Policy

DNA adopts the repository-agnostic [[governance/documentation-architecture]]. The
front matter above is DNA's profile of that standard (Appendix A there); this
document records only what is DNA-specific. Do not restate the standard here.

## Areas in use

| Area | DNA-specific responsibility |
|---|---|
| [requirements](../requirements/README.md) | normative DNA SRS, quality attributes, and constraints |
| [architecture](../architecture/README.md) | system context, module boundaries, flow, interfaces, invariants |
| [design](../design/README.md) | current scientific and algorithmic mechanisms |
| [decisions](../decisions/README.md) | durable ADR history |
| [proposals](../proposals/README.md) | changes under review, plus the non-normative roadmap |
| [research](../research/README.md) | exploratory scientific and engineering evidence |
| [validation](../validation/README.md) | acceptance, datasets, limitations, benchmarks, traceability |
| [engineering](../engineering/README.md) | development, testing, CI/CD, release, dependencies, code quality |
| [operations](../operations/README.md) | readiness, observability, investigation playbook |
| [security](../security/README.md) | local-CLI trust boundaries, threat model, supply chain |
| [reference](../reference/README.md) | JSON schemas, result and configuration semantics, coordinates, glossary |
| [governance](README.md) | this policy, versioning, data, and external-reference validation |

DNA creates no service-only areas (deployment topology, SLOs, incident history,
disaster recovery) while those boundaries do not exist.

## Identity and metadata

- ADRs are `ADR-NNNN`, proposals `PROP-NNNN`, research topics `RES-<topic>`.
  Requirement items are bold list items such as `SRS-ALN-001`; the part before the
  number (`SRS-ALN`) is the namespace. Identifiers are never reused.
- A decision's front matter carries `status`, `date`, `family`, `supersedes`,
  `related`; its header lines are a generated `lifecycle` block.
- A current document declares the decisions and requirement namespaces it rests on
  in `decided-by` and `implements`, and the source it describes in `covers`. These
  are the edges `trace` follows; prose citations without a declaration fail
  `DOC-LNK-006`.
- Cross-references to documents and IDs use `[[wikilinks]]`; relative Markdown
  links stay for source files, schemas, router READMEs, and generated blocks.

## Source-local documentation

Every directory under `src/` and `crates/*/src/` has a `README.md` describing its
responsibility boundary. File-only Rust modules use rustdoc. DNA keeps no
`docs/src/` mirror.

## Commands and CI

All commands run from `tools/python`.

| Command | Purpose |
|---|---|
| `uv run python scripts/validate_docs_structure.py` | `DOC-STR-*`, `DOC-LNK-001`: folder READMEs, source READMEs, forbidden paths, Markdown links |
| `uv run python scripts/docs_graph.py check` | `DOC-ID-*`, `DOC-LNK-002…006`, `DOC-GRF-*`, `DOC-GEN-*`, `DOC-RET-*`, `DOC-FRS-*` |
| `uv run python scripts/docs_graph.py sync` | write missing edge declarations and refresh generated blocks |
| `uv run python scripts/docs_graph.py trace <path\|ID>` | trace closure of a source path, requirement item, decision, or document |
| `uv run python scripts/docs_graph.py impact <base>` | `DOC-IMP-001` against a base ref; CI runs it on pull requests |
| `uv run python scripts/docs_graph.py candidates` | retirement candidates |
| `uv run python scripts/docs_graph.py retire <doc> --outcome "…"` | move a document into the [ledger](retired.md) |

CI runs `validate_docs_structure.py`, `docs_graph.py check`, `docs_graph.py sync
--check`, and, on pull requests, `docs_graph.py impact`. To waive an impact finding,
write `docs-unchanged: <document> — <reason>` in the pull-request description.

Structural validation cannot prove prose accuracy. Reviewers still own semantic
correctness.

## ADR policy

ADRs are append-only while live and are retired through the ledger, not deleted
(decision: [[ADR-0081]]).

- Search existing documents and decisions before creating one; do not record
  ordinary implementation progress.
- Refine a decision without replacing its core choice by amending the record; a
  material change creates a successor whose `supersedes` names the old record (with
  a section qualifier for a partial supersession, for example `ADR-0058 §4`).
- Do not leave two accepted ADRs claiming one decision scope; the generated family
  map in the [ADR index](../decisions/adr/README.md) shows each concern's chain
  and the current documents it governs.
- Retire a fully superseded record after folding anything still true into its
  successor (standard section 8).

## Cleanup

Do not keep compatibility pointers, old-and-current parallel documents, or a
`legacy/` folder. Retired knowledge lives in the ledger and in git.

Versioning of contracts is in [[governance/versioning|versioning]]; ownership of areas is the
table above and the CODEOWNERS file when present.
