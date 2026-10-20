# Documentation Architecture Standard

A repository-agnostic standard for a documentation **knowledge system**: typed
documents, a checked link graph, generated indexes, and managed retention.

**To adopt it in another repository, give an agent this file.** Nothing else is
required. Section 12 says what to create, Appendix A is the repository profile to
fill in, Appendix B is the agent contract, Appendix C is the check catalogue, and
Appendix D is the generator catalogue. No product, language, or toolchain is
assumed.

Key words MUST, SHOULD, and MAY are used as in RFC 2119. Each normative rule has
an ID (`DOC-…`) so a checker, a review comment, and a waiver can cite it.

## 1. Purpose

The system lets a reader or coding agent answer, without guessing:

1. What must be true? 2. Why was a durable choice made? 3. What change is
proposed? 4. What has been explored? 5. How is the system structured? 6. How does
the mechanism work? 7. What is exact? 8. Where is it implemented? 9. How is it
verified? 10. How is it built, released, operated, secured, and evolved?

and, for any file or requirement, **everything that must change with it**.

## 2. Principles

- **One fact, one home.** Other documents link to a fact; they do not restate it.
- **Documentation changes with the system**, in the same change set.
- **Authority is explicit.** Current truth, history, proposed change, and
  exploration are different classes (section 3); a lower-authority artifact never
  silently becomes production truth.
- **Edges are authored once and checked.** Relations between documents are
  declared in one direction; the inverse is computed (section 6).
- **Whatever can be derived is generated.** Indexes, status tables, coverage
  tables, and backlinks are produced by a tool and verified current by CI. A
  human never maintains a derived table by hand (section 7).
- **Documentation grows, then is pruned with evidence.** Retired knowledge leaves
  the tree but stays findable (section 8).

## 3. Knowledge classes and areas

| Class | Meaning | Retention |
|---|---|---|
| canonical / living | what is true now | update in place; retire when obsolete |
| historical / durable | why a decision was made, what happened | append-only while live; retire per section 8 |
| evolutionary | a change being proposed or implemented | retire after promotion |
| exploratory | evidence, experiments, ideas | promote, or retire with a summary |
| executable reality | source and tests | not documentation; linked from it |

An **area** is a top-level folder under `docs/`. Create an area only when the
repository has real knowledge for it; do not invent service-only areas (SLOs,
disaster recovery, incident history) a repository does not have.

| Area (`type`) | Answers | Class |
|---|---|---|
| `requirements/` (`requirement`) | what must be true | canonical |
| `architecture/` (`architecture`) | where responsibilities and invariants belong | canonical |
| `design/` (`design`) | how current mechanisms work | canonical |
| `reference/` (`reference`) | exact interfaces, schemas, formats, glossary | canonical |
| `validation/` (`validation`) | how claims are demonstrated; limitations | canonical |
| `engineering/` (`engineering`) | build, test, CI/CD, dependencies, release | canonical |
| `operations/` (`operations`) | readiness, observability, runbooks, playbooks | canonical |
| `security/` (`security`) | trust boundaries, threats, supply chain | canonical |
| `governance/` (`governance`) | this system's own rules, ownership, versioning | canonical |
| `decisions/` (`decision`) | why a durable choice was made | historical |
| `proposals/` (`proposal`) | a change under review | evolutionary |
| `research/` (`research`) | what has been learned, not yet accepted | exploratory |

Boundaries: a proposal is not a specification; a decision record is not a design
guide; research is not production authority. Runbook = known procedure;
playbook = guided investigation of a known symptom; postmortem = learning from a
real incident.

**Implementation documentation lives with implementation.** Do not keep a
`docs/src/` shadow tree. Each real source package directory has a `README.md`
stating what it owns, what it does not, its entry points, its dependency
direction, and the `covers`-linked documents that govern it. File-local API
detail belongs in the language's own documentation system.

## 4. Document model

### 4.1 Identity

- Every document has a **path identity**: its path under `docs/` without `.md`;
  a folder `README.md` is identified by its folder (`design/README.md` → `design`).
- Documents of a type with an **ID prefix** (profile, Appendix A) also have an
  **ID** `PREFIX-NNNN` (for example `ADR-0012`) that equals their file-name
  prefix, is unique, and **is never reused**, even after retirement.
- A list item that begins with a bold ID (`- **SRS-ALN-001:** …`) **defines that
  item** wherever it appears; the item resolves to its defining document. The
  ID without its trailing number (`SRS-ALN`) is the item's **namespace**.

`DOC-ID-001` Every document has a type: its front matter `type`, or else the type
of its area.
`DOC-ID-002` IDs are unique and match their file-name prefix.
`DOC-ID-003` `status` is a value allowed for the type (4.3); when absent it is
`active` (living types); typed-ID types require `id` and `status`.

### 4.2 Front matter

Front matter is **required where it carries information** (an ID, a lifecycle
status, or a declared edge) and omitted otherwise. Do not add it mechanically.

Front matter is a strict YAML subset so any language can parse it without a
library: flat `key: value` lines; a value is a scalar or a list of scalars,
written inline (`[a, b]`) or as `- item` lines. No nesting, anchors, or
multi-line scalars. Dates are `YYYY-MM-DD`.

| Field | Used by | Meaning |
|---|---|---|
| `type`, `status` | all | see 4.1, 4.3 |
| `id` | typed-ID types | the ID |
| `date` / `created` | decision / proposal | creation date |
| `summary` | proposal, research | one sentence for generated indexes |
| `owners` | any | people or teams accountable |
| `reviewed` | living types | date the content was last confirmed current |
| `implements` | design, architecture, reference | requirement items or namespaces it satisfies |
| `decided-by` | any current doc | decisions whose rationale governs it |
| `covers` | any | source/test/config path globs it describes or governs |
| `verified-by` | requirement, validation | tests or validation documents that demonstrate it |
| `supersedes` | decision, proposal | IDs this record replaces (wholly or in part) |
| `evidence` | proposal, decision | research documents behind it |
| `family` | decision | decision-family labels (for generated family maps) |
| `related` | any | a deliberate non-typed link; use sparingly |
| `implementation` | proposal | external references (PRs, commits) that delivered it |
| `unverified` | requirement | items or namespaces knowingly without a verifier, waiving `DOC-GRF-004` |
| `family-labels` | a folder README | `slug=Label` pairs naming decision families for `families:decision` |

The relation fields (`implements`, `decided-by`, `verified-by`, `covers`,
`supersedes`, `evidence`, `related`) are the **closed relation vocabulary**.
Unknown fields are errors, so a typo cannot silently drop an edge.

### 4.3 Status vocabulary

| Types | Allowed `status` |
|---|---|
| living types | `draft`, `active`, `deprecated` |
| `decision` | `proposed`, `accepted`, `rejected`, `superseded`, `superseded-in-part`, `deprecated` |
| `proposal` | `draft`, `proposed`, `accepted`, `rejected`, `implementing`, `implemented`, `superseded`, `withdrawn` |
| `research` | `active`, `concluded`, `promoted` |

A `deprecated` document is a retirement candidate and MUST be retired (section 8)
or restored.

## 5. Links

### 5.1 Syntax

A cross-reference to another document or ID is written as an Obsidian-style
**wikilink**:

- `[[design/alignment]]` — path identity; `[[design/alignment#Gap placement]]`
  with a heading;
- `[[ADR-0012]]`, `[[PROP-0004]]` — a document ID;
- `[[SRS-ALN-003]]` — a requirement item (resolves to its defining document);
- `[[ADR-0012|the gap decision]]` — with display text.

A wikilink target is resolved by ID, then path identity, then unique base name.
An ambiguous or unresolved target is an error.

Relative Markdown links remain **valid and equivalent** for targets that are not
documents (source files, schemas, external URLs) and in router READMEs and
generated blocks, where rendering on the hosting site matters. A bare ID in prose
(`ADR-0012`) is also resolved and checked, but a wikilink is preferred.

`DOC-LNK-001` Relative Markdown links resolve to existing repository paths.
`DOC-LNK-002` Wikilinks resolve to exactly one live document or ledger entry.
`DOC-LNK-003` Bare IDs in prose resolve to a live document or ledger entry.
`DOC-LNK-004` Repository paths named in code spans exist (warning where a path is
illustrative).

### 5.2 Typed edges

Prose links carry no meaning to a machine. **Typed edges are declared in front
matter** (4.2) and authored **in one direction only** — the downstream document
declares its upstream:

| Authored in | Field | Points to | Computed inverse |
|---|---|---|---|
| design / architecture / reference | `implements` | requirement | `implemented-in` |
| any current doc | `decided-by` | decision | `governs` |
| any doc | `covers` | source/test/config globs | `documented-by` (for a path) |
| requirement / validation | `verified-by` | tests, validation docs | `verifies` |
| successor decision / proposal | `supersedes` | earlier record | `superseded-by` |
| proposal / decision | `evidence` | research | `informs` |

`DOC-LNK-005` Relation fields use only the vocabulary above and every target
resolves (to a document, an ID, or at least one existing path for `covers` and
`verified-by`).
`DOC-LNK-006` Every decision ID or requirement item a **living** document cites
in prose is also declared in its front matter (`decided-by`, `implements`; a
namespace entry covers all its items), so no edge lives only in prose. A tool MAY
write the missing declarations. Historical, evolutionary, and exploratory
documents cite freely and are exempt.

### 5.3 The trace closure

For a document, an ID, or a source path, the **trace closure** is everything
connected through typed edges in both directions: the requirement items it
implements, their other implementers, governing decisions and their supersession
chain, covering documents, verifiers, owning source README, and generated
backlinks. A tool MUST be able to print it (`trace`, Appendix D). It is how an
agent answers "what else must change".

`DOC-GRF-001` Supersession is symmetric: a record named by `supersedes` has status
`superseded` or `superseded-in-part`, and no such status exists without a
successor.
`DOC-GRF-002` Every document is reachable from `docs/README.md` through router
links within a bounded depth (profile; default 3).
`DOC-GRF-003` `covers` globs match at least one path.
`DOC-GRF-004` (warning) A requirement item is implemented by at least one
document and verified item by item: cited by a test, code, or enforcing
workflow (`SRS-X-001` in a comment), or listed in `verified-by`, or waived with
`unverified:` (item or namespace) with a stated reason (`SRS-X-001 process rule`). A `verified-by` path on a document is context
for readers and the trace closure; it does not verify its items.
`DOC-GRF-005` (warning) An accepted decision is governing at least one current
document; otherwise it is a retirement candidate.

## 6. README routers

A README orients; it does not duplicate. The repository README routes users and
contributors; `docs/README.md` routes by knowledge role and authority; a folder
README states what belongs, what does not, its canonical entry points, and any
lifecycle ambiguity; a source README routes from implementation ownership to
canonical knowledge. The list of a folder's members in a router SHOULD be a
generated block (section 7), not hand-maintained.

`DOC-STR-001` Every documentation folder has exactly one `README.md`.
`DOC-STR-002` Every source package directory named by the profile has a `README.md`.
`DOC-STR-003` No `docs/src/` mirror and no `topic.md` beside `topic/README.md`.
`DOC-STR-004` No file or folder named by the profile's forbidden names
(`legacy`, `archive`, `old`, `temp`, …) exists as a current-state path.

## 7. Generated content

A **generated block** is delimited by

```text
<!-- docs:begin NAME -->
…generated…
<!-- docs:end NAME -->
```

The tool rewrites only the text between markers, and generated text uses
Markdown links so it renders everywhere. Generator contracts are in Appendix D.
`sync` rewrites all blocks; `sync --check` fails when any block is stale. CI runs
`sync --check`.

`DOC-GEN-001` Every generated block is current.
`DOC-GEN-002` A table derivable from front matter (a status index, a family map,
a coverage table) is generated, never hand-written.

## 8. Retention: grow, then prune with evidence

Git history alone is a poor archive: nobody knows which commit to search, and a
deleted file leaves dangling links. Retention therefore has three states:

```text
Live (full text in the tree)
  → Tombstone (one ledger row; full text in a named commit)
  → Git (the commit)
```

### 8.1 The ledger

One document (profile `ledger`, default `governance/retired.md`) holds a table,
one row per retired document:

| id | title | outcome | superseded-by | retired | full-text |
|---|---|---|---|---|---|
| `ADR-0004` | Poll the upstream service | replaced by an event stream | `ADR-0010` | 2026-10-19 | `a1b2c3d:docs/decisions/adr/0004-….md` |

`outcome` is one sentence: what was decided or learned and why it stopped being
live. `full-text` is `<commit>:<path>` of the last commit containing the file.
The ledger is append-only. Create a git tag (`docs-retired/<date>`) at that
commit so history rewriting and garbage collection cannot lose it.

### 8.2 When a document may be retired

A document MAY be retired only if all hold:

1. it is `superseded` (not merely in part), `rejected`, `withdrawn`,
   `implemented` (proposal), `promoted` or `concluded` (research), or
   `deprecated`;
2. everything still true in it lives in a current document — **fold before
   retire**; a record superseded only in part stays live until its surviving
   part is folded into the successor or a design rationale;
3. no current document links to its path except through its supersession chain
   or the ledger — links are rewritten to the successor or to the ledger row
   (`[[ADR-0004]]` keeps working because the ID resolves to the ledger).

`DOC-RET-001` A retired ID resolves through the ledger; its path does not exist.
`DOC-RET-002` Ledger rows are well-formed and their `superseded-by` IDs resolve.
`DOC-RET-003` (warning) An ID absent from the tree and the ledger is a numbering
gap; record it as `never-committed` or retire it properly.
`DOC-RET-004` IDs are never reused.

### 8.3 Keeping growth bounded

- Retire on a cadence (each release) using the retirement candidates the tool
  lists: `deprecated`, fully superseded, implemented proposals, concluded
  research, and decisions governing nothing (`DOC-GRF-005`).
- Prefer **folding a chain** — several micro-decisions in one family — into one
  successor record over keeping every step live.
- Keep documents under the profile's size budget (`DOC-FRS-002`); split or
  distill instead of letting one file become the repository's encyclopedia.
- Mark living documents with `reviewed`; stale ones are surfaced
  (`DOC-FRS-001`) for confirmation or retirement.

`DOC-FRS-001` (warning) A living document's `reviewed` date is older than the
profile's review window.
`DOC-FRS-002` (warning) A living document exceeds the profile's maximum line count.
Evolutionary and historical records are exempt; they are retired or distilled
instead.

## 9. Change impact

The change that alters behavior owns the documentation update. The tool checks it:
given a diff, every document whose `covers` globs match a changed path MUST
appear in the diff or be waived with a reason in the change description
(`docs-unchanged: <path> — <reason>`).

`DOC-IMP-001` A changed path covered by a document leaves that document changed
or waived.

Typical impact: behavior → requirements, design, tests, validation; module
responsibility or dependency → source README, architecture; public contract →
reference, requirements, examples, tests; durable choice → decision record or a
successor; unaccepted change → proposal; exploration → research.

## 10. Promotion lifecycle

```text
research → proposal / review → decision (when durable rationale is needed)
  → requirements + architecture + design + reference
  → source + tests → validation → release → operations feedback
```

Not every change needs every stage. The rule is that a lower-authority artifact
cannot silently become production truth.

## 11. Conformance levels

| Level | Rules enforced as errors |
|---|---|
| 1 Minimum | `DOC-STR-*`, `DOC-ID-*`, `DOC-LNK-001` |
| 2 Standard | level 1 + `DOC-LNK-002…006`, `DOC-GRF-001…003`, `DOC-GEN-*` |
| 3 Enterprise | level 2 + `DOC-RET-*`, `DOC-IMP-001`; warnings `DOC-GRF-004/005`, `DOC-FRS-*` |

The profile states the level. A small repository SHOULD start at level 1.

## 12. Adopting in another repository

Given only this file, an agent:

1. Reads Appendix A and writes the repository profile (the front matter of
   `docs/governance/documentation.md`) with the repository's areas, ID prefixes,
   source roots, and level.
2. Creates `docs/README.md` and a `README.md` in each area that has content;
   creates only areas with real knowledge (section 3).
3. Writes `AGENTS.md` from Appendix B, filling the repository-specific lines.
4. Inventories existing documents, classifies each by area and class, chooses one
   home per fact, and moves current truth before deleting duplicates.
5. Adds IDs and front matter where required (4.2); migrates prose citations into
   typed edges (`DOC-LNK-006`).
6. Implements the checks of Appendix C and generators of Appendix D in the
   repository's own tooling language (or reuses an existing implementation), and
   adds `check`, `sync --check`, and the impact check to CI.
7. Replaces every hand-written index or status table with a generated block.
8. Starts the ledger (empty) so retirement is available from day one.

Do not keep compatibility pointers for old paths unless stable published URLs
are a stated requirement; the ledger and git are the history.

## Appendix A — Repository profile

The profile is the front matter of `docs/governance/documentation.md`:

```yaml
---
level: 2
areas: [requirements=requirement, architecture=architecture, design=design, reference=reference, decisions=decision, proposals=proposal, research=research, validation=validation, engineering=engineering, operations=operations, security=security, governance=governance]
id-prefixes: [decision=ADR, proposal=PROP, research=RES]
source-roots: [src]
forbidden-names: [archive, deprecated, history, legacy, old, temp, temporary]
ledger: governance/retired.md
reachability-depth: 3
review-days: 365
max-lines: 500
---
```

Below the front matter the same document records the repository-specific
application: the actual areas in use, the commands that run the checks, and local
conventions. It does not restate this standard.

## Appendix B — Agent contract (`AGENTS.md` template)

Agent instructions route and constrain; they are not an encyclopedia. Use:

````markdown
# AGENTS.md

Repository-level routing and invariants for coding agents. Product knowledge lives
in `docs/`, not here.

## Read before changing

1. `README.md` for scope and entry points.
2. `docs/README.md` to find requirements, architecture, decisions, design,
   reference, validation, and operations by role.
3. Run `<trace command> <path-or-ID>` for every file or requirement you will
   change, and read the closure it prints: requirements, design, decisions,
   reference, source README, tests, validation. Do not rely on memory.
4. The nearest source `README.md`, then the source and tests.
5. Research and proposals are context only, never current authority.

If source and normative documentation disagree, surface and reconcile; do not
silently pick one.

## Change invariants

- Smallest coherent production-quality change.
- One authoritative implementation; remove obsolete and duplicate paths.
- <repository-wide invariants: dependency direction, determinism, error typing,
  resource bounds, domain rules>.
- No speculative abstractions, unused configuration, or unrequested dependencies.

## Documentation is part of the change

| Change | Update |
|---|---|
| behavior | requirements, design, tests, validation |
| module responsibility/dependency | source README, architecture |
| public API/schema/config/CLI | reference, requirements, examples, tests |
| durable choice | decision record or successor |
| unaccepted idea | proposal; exploration → research |

Edit typed edges in front matter, never only in prose. Run `<sync command>`; never
hand-edit a generated block. To drop a document, use `<retire command>`, not
`rm`.

## Verification

Discover the gates from CI and the engineering docs. Run those suited to the
change, then the full required set. Do not weaken, skip, or suppress a gate to
pass.

## Final review

Inspect the whole diff; remove stale paths and accidental compatibility layers;
confirm source, tests, contracts, and docs agree; confirm `<check command>` and
`<sync command> --check` pass; report what changed, what was verified, and what
remains uncertain.
````

## Appendix C — Check catalogue

A conforming checker reports `error` or `warning` with the rule ID, the file, and
the line. Rules and levels are in sections 4–9 and 11. Checks that need no graph
(`DOC-STR-*`, `DOC-LNK-001`) SHOULD run even when front matter is malformed.
`DOC-IMP-001` additionally takes a base reference to diff against.

## Appendix D — Generator and command catalogue

| Block name | Content |
|---|---|
| `lifecycle` | a decision's header lines: status (with successors and section qualifiers), date, supersedes, related |
| `index:<type>[:<group>]` | every document of the type: ID/path, title, status, summary; retired members come from the ledger. Groups: proposals and research have `open` (draft, proposed, accepted, implementing / active) and `done` (the rest) |
| `families:decision` | per `family` named in the containing README's `family-labels`, the decisions (superseded ones marked) and the current documents that name them in `decided-by` |
| `coverage:requirement` | per requirement namespace: defining documents, item count, implementing documents, governing decisions, covered source, items cited in code or tests, items waived as unverified |
| `backlinks` | the computed inverse edges of the containing document: `Governs` (decisions) and `Implemented in` (requirements) |

| Command | Behavior |
|---|---|
| `check` | run the catalogue for the profile's level |
| `sync [--check]` | rewrite generated blocks and `DOC-LNK-006` declarations; `--check` only reports staleness |
| `trace <path\|ID>` | print the trace closure (5.3) |
| `impact <base-ref>` | `DOC-IMP-001` for the diff against `base-ref` |
| `retire <path> --outcome "…"` | verify 8.2, append the ledger row, delete the file, rewrite inbound links to the successor, print the tag command |
| `candidates` | list retirement candidates (8.3) |

## Appendix E — Influences

Small, current, README-routed documentation updated with code (Google); ADRs as
an append-only log with superseding records (Azure Well-Architected); explicit
proposal states separate from implementation (Swift Evolution); lifecycle
metadata and production-readiness review (Kubernetes KEPs); runbook versus
playbook (AWS Well-Architected); wikilink and property conventions
(Obsidian). Influences, not a claim that any one source uses this exact model.
