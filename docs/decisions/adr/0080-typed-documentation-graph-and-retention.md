---
id: ADR-0080
status: accepted
date: 2026-10-19
family: [documentation]
supersedes: [ADR-0022 §4/§6/§8]
related: [ADR-0006]
---

# ADR-0080: Document as a typed graph with generated indexes and ledgered retention

<!-- docs:begin lifecycle -->
- **Status:** Accepted
- **Date:** 2026-10-19
- **Supersedes:** [ADR-0022](0022-documentation-knowledge-system.md) (§4/§6/§8)
- **Related:** [ADR-0006](../../governance/retired.md)
<!-- docs:end lifecycle -->

<!-- docs:begin backlinks -->
- **Governs:** [governance/documentation](../../governance/documentation.md)
<!-- docs:end backlinks -->

## Context

[[ADR-0022]] made documentation an executable knowledge system: roles, README
routers, one canonical home per fact, and CI checks for structure and Markdown
links. Growth exposed what that did not cover.

- **Links are untyped.** Relations between requirements, designs, decisions, and
  source live in prose and in one hand-written traceability table. 46 of 58 ADRs
  were cited by no requirement, design, architecture, or reference document, and a
  document could depend on a decision without saying so. The table named a source
  file that no longer exists.
- **Derived tables are hand-maintained.** The ADR status table repeats each
  record's own status line; the decision-family map and the proposal lists are
  written by hand and drift.
- **Retention contradicts growth.** ADRs are append-only and "git is the archive",
  so superseded records stay in the tree forever, while deleting one leaves
  dangling links and no way to know which commit holds it. The documentation
  standard was an essay, not something another repository could adopt.

## Decision

1. **A repository-agnostic standard.**
   [[governance/documentation-architecture]] is a single self-contained
   specification with numbered rules, per-type front matter, conformance levels,
   an agent-contract template, and a profile template. DNA's profile is the front
   matter of [[governance/documentation]].
2. **Typed edges, authored once.** `implements`, `decided-by`, `covers`,
   `verified-by`, `supersedes`, `evidence`, and `related` are declared in front
   matter in one direction; inverses are computed. Prose citations of decisions and
   requirement items in current documents must be declared (`DOC-LNK-006`).
3. **`[[wikilinks]]`** are the cross-reference form for documents and IDs, and
   resolve through the same graph. Relative Markdown links stay for source files,
   routers, and generated text.
4. **Generated blocks.** Status indexes, the decision-family map, the proposal and
   research lists, ADR headers, backlinks, and the requirement-coverage table are
   produced by `docs_graph.py sync` and verified by `sync --check` in CI.
5. **Ledgered retention.** A fully superseded, rejected, implemented, or concluded
   document may leave the tree once its surviving content is folded into current
   documents. One ledger row keeps its ID, outcome, successor, and the commit and
   path of its full text. Retired IDs keep resolving.
6. **Impact check.** A pull request that changes a path covered by a document must
   change that document or state why not (`DOC-IMP-001`).

This replaces the parts of [[ADR-0022]] that said superseded knowledge is simply
removed with git as the only archive (§4), that `traceability.md` is a hand-written
map (§6), and that CI validates only structure and links (§8). Roles, routers,
colocated source READMEs, and one canonical home per fact stand.

## Alternatives considered

- **Git history only.** Rejected: nothing tells a reader which commit to search, and
  deletion breaks links.
- **An `archive/` or `legacy/` folder.** Rejected: it keeps the volume and invites
  readers and agents to treat it as current.
- **A stub file per retired record.** Rejected: the stubs accumulate as files; one
  ledger table is searchable, diffable, and cannot be mistaken for a record.
- **Front matter on every document.** Rejected: it adds noise where there is no
  state or edge. It is required only where it carries information.
- **Wikilinks everywhere, converted in bulk.** Rejected: churn and worse GitHub
  rendering for no new information; existing Markdown links stay valid.

## Consequences

- An agent can print everything connected to a path, requirement, or decision with
  one command instead of relying on memory.
- Indexes and coverage tables cannot drift; changing a status or an edge in one
  place updates every view.
- Decisions that govern no current document, documents past their size budget, and
  stale reviews surface as warnings, which makes pruning routine rather than
  occasional.
- Authors maintain front matter edges and run `sync`; the cost is small compared
  with hand-maintained tables, and the check tells them when they forgot.
- Wikilinks read as text on GitHub; Obsidian and the checker resolve them.
