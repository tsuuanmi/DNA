# AGENTS.md

This file is the repository-level routing and invariant guide for coding agents.
It follows the agent contract of the
[documentation architecture standard](docs/governance/documentation-architecture.md)
(Appendix B). Product knowledge belongs in canonical documentation and
source-local README files, not here.

## Read before changing

For any non-trivial change:

1. Read the repository `README.md` for scope and development entry points.
2. Read `docs/README.md` to find the authoritative requirements, architecture,
   decisions, design, reference, validation, and operations documentation.
3. **Trace every file, requirement, or decision you will change** and read the
   closure it prints; do not rely on memory of what else is connected. From
   `tools/python`:
   `uv run python scripts/docs_graph.py trace <path | SRS-… | ADR-… | document>`.
4. Read the nearest `README.md` in the affected source directory to understand
   module responsibility and dependency boundaries.
5. Read the affected source and tests.
6. Read research or roadmap material only for context; do not treat it as current
   production authority.

Resolve documentation by role rather than filename. If source and normative
documentation disagree, surface and reconcile the mismatch instead of silently
choosing one.

## Change invariants

- Make the smallest coherent production-quality change.
- Keep one authoritative implementation; remove obsolete, duplicate, fallback,
  legacy, and compatibility paths unless the current specification requires them.
- Preserve clear ownership, dependency direction, typed failures, deterministic
  behavior, bounded resource use, and domain invariants.
- Do not introduce speculative abstractions, unused configuration, or dependencies
  without a concrete requirement.
- Research is non-normative until explicitly promoted into current requirements,
  decisions, design, reference, implementation, tests, and validation.

Documentation is part of the change:

- code behavior change -> update affected canonical docs and tests;
- module responsibility/dependency change -> update the nearest source README;
- public API/schema/config/CLI change -> update reference/contracts, examples, requirements,
  and tests;
- algorithm/scientific behavior change -> update requirements, design, tests, and
  validation implications;
- architectural decision change -> update architecture and the owning ADR, or
  create a successor only when the decision materially changes;
- research-only work -> keep it under research until promotion.

Keep the knowledge graph in step with the change:

- declare relations in front matter (`implements`, `decided-by`, `covers`,
  `verified-by`, `supersedes`), not only in prose; refer to documents and IDs with
  `[[wikilinks]]`;
- run `docs_graph.py sync` after changing statuses or edges, and never edit a
  generated block by hand;
- drop a document with `docs_graph.py retire`, not `rm`, so its ID, outcome, and
  full-text commit stay in the ledger;
- when a pull request changes a path that a document `covers`, change that
  document or state `docs-unchanged: <document> — <reason>` in the description.

Do not duplicate the same fact across documents. Keep one canonical home and link
to it from routers and related material.

## Verification

Discover required checks from CI, build metadata, and operations documentation.
Run checks appropriate to the changed failure modes, then run the repository's
required gates before finalizing when possible.

Do not weaken, skip, or suppress a required gate merely to make a change pass.

## Final review

Before completion:

- inspect the complete diff;
- remove stale paths, dead documentation, and accidental compatibility layers;
- confirm source, tests, contracts, and documentation agree;
- confirm every source directory still has an up-to-date `README.md`;
- confirm `docs_graph.py check` and `docs_graph.py sync --check` pass;
- report what changed, what was verified, and any remaining uncertainty.

## Canonical navigation

Use [docs/README.md](docs/README.md) as the knowledge router. The reusable taxonomy,
link model, and retention rules are defined by the [documentation architecture
standard](docs/governance/documentation-architecture.md); DNA's profile, commands,
and CI enforcement are in the [documentation policy](docs/governance/documentation.md).
Development/verification commands belong in [engineering](docs/engineering/README.md),
not in this file.
