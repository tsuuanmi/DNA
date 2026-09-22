# AGENTS.md

This file is the repository-level routing and invariant guide for coding agents.
It is intentionally concise. Product knowledge belongs in canonical documentation
and source-local README files, not here.

## Read before changing

For any non-trivial change:

1. Read the repository `README.md` for scope and development entry points.
2. Read `docs/README.md` to find the authoritative requirements, architecture,
   decisions, design, reference, validation, and operations documentation.
3. Read the nearest `README.md` in the affected source directory to understand
   module responsibility and dependency boundaries.
4. Read the affected source and tests.
5. Read research or roadmap material only for context; do not treat it as current
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
- report what changed, what was verified, and any remaining uncertainty.

## Canonical navigation

Use [docs/README.md](docs/README.md) as the knowledge router. Development/verification commands belong in [engineering](docs/engineering/README.md), not in this file.
