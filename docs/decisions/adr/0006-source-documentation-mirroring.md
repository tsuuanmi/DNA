# ADR-0006: Mirror Rust Source Documentation One-to-One

- **Status:** Superseded by ADR-0022
- **Date:** 2026-08-22
- **Superseded:** 2026-09-22

## Context

A production scientific tool needs module documentation that stays aligned with
implementation. A separate prose hierarchy easily accumulates missing, orphaned,
or stale module pages.

## Options

1. Rely only on rustdoc comments.
2. Maintain hand-written docs without structural enforcement.
3. Mirror every `src/**/*.rs` file under `docs/src/` and check it in CI.

## Decision

Choose option 3, while retaining useful rustdoc. `src/x.rs` maps to
`docs/src/x.md`, and `src/x/mod.rs` maps to `docs/src/x/mod.md`. CI rejects a
missing or orphaned counterpart.

Each manual page records purpose, boundaries, inputs/outputs, invariants, errors,
Apollo mapping, SRS/ADR links, tests, and status. Semantic changes update source
and its manual page together.

## Consequences

Navigation and review ownership are explicit, and source-layout drift is caught
automatically. The repository carries more documentation files, and structural
checks cannot prove prose accuracy; review discipline remains necessary.

## Supersession

ADR-0022 replaces this one-to-one shadow tree with source-directory `README.md`
routers colocated with the code. The replacement preserves module ownership and CI
coverage while avoiding duplicated per-file documentation.
