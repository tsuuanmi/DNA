# ADR-0006: Mirror Rust Source Documentation One-to-One

- **Status:** Accepted
- **Date:** 2026-08-22

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

Any replacement must provide equivalent automatic coverage and traceability.
