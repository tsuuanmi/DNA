# ADR-0002: Use a Single Layered Crate for the MVP

- **Status:** Accepted
- **Date:** 2026-08-22

## Context

DNA needs clear boundaries between binary input, scientific algorithms,
orchestration, and output. A premature multi-crate workspace would add versioning
and dependency overhead before any external library consumer exists.

## Options

1. One monolithic command module.
2. One library-plus-binary crate with internal modules.
3. A workspace with separate model, format, algorithm, and CLI crates.

## Decision

Choose option 2. `main` is a thin process boundary; `pipeline` orchestrates;
scientific modules operate independently of CLI and output side effects; `model`,
`config`, and `error` are shared inward dependencies.

## Consequences

The design has production boundaries without workspace overhead. Internal module
APIs can evolve during MVP work. A later crate split requires evidence of reuse,
build isolation, or independent release needs.

## Supersession

A workspace ADR may supersede this decision after identifying concrete consumers
and a cycle-free dependency graph.
