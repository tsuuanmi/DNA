# ADR-0003: Target Behavioral Compatibility with Apollo C++

- **Status:** Superseded in part by ADR-0009
- **Date:** 2026-08-22

## Context

Apollo documentation, C++ source, and the previous Rust port disagree in several
places. Exact CLI and serialized-byte compatibility would also preserve legacy
interface constraints that are unnecessary for DNA.

## Options

1. Match the previous Rust port.
2. Match Apollo CLI and files byte-for-byte.
3. Treat C++ headers as the scientific oracle while versioning DNA interfaces.

## Decision

Choose option 3. Exact deterministic stage results are preferred. Equivalent
indel placement is compared after normalization. Any intentional scientific
divergence requires its own ADR and validation evidence.

## Consequences

DNA can use clean typed interfaces and a stable schema without inheriting known
Rust-port shortcuts. Compatibility work requires source mapping and carefully
provenanced fixtures; existing loose goldens are not sufficient by themselves.

## Supersession

A future release may add a separate Apollo compatibility adapter without changing
the scientific core.
