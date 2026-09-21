# DNA Documentation

This is the documentation entry point for humans and coding agents.

DNA documentation is organized by **authority and responsibility**. Do not infer current production behavior from a research note or roadmap item.

## Read order for a change

1. [SRS](srs/README.md) — what the system is intended to do.
2. [Architecture](architecture/README.md) and [system invariants](architecture/invariants.md) — where the behavior belongs and what must always remain true.
3. [Accepted/proposed ADRs](adr/README.md) — why relevant decisions were made.
4. [Methods](methods/README.md) — current scientific/algorithmic behavior.
5. [Contracts](contracts/README.md) — public/machine-visible interfaces and coordinate semantics.
6. [Source mirror](source-mirror.md) and affected source — implementation ownership and current executable behavior.
7. [Traceability](traceability.md) and CI/release operations — tests and evidence that protect the behavior.

## Authority

| Documentation | Role |
|---|---|
| [SRS](srs/README.md) | normative intended behavior |
| JSON schemas / config contract | exact machine-visible contract for the named version |
| [ADRs](adr/README.md) | decision and rationale |
| [Architecture + invariants](architecture/README.md) | boundaries and cross-cutting truths |
| [Methods](methods/README.md) | detailed current scientific/algorithmic semantics |
| source code | actual behavior executed by the current revision |
| [docs/src](source-mirror.md) | descriptive module ownership; must track source |
| [Roadmap](roadmap.md) | future direction; non-normative |
| [Research](research/README.md) | exploratory work; non-normative |

If source and normative production documentation disagree, surface the mismatch. Do not silently choose whichever artifact is convenient. See [documentation governance](governance/documentation.md).

## Product and requirements

- [SRS](srs/README.md)
- [Roadmap](roadmap.md)
- [Development readiness](development-readiness.md)
- [Glossary](glossary.md)

## Architecture and decisions

- [Architecture index](architecture/README.md)
- [System overview](architecture.md)
- [System invariants](architecture/invariants.md)
- [Source layout](source-layout.md)
- [ADR index](adr/README.md)
- [Traceability](traceability.md)

## Current methods

- [Method index](methods/README.md)
- [Pipeline](pipeline.md)
- [DNA processing](signal-processing.md)

## Contracts

- [Contract index](contracts/README.md)
- [Coordinate contract](contracts/coordinates.md)
- [Configuration](configuration.md)
- [Reference-free basecall output](basecall-output.md)
- [Sample evidence output](sample-output.md)
- [Reference analysis output](json-output.md)
- [Schemas](schemas/)
- [Examples](examples/)

## Implementation manuals

- [Source mirror](source-mirror.md) — mirrors `src/**/*.rs`.

Every mirrored source file has a same-relative-path manual. Manuals describe ownership, inputs/outputs, invariants, dependencies, failure modes, and traceability rather than translating code line by line.

## Validation and operations

- [Development/release operations](operations/README.md)
- [CI lanes](operations/ci.md)
- [Release operations](operations/release.md)
- [Security and trust boundaries](operations/security.md)
- [Data policy](data.md)
- [Delivery record](delivery-record.md)

## Governance

- [Governance index](governance/README.md)
- [Documentation governance](governance/documentation.md)
- [Compatibility](compatibility.md)
- [Changelog](../CHANGELOG.md)

## Research

Research lives only under [`docs/research/`](research/README.md) and does not change production behavior until promoted through the root ADR/SRS/contract process.

- [DNA research](research/DNA/README.md)
- Tracy research follows the same model under `docs/research/Tracy/` once that research subtree is integrated.

The old catch-all DNA research files were moved into `docs/research/DNA/` so they cannot be mistaken for production requirements.
