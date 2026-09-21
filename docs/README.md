# DNA Documentation

This is the documentation entry point for humans and coding agents.

DNA documentation is organized by **authority and responsibility**. Do not infer current production behavior from a research note or roadmap item.

## Read order for a change

1. [SRS](srs/README.md) — what the system is intended to do.
2. [Architecture](architecture/README.md) and [system invariants](architecture/invariants.md) — where behavior belongs and what must always remain true.
3. [ADRs](adr/README.md) — why durable architectural or scientific choices were made.
4. [Methods](methods/README.md) — current scientific and algorithmic behavior.
5. [Contracts](contracts/README.md) — public and machine-visible interfaces.
6. [Implementation mirror](architecture/source-mirror.md) and affected source — implementation ownership and executable behavior.
7. [Traceability](traceability.md) and [operations](operations/README.md) — tests and evidence that protect the behavior.

## Authority

| Documentation | Role |
|---|---|
| [SRS](srs/README.md) | normative intended behavior |
| [Contracts](contracts/README.md) | exact public/configuration/serialization contract |
| [ADRs](adr/README.md) | durable decision and rationale |
| [Architecture](architecture/README.md) | boundaries and cross-cutting invariants |
| [Methods](methods/README.md) | detailed current scientific/algorithmic semantics |
| source code | actual behavior executed by the current revision |
| [docs/src mirror policy](architecture/source-mirror.md) | descriptive implementation ownership; `docs/src/` must track source one-to-one |
| [Governance](governance/README.md) | change, compatibility, data, and documentation policy |
| [Operations](operations/README.md) | development, batch, CI, release, and delivery procedures |
| [Roadmap](roadmap.md) | future direction; non-normative |
| [Research](research/README.md) | exploratory work; non-normative |

If source and normative production documentation disagree, surface the mismatch. Do not silently choose whichever artifact is convenient. See [documentation governance](governance/documentation.md).

## Directory layout

```text
docs/
├── README.md
├── srs/            # normative system requirements
├── architecture/   # system structure and invariants
├── adr/            # durable design/science decisions
├── methods/        # current algorithms and scientific behavior
├── contracts/      # config, result semantics, schemas, examples
├── src/            # exact source-module implementation mirror
├── operations/     # development, batch, CI, release, delivery
├── governance/     # documentation, compatibility, data policy
├── research/       # non-normative exploration
├── glossary.md
├── roadmap.md
└── traceability.md
```

Root-level documents are intentionally limited to cross-cutting entry points.

## Key entry points

- [SRS](srs/README.md)
- [Architecture](architecture/README.md)
- [ADR index](adr/README.md)
- [Methods](methods/README.md)
- [Contracts](contracts/README.md)
- [Implementation mirror](architecture/source-mirror.md)
- [Operations](operations/README.md)
- [Governance](governance/README.md)
- [Research](research/README.md)
- [Traceability](traceability.md)
- [Glossary](glossary.md)
- [Roadmap](roadmap.md)
