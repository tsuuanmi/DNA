# ADR-0022: Govern documentation as an executable knowledge system

## Status

Accepted

## Context

DNA is becoming a scientific software system with several kinds of documentation:

- normative requirements;
- architecture and design decisions;
- scientific/method descriptions;
- public schemas and configuration contracts;
- implementation manuals;
- validation evidence;
- research notes;
- roadmap material;
- agent instructions.

These documents do not have equal authority.

Without an explicit documentation model, humans and coding agents can select the wrong source of truth, treat research as production behavior, copy stale implementation notes, or make a locally reasonable refactor that violates a cross-cutting scientific invariant.

The problem is therefore not "write more documentation". The problem is to make documentation navigable, authoritative, traceable, and hard to misinterpret.

## Decision

DNA documentation will be governed as a layered knowledge system.

### 1. Root production documentation is authoritative by role

The production documentation set has distinct responsibilities:

- **SRS** defines what the current intended system MUST/SHOULD/MAY do.
- **Accepted ADRs** record why architectural or scientific decisions were made.
- **Architecture** defines component boundaries, dependencies, data flow, and cross-cutting invariants.
- **Method documentation** defines the current scientific and algorithmic behavior.
- **Contracts** define user-visible and machine-visible interfaces such as CLI, configuration, coordinates, schemas, and serialization semantics.
- **docs/src** mirrors implementation ownership and module responsibilities.
- **Validation** defines how software and scientific claims are verified.
- **Operations/governance** define development, release, security, data, and documentation processes.
- **Roadmap** describes future direction and is non-normative.

### 2. Research is explicitly non-normative

All exploratory work belongs under:

```text
docs/research/<topic>/
```

A research subtree may mirror production documentation with its own requirements, ADRs, architecture, validation, and roadmap, as done by the Tracy research work.

Research material does not change production behavior merely by existing.

Promotion requires an explicit path:

```text
research evidence
  ↓
root ADR when a decision is architectural/scientific
  ↓
root SRS/contract change when behavior changes
  ↓
implementation + tests
  ↓
validation
```

### 3. Documentation authority is explicit

When documents disagree:

1. a machine-readable public schema is authoritative for the serialized shape of that schema version;
2. root SRS is authoritative for intended normative behavior;
3. accepted ADRs are authoritative for the rationale and decision they govern;
4. current method/contracts define the intended detailed semantics;
5. source code is authoritative for what the current revision actually executes;
6. docs/src describes implementation ownership and must track source;
7. roadmap and research are never authority for current production behavior.

A disagreement between source and the normative production docs is a defect, incomplete implementation, or stale documentation. It must be surfaced rather than silently resolved by choosing whichever file is convenient.

### 4. Cross-cutting invariants have one explicit home

Stable invariants that span modules belong in `docs/architecture/invariants.md`.

Examples include:

- coordinate bases and units;
- source-evidence immutability;
- strand mapping;
- output atomicity;
- unresolved biological state semantics;
- deterministic ordering.

Implementation manuals should reference these invariants rather than restating them inconsistently.

### 5. Traceability connects intent to evidence

`docs/traceability.md` maps requirement families to:

- architecture/method documentation;
- owning source modules;
- tests;
- public contracts;
- validation evidence.

The traceability map is a navigation aid, not a duplicated specification.

### 6. AGENTS.md is a reusable workflow protocol, not a repository encyclopedia

Root `AGENTS.md` should define a generic development workflow that can transfer across repositories.

It should teach an agent how to:

- discover repository instructions and documentation roles;
- resolve authority by role rather than by hard-coded filenames;
- understand intent and implementation before editing;
- plan the smallest coherent change;
- choose verification based on failure modes;
- reconcile source, contracts, tests, and documentation;
- review the final diff and report residual risk.

Repository-specific commands, scientific invariants, file paths, and product semantics remain in the repository's own docs, CI configuration, build metadata, and contracts. The agent discovers those sources rather than having them duplicated into AGENTS.md.

### 7. Documentation mirrors responsibility, not line-by-line code

`docs/src` mirrors the source/crate/module layout so agents can discover ownership.

A module manual should describe:

- responsibility;
- inputs/outputs;
- invariants;
- dependencies;
- failure modes;
- related SRS/ADR/contracts/tests.

It should not translate implementation line by line.

### 8. ADRs are deduplicated by decision boundary

One accepted ADR is the canonical rationale for one durable decision scope.

Before adding a production ADR, contributors and agents must search the SRS,
architecture/invariants, method and contract documentation, and the ADR index for
an existing decision that already owns the concern.

A new implementation step, schema version, validation observation, or documentation
clarification does not by itself require a new ADR. If the underlying choice is
unchanged, update the authoritative SRS/method/contract and, when useful, add a
small dated clarification or follow-up to the existing ADR.

When a material choice changes, the replacement decision must explicitly supersede
the prior ADR in whole or in part, and the prior ADR must point to the replacement.
Two accepted ADRs must not claim authority over the same decision scope.

Superseded ADRs remain in the repository as historical provenance; their numbers
are never reused.

## Consequences

### Positive

- agents can navigate the repository without guessing which document is authoritative;
- research can be rich without contaminating current production semantics;
- cross-cutting invariants are less likely to drift between modules;
- source changes have a clear documentation impact path;
- the repository becomes easier to audit and hand over.

### Cost

- documentation changes require discipline when behavior crosses multiple layers;
- stale links or traceability entries become maintenance issues;
- some existing flat documents need classification or gradual relocation.

## Non-goal

This ADR does not require moving every existing Markdown file immediately. Structure should improve incrementally without creating churn that provides no additional clarity.
