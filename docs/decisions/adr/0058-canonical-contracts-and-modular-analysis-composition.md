# ADR-0058: Canonical contracts and modular analysis composition

- **Status:** Superseded in part by [ADR-0069](0069-plugin-first-modality-core-post-calling.md) (§4) and [ADR-0070](0070-workspace-split-by-plugin-family.md) (§9)
- **Date:** 2026-10-03
- **Related proposal:** [PROP-0001](../../proposals/0001-modular-dna-analysis-platform.md)

## Context

DNA currently has one production-oriented Sanger path, but the intended system
must support independently evolving input modalities, scientific algorithms,
downstream analysis capabilities, and scientific data sources.

Examples include Sanger and NGS input, alternative alignment methods, alternative
variant callers, haplogroup classification against different databases, and
nomenclature against different rule or knowledge sources.

If adjacent components exchange implementation-specific representations, a
change in one implementation forces changes through the rest of the pipeline.
If all variation is hidden behind one generic plugin interface, scientifically
different operations lose their domain-specific contracts and validation
requirements.

DNA therefore needs stable semantic boundaries without prematurely turning every
module into a trait, crate, or dynamically loaded plugin.

ADR-0002 remains authoritative for the current single-crate layout: a workspace
split still requires concrete evidence of reuse, build isolation, ownership, or
independent release needs.

## Decision

DNA will evolve as a modular analysis platform centered on canonical typed
contracts.

### 1. Canonical contracts define component boundaries

A replaceable component depends on the semantic input and output contract it
needs, not on the concrete implementation before or after it.

Canonical contracts own shared meaning such as coordinate conventions, strand,
reference identity, unresolved states, result ordering where relevant, required
evidence, and provenance.

Implementation-specific private state must not become an implicit downstream
contract.

### 2. External formats enter through adapters

Sanger, NGS, and future external representations are input modalities, not
platform-wide domain models.

Adapters validate and translate external data into canonical evidence while
preserving modality-specific evidence that remains scientifically meaningful.

Canonicalization must not flatten scientifically relevant Sanger or NGS evidence
to a lowest-common-denominator sequence string.

### 3. Public Rust APIs are capability-oriented

The public Rust library API exposes stable capabilities and canonical result
types rather than requiring callers to orchestrate private implementation
stages.

Raw sequencing to canonical variant evidence is treated as one high-level
Variant Analysis capability even though its internal implementation remains
split into independently owned and testable scientific stages.

Future capabilities such as haplogroup and nomenclature consume canonical
artifacts instead of depending directly on Sanger, NGS, a particular aligner, or
a particular caller.

### 4. Replaceability is introduced at real variation seams

Alignment methods, variant-calling algorithms, input adapters, and scientific
data providers may have replaceable implementations when independent variation
is required.

A module is not automatically a plugin. Traits, separate crates, runtime
selection, and provider abstractions are introduced only when they remove a real
coupling, dependency, testing, ownership, or configuration problem.

> **Superseded in part (2026-10-09):** [ADR-0069](0069-plugin-first-modality-core-post-calling.md)
> declares the modality, core-caller, and post-calling families and their data
> contracts as variation seams. Stages inside a plugin remain ordinary modules.

### 5. Algorithm and data-provider variation remain independent

When an analysis algorithm and its knowledge source can vary separately, DNA
models them as separate contracts.

For example, a haplogroup classifier and a haplogroup database are distinct
boundaries; the same principle applies to nomenclature engines and nomenclature
data sources.

### 6. Composition chooses concrete implementations

Concrete implementation selection belongs at an application/composition
boundary.

Canonical domain models and contracts do not depend on concrete Sanger, NGS,
alignment, calling, haplogroup-database, or nomenclature-database
implementations.

The intended dependency direction is:

```text
implementation / provider
          |
          v
canonical contract / domain model
```

### 7. Provenance records implementation identity

Scientifically meaningful results retain enough provenance to identify the
source artifact, reference, configuration, algorithm/component, database or
knowledge source where relevant, and parent artifact identities required for
audit and reproduction.

Downstream components may inspect provenance for audit purposes, but provenance
does not replace an explicit semantic input contract.

### 8. Software compatibility is not scientific equivalence

Two implementations that satisfy the same software contract are
interface-compatible. That fact alone does not establish that they are
scientifically equivalent or interchangeable for a validated claim.

Alternative algorithms and databases require validation appropriate to their
scientific behavior even when swapping them requires no downstream source
changes.

### 9. Single-crate layering remains current until a split is justified

This ADR does not supersede ADR-0002.

The current source may introduce these boundaries inside one crate. A future
workspace ADR may supersede ADR-0002 only after identifying concrete crate
consumers and a cycle-free dependency graph with a justified dependency,
ownership, reuse, build, or release boundary.

> **Superseded (2026-10-09):** [ADR-0070](0070-workspace-split-by-plugin-family.md)
> splits the workspace by plugin family.

## Alternatives considered

### Keep one fixed implementation pipeline

This minimizes abstraction today but couples future input, algorithm, and
database changes to adjacent stages.

### Introduce one universal plugin trait

Rejected because input adaptation, alignment, variant calling, haplogroup
classification, nomenclature, and data provision have different semantics,
inputs, outputs, and validation requirements.

### Turn every current module into a plugin or crate immediately

Rejected as premature abstraction and inconsistent with ADR-0002.

### Use external file formats as internal contracts

Rejected because AB1, BAM/CRAM, VCF, JSON, and similar formats are transport or
storage representations rather than the canonical in-process scientific model.

## Consequences

### Positive

- Input modalities can evolve without teaching downstream analysis about their
  file formats.
- Alternative algorithms can be introduced without changing unrelated modules.
- Haplogroup, nomenclature, and future analysis capabilities can consume the same
  canonical upstream artifacts.
- Rust types can enforce semantic boundaries at compile time.
- Public API evolution is separated from private implementation topology.
- Scientific provenance remains auditable across replaceable implementations.
- A future crate split can follow already-established semantic boundaries.

### Costs

- Canonical contract design becomes a first-class engineering task.
- Public types require deliberate compatibility discipline.
- Adapters may need to preserve both shared semantics and modality-specific
  evidence.
- Alternative implementations require contract tests plus separate scientific
  validation.
- Poorly chosen abstractions can still create accidental coupling, so traits and
  crates must remain evidence-driven.

## Implementation rule

Adopt the architecture incrementally.

For each new boundary:

1. characterize existing behavior before structural refactoring;
2. add or update focused tests first for behavioral changes;
3. define the smallest canonical contract needed by the current use case;
4. keep concrete implementation selection outside the canonical domain model;
5. update architecture/design/reference/source-local documentation in the same
   change when current behavior or public contracts change;
6. validate scientific behavior separately from interface conformance.

Behavioral changes follow Red -> Green -> Refactor TDD. Boundary-preserving
refactors require sufficient characterization tests before restructuring.
