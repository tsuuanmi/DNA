---
id: PROP-0001
type: proposal
status: implementing
owners: []
created: 2026-10-03
related-requirements: []
related-decisions: [ADR-0058]
implementation: [PR-9, PR-10, PR-11, PR-12, PR-15, PR-16]
---

# Proposal: Modular DNA Analysis Platform

## Problem

DNA currently implements a production-oriented Sanger AB1 analysis path whose
stages are deliberately explicit: trace decoding, base calling, signal
processing, quality control, alignment, variant calling, sample evidence, and
reporting.

The project is expected to grow beyond one fixed Sanger pipeline. Future use
cases include:

- multiple input modalities such as Sanger and NGS;
- multiple alignment methods;
- multiple variant-calling algorithms;
- downstream DNA-analysis capabilities such as haplogroup and nomenclature;
- alternative databases or knowledge sources used by those downstream
  capabilities;
- a stable public Rust API in addition to the CLI.

If implementations exchange implementation-specific types directly, replacing
one component will force changes into adjacent components. The same problem
appears when external file formats, algorithms, databases, and public API
contracts are allowed to become one boundary.

The architecture therefore needs stable semantic contracts between capabilities
while preserving the evidence and provenance needed for scientific auditability.

## Goals

- Make DNA usable as a Rust library through a deliberately versioned public API.
- Treat raw-to-variant analysis as one high-level capability while preserving
  smaller internal scientific stages.
- Allow input adapters, alignment methods, variant callers, and scientific data
  providers to vary independently where there is a real replacement need.
- Let downstream analysis modules consume canonical DNA data rather than depend
  on Sanger, NGS, a particular aligner, or a particular caller.
- Preserve modality-specific evidence instead of reducing all inputs to the
  lowest common denominator.
- Record enough provenance to identify the input, reference, implementation,
  configuration, and database that produced a scientific result.
- Keep component compatibility distinct from scientific equivalence: satisfying
  a software contract does not establish that two algorithms are scientifically
  interchangeable.
- Allow today's single-crate implementation to evolve toward stronger crate
  boundaries only when those boundaries are justified.

## Non-goals

- Converting every current source module into a trait or crate immediately.
- Runtime dynamic library loading.
- Claiming that different aligners, variant callers, or databases produce
  equivalent scientific results.
- Defining final NGS, haplogroup, nomenclature, or population-analysis methods in
  this proposal.
- Making experimental capabilities part of the current production contract.
- Replacing the current validated Sanger path before equivalent requirements,
  tests, and validation evidence exist.

## Evidence and research

The current source already separates scientific responsibilities behind small
module APIs and routes dependencies toward shared model/config/error boundaries.
The existing documentation architecture also distinguishes current architecture
from proposals and requires exact interface shapes to live in reference
documentation once they become current.

This proposal builds on those boundaries rather than replacing them with a
plugin framework pre-emptively.

## Proposed design

### Architectural principle

The central rule is:

> Implementations may vary; the semantic contract between capabilities remains
> stable.

A component depends on the input and output contracts it needs, not on the
concrete implementation that produced the input or consumes the output.

Provenance may record those concrete implementations for auditability without
turning that information into coupling between modules.

### Platform model

The long-term system is an analysis platform centered on canonical DNA data:

```text
External inputs
    |
    +-- Sanger
    +-- NGS
    +-- future formats
    |
    v
Input adapters
    |
    v
Canonical evidence
    |
    v
Variant Analysis
    |
    v
Canonical variant evidence
    |
    +----------------+-------------------+
    |                |                   |
    v                v                   v
Haplogroup       Nomenclature       Future analysis
    |                |                   |
    v                v                   v
Canonical results / artifacts
```

The platform is not defined by one linear pipeline. As capabilities grow, it may
form a graph of analysis modules operating on canonical artifacts.

### Capability, module, plugin, and provider

These terms have different meanings:

- **capability**: a stable operation the platform offers, such as variant
  analysis, haplogroup classification, or nomenclature;
- **module**: the product/domain boundary that owns one capability and its
  scientific rules;
- **plugin implementation**: a replaceable implementation behind an explicit
  capability contract when independent variation is actually required;
- **provider**: a replaceable source of external or reference data, such as a
  haplogroup database.

A module does not need a plugin interface merely because it is modular.
Abstractions are introduced when independent implementations, dependency
isolation, testing seams, or configuration-driven selection make them useful.

### High-level variant analysis

Raw sequencing data to canonical variant evidence is exposed as one high-level
capability:

```text
raw sequencing
      |
      v
+-------------------------+
| Variant Analysis        |
|                         |
| decode                  |
| signal processing       |
| noise handling          |
| base calling            |
| quality control         |
| alignment               |
| variant calling         |
+------------+------------+
             |
             v
canonical variant evidence
```

Internally these stages remain independently owned and testable. Public callers
are not required to orchestrate every stage themselves.

The current alignment and variant-calling logic remain a scientific kernel
inside this capability. Earlier signal/noise/base-calling stages and later
analysis capabilities may evolve without redefining that kernel as the whole
platform core.

### Canonical contracts

Each replaceable boundary has an explicit typed input and output contract.

Conceptually:

```text
InputAdapter:
external input -> canonical evidence

Aligner:
alignment input -> canonical alignment

VariantCaller:
variant-calling input -> canonical variant set

HaplogroupClassifier:
canonical sample/variant evidence -> canonical haplogroup result

NomenclatureEngine:
canonical variant evidence -> canonical nomenclature result
```

Concrete names and exact Rust shapes are intentionally deferred until
implementation. Once public, exact shapes belong under `docs/reference/`.

Canonical contracts define semantic meaning, including:

- coordinate and strand conventions;
- required and optional evidence;
- unresolved/unknown states;
- reference identity;
- result cardinality and ordering where applicable;
- provenance requirements.

Implementation-specific private state must not leak into downstream contracts.

### Preserve modality-specific evidence

Canonicalization does not mean flattening Sanger and NGS into an impoverished
common representation.

For example, Sanger may carry chromatogram peaks and locus evidence while NGS
may carry read depth, per-base quality, mapping quality, and strand support.
Shared contracts expose common semantics while retaining source-specific
evidence where it remains scientifically meaningful.

A downstream component should request only the capability it needs. An aligner
should not need to know that its sequence evidence came from AB1 or BAM/CRAM
unless that distinction is scientifically part of the contract.

### Independent algorithm selection

Alignment and variant calling are independently replaceable dimensions.

Conceptually:

```text
                    Aligner contract
                    /      |       \
               Gotoh      WFA      future

                VariantCaller contract
                  /       |        \
              caller A  caller B   future
```

A pipeline or composition layer chooses implementations. Scientific modules do
not select concrete dependencies internally unless the implementation is itself
part of the module's fixed scientific policy.

### Separate algorithms from data providers

Downstream analysis may also vary in two independent dimensions: algorithm and
knowledge source.

For example:

```text
HaplogroupClassifier
        |
        v
HaplogroupDatabase
   /          \
database A   database B
```

The classifier contract and database/provider contract should remain distinct
when they can vary independently.

The same rule applies to nomenclature engines and nomenclature databases.

### Public API

DNA should expose a stable Rust library API organized around capabilities rather
than internal implementation stages.

A future high-level API may conceptually resemble:

```rust
let analysis = dna::variant_analysis::analyze(input, reference, config)?;
let haplogroup = dna::haplogroup::analyze(&analysis, haplogroup_config)?;
let nomenclature = dna::nomenclature::analyze(&analysis, nomenclature_config)?;
```

This example is illustrative, not an accepted exact API.

The public API should expose:

- high-level analysis capabilities;
- stable canonical domain/result types needed by callers;
- configuration required to control supported behavior;
- typed public errors;
- selected expert-level APIs only when there is a demonstrated use case.

Internal helpers and implementation topology remain private by default.

### Composition layer

Concrete implementation selection belongs at a composition boundary rather than
inside the domain model.

Conceptually:

```text
                   Public API / application
                            |
                            v
                     Composition layer
                            |
          +-----------------+-----------------+
          |                 |                 |
          v                 v                 v
     input adapter        aligner       variant caller
          |                                   |
          +-----------------+-----------------+
                            |
                            v
                    canonical contracts
```

This keeps the canonical model independent of concrete Sanger, NGS, alignment,
calling, or database implementations.

### Dependency direction

The desired dependency direction is:

```text
implementation/provider
          |
          v
canonical contracts / domain model
```

A future workspace may eventually separate crates, but `dna-core`-style
contract/domain code must not depend on concrete implementation crates.

The current single-crate layout may remain while these boundaries mature.
Module boundaries should be designed so a justified future crate extraction does
not require redesigning the semantic contracts.

### Provenance

Every scientifically meaningful result must retain enough identity to reproduce
and audit how it was produced.

Depending on the capability, provenance may include:

- source artifact identity and checksum;
- reference identity and checksum;
- component/algorithm identifier and version;
- configuration identity and checksum;
- database/knowledge-source identifier and version;
- relevant parent artifact identities.

Provenance records implementation identity; downstream modules must not use it
as an implicit substitute for an explicit input contract.

### Compatibility versus scientific equivalence

Two implementations that satisfy the same Rust/API contract are interface
compatible. They are not automatically scientifically equivalent.

New implementations require validation appropriate to the claims they make.
Selection of another aligner, caller, or database may change scientific results
even when no downstream source code changes.

Validation and release evidence therefore remain separate from software
substitutability.

### Possible future crate shape

If independent dependency or release boundaries justify it, the repository may
evolve toward a workspace such as:

```text
crates/
├── dna-core/
├── dna-variant-analysis/
├── dna-input-sanger/
├── dna-input-ngs/
├── dna-haplogroup/
├── dna-nomenclature/
└── dna/                    # public facade
```

This is a direction, not a required migration plan. Crates are introduced only
when they provide a real dependency, ownership, reuse, release, or compilation
boundary.

## Alternatives considered

### Keep one fixed pipeline

This is simpler while DNA only serves the current Sanger path, but it makes
future modality and algorithm changes progressively more coupled.

### Make every current stage a trait/plugin now

Rejected as premature abstraction. Most current stages have one authoritative
implementation. Traits and separate crates should be earned by real independent
variation or boundary needs.

### Give every implementation its own result model

Rejected because it transfers integration complexity downstream. Canonical
contracts should absorb representation differences while preserving scientifically
important source evidence.

### Flatten all inputs to one minimal sequence representation

Rejected because it would discard modality-specific evidence that may be needed
for scientific interpretation, validation, and auditability.

### Use file formats as the internal integration contract

Rejected because public/external formats such as AB1, BAM/CRAM, VCF, or JSON are
transport/storage concerns. Internal scientific contracts should express typed
domain semantics directly.

## Validation plan

Architecture implementation should be introduced incrementally.

For each new replaceable boundary:

1. characterize the existing behavior with focused tests before refactoring;
2. define the canonical contract and invariants;
3. add contract tests that every implementation must satisfy;
4. verify that swapping implementations does not require changes to unrelated
   modules;
5. separately validate scientific behavior and equivalence claims using approved
   datasets where applicable;
6. keep serialized public contracts versioned and validated.

Behavioral source changes follow the repository's Red -> Green -> Refactor TDD
workflow. Pure boundary-preserving refactors require adequate characterization
tests before restructuring.

## Rollout and rollback

The current Sanger implementation remains the authoritative production path
while the architecture is introduced.

Recommended rollout sequence:

1. accept the architectural direction and public-contract principles;
2. define the initial public Rust API without changing scientific behavior;
3. make raw-to-variant analysis a high-level capability boundary;
4. introduce canonical contracts at demonstrated replacement seams;
5. add additional implementations such as NGS input, alternative alignment, or
   alternative calling independently;
6. add downstream analysis modules such as haplogroup and nomenclature;
7. extract crates only when their dependency or lifecycle boundary is clear.

Each implementation PR should be independently revertible and must not require a
big-bang migration.

## Observability and operations

Operational logs should identify the selected implementations and providers when
selection becomes variable. Scientific result provenance remains the canonical
audit trail; operational logs do not replace result provenance.

No new runtime service or distributed operational boundary is introduced by
this proposal.

## Security and privacy

New input adapters must preserve the existing bounded-input and untrusted-data
rules. External databases/providers introduce new dependency and supply-chain
boundaries and require review before becoming production dependencies.

Canonical contracts should avoid unnecessarily copying or exposing sensitive
sample metadata. Scientific identifiers and provenance should contain only the
information required for reproducibility and auditability.

## Compatibility and migration

No compatibility layer is required solely to preserve today's private internal
module topology.

Public Rust APIs introduced under this proposal must be treated as intentional
contracts. Incompatible public API or serialized-result changes require the
repository's normal versioning and migration policy.

The current CLI and JSON contracts remain unchanged until a later implementation
explicitly updates their requirements/reference documentation.

## Decision

Accepted on 2026-10-03.

The durable architectural choices around canonical contracts, dependency
direction, public API ownership, replaceable implementations, provider
boundaries, and provenance are recorded by
[ADR-0058](../decisions/adr/0058-canonical-contracts-and-modular-analysis-composition.md).

Acceptance establishes architectural direction; it does not make unimplemented
capabilities current production behavior.

## Implementation

Implementation is in progress through focused PRs.

- [PR #9](https://github.com/tsuuanmi/DNA/pull/9) decouples the pipeline operation
  boundary from CLI/`clap` argument types without changing scientific behavior.

- [PR #10](https://github.com/tsuuanmi/DNA/pull/10) makes configuration-path
  selection explicit at the outer application boundary so pipeline operations
  no longer depend on process-global environment selection.

- [PR #11](https://github.com/tsuuanmi/DNA/pull/11) separates single-read
  scientific input loading from JSON publication-path validation so reusable
  analysis is not coupled to an existing `results/*.json` target.

- [PR #12](https://github.com/tsuuanmi/DNA/pull/12) introduces the minimal
  internal `StageLog` capability so shared read/observation science no longer
  depends on the concrete file-backed operational logger.

- [PR #15](https://github.com/tsuuanmi/DNA/pull/15) introduces the first public
  capability-oriented Rust API: Sanger input to canonical typed Variant Analysis
  results without CLI logging or JSON publication side effects.

- [PR #16](https://github.com/tsuuanmi/DNA/pull/16) moves shared scientific
  ownership out of the CLI-oriented pipeline: reference-free read processing is
  crate-internal shared science and reference-guided observation belongs to
  Variant Analysis.

Each implementation PR must update current
architecture/design/reference/source-local documentation in the same change when
it changes current truth.
