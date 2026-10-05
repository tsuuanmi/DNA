---
id: PROP-0001
type: proposal
status: implementing
owners: []
created: 2026-10-03
related-requirements: []
related-decisions: [ADR-0058, ADR-0059, ADR-0060]
implementation: [PR-9, PR-10, PR-11, PR-12, PR-15, PR-17, PR-18, PR-19, PR-20, PR-21, PR-22, PR-23, PR-25, PR-26, PR-32]
---

# Proposal: Modular DNA Analysis Platform

## Problem

DNA currently implements a production-oriented Sanger ABIF analysis path whose
stages are deliberately explicit: trace decoding, base calling, signal
processing, quality control, alignment, variant calling, sample evidence, and
reporting.

The project is expected to grow beyond one fixed Sanger pipeline and beyond one
biological target. Future use cases include:

- biological targets such as mitochondrial DNA, nuclear/genomic DNA, and
  targeted loci or panels;
- SNP/genotyping analysis as a variant-focused use case rather than a sequencing
  modality;
- multiple input modalities such as Sanger and NGS;
- multiple external representations such as ABIF, FASTA/FASTQ, BAM/CRAM, and
  VCF/BCF;
- multiple alignment or mapping methods;
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
- Keep the platform independent of one biological target: mtDNA is one current
  application context, not the definition of the DNA platform.
- Model biological target, sequencing modality, external format, and analysis
  capability as separate architectural dimensions.
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
- Reuse maintained ecosystem parsers, algorithms, and HTS infrastructure when
  they satisfy DNA-owned semantic contracts; custom implementations require a
  concrete unmet requirement.

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

DNA separates four dimensions that must not be collapsed into one plugin axis:

```text
biological target   mtDNA | nuclear/genomic DNA | targeted loci/panels
sequencing modality Sanger | NGS | future modalities
external format     ABIF | FASTA/FASTQ | BAM/CRAM | VCF/BCF | ...
capability          basecall | align/map | variant | SNP/genotype | haplogroup | ...
```

SNP is a variant class or analysis use case, not another sequencing modality.

The long-term analysis platform preserves modality-specific evidence until a
semantic boundary genuinely converges:

```text
Sanger formats --> Sanger evidence --> Sanger analysis --> Sanger caller --+
                                                                           |
NGS formats ----> NGS evidence ----> NGS analysis ----> NGS caller --------+--> called variants
                                                                           |         |
future source --> modality evidence --> source-specific caller/adapter -----+         v
                                                                           canonicalization
                                                                                  |
                                                                                  v
                                                                         target nomenclature
                                                                                  |
                                                                                  v
                                                                          canonical variants
                                                                           /      |       \
                                                                          v       v        v
                                                                     sample   haplogroup  SNP/
                                                                   reconcile             targeted
```

Called biological variants are the first likely convergence point because
upstream Sanger chromatogram evidence and future NGS read/alignment/depth
evidence remain scientifically different. Normalization/canonicalization and
target nomenclature are optional capabilities selected by the workflow rather
than mandatory stages of every variant path.

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

### Reuse-first implementation

DNA owns biological semantics, canonical contracts, provenance, and scientific
policy. Commodity bioinformatics machinery should come from maintained
ecosystem implementations when they meet those contracts.

The preferred implementation order is:

```text
reuse -> adapt -> extend -> custom implementation
```

In particular:

- FASTA/FASTQ support should first evaluate `rust-bio` or `noodles`;
- SAM/BAM/CRAM and VCF/BCF support must use maintained HTS implementations such
  as `noodles` or `rust-htslib` rather than bespoke DNA parsers when those
  libraries satisfy the required contract;
- pairwise alignment and related sequence algorithms should evaluate
  `rust-bio` before expanding DNA's custom algorithm surface.

External library structures remain implementation details behind DNA-owned
contracts. Dependencies are added only when a current production capability uses
them; no package is added speculatively.

This policy is recorded by
[ADR-0059](../decisions/adr/0059-reuse-ecosystem-machinery-behind-dna-contracts.md).

### High-level variant analysis

Raw sequencing data to evidence-backed called variants is exposed as one
high-level capability:

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
       called variants
        /      |       \
       v       v        v
 direct use  optional  optional
             normalize nomenclature
                |
                +----> optional nomenclature
```

Internally these stages remain independently owned and testable. Public callers
are not required to orchestrate every upstream read-processing stage themselves.

The current production Variant Analysis capability ends at the called-variant
result. Post-calling normalization/canonicalization and target nomenclature are
explicit optional downstream capabilities defined by ADR-0060 and must be
implemented separately rather than hidden inside the caller.

### Canonical contracts

Each replaceable boundary has an explicit typed input and output contract.

Conceptually:

```text
InputAdapter:
external input -> canonical evidence

Aligner:
alignment input -> canonical alignment

VariantCaller:
modality-specific interpreted evidence -> called variant set

VariantNormalizer:
called variants + reference + explicit policy -> normalization result
(source variants + reconstructed haplotype + normalized variants)

NomenclatureEngine:
canonicalization result + target policy -> target canonical variants

HaplogroupClassifier:
canonical sample/variant evidence -> canonical haplogroup result
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
may carry reads, per-base quality, depth, mapping quality, CIGAR state, and
strand support. These evidence models should remain different until a downstream
operation needs semantics they genuinely share.

The rule is:

> Normalize external formats early. Normalize biological evidence only where the
> biology actually converges.

A downstream component should request only the smallest semantic contract it
needs. Pairwise Sanger alignment, NGS read mapping, and canonical variant
analysis may therefore use different input contracts rather than pretending to
be one universal sequence-evidence interface.

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
let called = dna::variant_analysis::analyze(input, reference, config)?;

// Optional, selected by the workflow/target.
let normalized = dna::variant_normalization::normalize(
    reference_path,
    &called,
    dna::variant_normalization::NormalizationPolicy::MtDnaRightAligned,
)?;
let represented = dna::nomenclature::apply(&normalized, nomenclature_config)?;

// Other consumers may use called or normalized variants directly.
let haplogroup = dna::haplogroup::analyze(&normalized, haplogroup_config)?;
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

Rejected because public/external formats such as ABIF, BAM/CRAM, VCF, or JSON are
transport/storage concerns. Internal scientific contracts should express typed
domain semantics directly.

### Reimplement standard formats and generic algorithms by default

Rejected. DNA should not spend its custom implementation surface on commodity
FASTA/FASTQ, BAM/CRAM, VCF/BCF, indexing, pileup, or generic sequence machinery
when maintained ecosystem implementations meet the required contract. DNA owns
the semantic adapter and scientific policy instead.

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
3. keep raw-to-called-variant analysis as a high-level capability boundary;
4. introduce optional post-calling haplotype normalization/canonicalization as an explicit capability;
5. introduce optional target-specific nomenclature as a separate representation capability;
6. add additional implementations such as NGS input, alternative alignment, or
   alternative calling independently;
7. add downstream analysis modules such as haplogroup and targeted/SNP analysis;
8. extract crates only when their dependency or lifecycle boundary is clear.

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
The reuse-first implementation policy is recorded by
[ADR-0059](../decisions/adr/0059-reuse-ecosystem-machinery-behind-dna-contracts.md).
The variant lifecycle boundary between calling, canonicalization, and
nomenclature is recorded by
[ADR-0060](../decisions/adr/0060-separate-variant-canonicalization-nomenclature.md).

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
  capability-oriented Rust API: Sanger input to typed Variant Analysis
  results without CLI logging or JSON publication side effects.

- [PR #17](https://github.com/tsuuanmi/DNA/pull/17) moves shared scientific
  ownership out of the CLI-oriented pipeline: reference-free read processing is
  crate-internal shared science and reference-guided observation belongs to
  Variant Analysis.

- [PR #18](https://github.com/tsuuanmi/DNA/pull/18) extracts source-specific
  Sanger input loading from the operation pipeline so Variant Analysis and CLI
  workflows consume the same adapter without coupling scientific inputs to
  result-path or publication concerns.

- [PR #19](https://github.com/tsuuanmi/DNA/pull/19) formalizes the format/modality
  boundary: ABIF parsing/decoding is owned by the Sanger input adapter, while the
  decoded `Chromatogram` remains canonical Sanger evidence independent of file
  extension and ABIF container internals.

- [PR #20](https://github.com/tsuuanmi/DNA/pull/20) canonicalizes internal
  Sanger domain vocabulary: ABIF `PLOC.2` is projected to generic locus
  positions at the input boundary, while current versioned JSON contracts retain
  their existing serialized `ploc_*` field names.

- [PR #21](https://github.com/tsuuanmi/DNA/pull/21) broadens the platform scope
  beyond one Sanger/mtDNA context and formalizes reuse-first implementation:
  biological target, sequencing modality, external format, and analysis
  capability are separate dimensions; maintained ecosystem machinery is
  preferred behind DNA-owned contracts.

- [PR #22](https://github.com/tsuuanmi/DNA/pull/22) applies reuse-first to the
  current production path by delegating commodity FASTA record parsing to
  `noodles-fasta` while retaining DNA-owned single-reference semantics and the
  Rust 1.88 MSRV.

- [PR #23](https://github.com/tsuuanmi/DNA/pull/23) records the complementary
  reuse-first outcome for current custom machinery: the evidence-profile Gotoh
  aligner and bounded ABIF decoder remain first-party because evaluated
  ecosystem implementations do not currently satisfy their scientific or
  trust-boundary contracts.

- [PR #25](https://github.com/tsuuanmi/DNA/pull/25) separates evidence-backed
  variant calling from future haplotype-preserving canonicalization and
  target-specific nomenclature, allowing future NGS callers to converge without
  depending on Sanger alignment topology.

- [PR #26](https://github.com/tsuuanmi/DNA/pull/26) implements the first
  cross-modality `CalledVariantSet` boundary plus optional
  haplotype-preserving human-mtDNA 3'/right-most post-calling normalization
  while preserving source variants and reconstructed alternate sequence.

- [PR #32](https://github.com/tsuuanmi/DNA/pull/32) adds the minimal optional
  nomenclature seam over normalization results, retaining reference identity and
  exposing source variants, reconstructed haplotype, and normalized variants
  without implementing target-specific representation rules.

Each implementation PR must update current
architecture/design/reference/source-local documentation in the same change when
it changes current truth.
