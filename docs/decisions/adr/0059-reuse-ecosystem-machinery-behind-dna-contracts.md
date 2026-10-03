# ADR-0059: Reuse ecosystem machinery behind DNA-owned contracts

- **Status:** Accepted
- **Date:** 2026-10-03
- **Related decision:** [ADR-0058](0058-canonical-contracts-and-modular-analysis-composition.md)
- **Related proposal:** [PROP-0001](../../proposals/0001-modular-dna-analysis-platform.md)

## Context

DNA is evolving from one production Sanger path into a modular analysis
platform. That growth introduces capabilities that already have mature Rust or
HTS ecosystem implementations: FASTA/FASTQ I/O, SAM/BAM/CRAM and VCF/BCF
handling, pairwise alignment machinery, pileup traversal, indexing, and related
infrastructure.

Reimplementing those mechanisms inside DNA would increase correctness,
maintenance, interoperability, performance, and supply-chain review burden while
distracting from the scientific semantics that DNA actually needs to own.

At the same time, exposing dependency-specific structures as DNA's domain or
public contracts would couple the platform to one library and make replacement
costly.

ADR-0058 already establishes canonical DNA-owned contracts and dependency
direction. This ADR defines how concrete implementations should be sourced.

## Decision

DNA follows a reuse-first implementation policy.

### 1. Own scientific semantics; reuse commodity machinery

DNA owns:

- canonical domain and result contracts;
- biological and scientific policy;
- coordinate/reference semantics;
- provenance requirements;
- validation rules and scientific claim boundaries;
- composition of capabilities.

DNA should reuse maintained ecosystem implementations for commodity algorithms,
file-format handling, indexing, traversal, and infrastructure when they satisfy
those contracts.

### 2. Evaluate reuse before custom implementation

Before adding a new parser, algorithm, data structure, or infrastructure
component, implementation work MUST evaluate maintained ecosystem options.

The preferred order is:

```text
reuse directly
    ↓
adapt behind a DNA-owned boundary
    ↓
extend a maintained implementation
    ↓
custom implementation only for a documented unmet requirement
```

A custom implementation requires a concrete gap such as a DNA-specific
scientific invariant, unsupported scoring semantics, unacceptable dependency or
platform constraints, or validation evidence showing the available
implementation is unsuitable.

### 3. Standard bioinformatics formats are external responsibilities

DNA MUST NOT create bespoke SAM/BAM/CRAM, VCF/BCF, FASTA, or FASTQ parsers when a
maintained standards-aware library can satisfy the required contract.

Current candidates include:

- `rust-bio` for sequence algorithms and FASTA/FASTQ support;
- `noodles` for pure-Rust FASTA/FASTQ and HTS formats including
  SAM/BAM/CRAM and VCF/BCF;
- `rust-htslib` when HTSlib-backed BAM/CRAM, VCF/BCF, pileup, indexing, or
  interoperability is the better production boundary.

These are evaluation candidates, not dependencies that must be added
speculatively.

### 4. Algorithm libraries remain implementation details

Pairwise alignment and other sequence algorithms SHOULD evaluate maintained
implementations such as `rust-bio` before extending DNA's custom algorithm
surface.

If a maintained algorithm cannot express a required scientific behavior, DNA
may retain or implement the missing mechanism, but the reason must be explicit
and validated.

External library types do not become DNA's canonical domain or public API merely
because the implementation uses that library.

### 5. Dependencies are added only for real consumers

No dependency is added solely to reserve a future architecture.

A dependency enters production only when a current capability uses it and the
change reviews maintenance, security, licensing, platform, performance, and
lockfile impact.

### 6. Reuse does not establish scientific equivalence

Replacing one implementation with another can preserve a software contract
without preserving scientific output.

Contract conformance and scientific validation remain separate obligations under
ADR-0058.

## Alternatives considered

### Implement all algorithms and formats internally

Rejected because it duplicates mature ecosystem work and expands DNA's
maintenance surface without adding DNA-specific scientific value.

### Expose ecosystem types directly as DNA contracts

Rejected because dependency-specific structures would become accidental public
or cross-module contracts.

### Add broad bioinformatics dependencies now for future use

Rejected as speculative dependency growth. Reuse remains demand-driven.

## Consequences

### Positive

- engineering effort concentrates on DNA-specific scientific semantics;
- standard-format interoperability benefits from maintained implementations;
- dependencies remain replaceable behind DNA-owned contracts;
- custom algorithms have an explicit justification threshold;
- future NGS work can use established HTS libraries rather than building format
  infrastructure from scratch.

### Costs

- adapters are required between ecosystem representations and DNA contracts;
- dependency selection still requires maintenance and supply-chain review;
- some DNA-specific requirements may justify keeping custom implementations;
- replacement of an implementation requires separate scientific validation.

## Implementation rule

For each new capability or substantial refactor:

1. define the smallest DNA-owned semantic contract;
2. inventory maintained ecosystem implementations;
3. test whether they satisfy required semantics and operational constraints;
4. adapt the selected implementation behind the DNA boundary;
5. document any remaining custom implementation and the unmet requirement that
   justifies it;
6. validate scientific behavior independently from interface compatibility.
