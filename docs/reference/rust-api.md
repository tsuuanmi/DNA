# Rust Public API

This document defines the current stable Rust library boundary exposed by the
`dna` crate.

## Variant Analysis

The first capability-oriented API is the Sanger adapter:

```rust
pub fn analyze_sanger(
    trace: &Path,
    reference: &Path,
    config: &Path,
) -> dna::error::Result<VariantAnalysisResult>
```

It is exposed as:

```rust
dna::variant_analysis::analyze_sanger(...)
```

The function validates and loads one AB1 trace, one single-record FASTA
reference, and one explicit DNA configuration, then runs the same canonical
read-processing, alignment, and variant-calling scientific path used by the CLI.

It does not derive or validate a CLI result path, create `results/`, open an
operational log, serialize JSON, or publish a file.

### VariantAnalysisResult

```rust
pub struct VariantAnalysisResult {
    pub input_sha256: String,
    pub reference: ReferenceIdentity,
    pub configuration_sha256: String,
    pub reference_segments: Vec<ReferenceSegment>,
    pub variants: Vec<Variant>,
}
```

`input_sha256`, `reference.sha256`, and `configuration_sha256` identify the
exact source artifact, reference sequence, and validated configuration content
used for the result.

### ReferenceIdentity

```rust
pub struct ReferenceIdentity {
    pub name: String,
    pub sha256: String,
}
```

### ReferenceSegment

```rust
pub struct ReferenceSegment {
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
}
```

Reference segments use zero-based half-open coordinates.

### Variant

```rust
pub struct Variant {
    pub contig: String,
    pub position_1based: usize,
    pub reference: String,
    pub alternate: String,
    pub kind: VariantKind,
}

pub enum VariantKind {
    Snv,
    Ins,
    Del,
}
```

Variant positions are explicitly one-based. Variants are normalized,
configuration-eligible primary-sequence differences produced by the current
scientific path.

## Boundary rules

The Rust API returns canonical typed data. The versioned JSON documents under
this reference directory remain separate serialization/publication contracts;
`dna.analysis/v7` is not the Rust API result model.

The Sanger adapter name is source-specific by design. Future input modalities
may provide additional adapters while preserving the canonical output semantics.
