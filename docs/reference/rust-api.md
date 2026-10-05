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

The function validates and loads one Sanger sequencing trace encoded as ABIF,
one single-record FASTA reference, and one explicit DNA configuration, then runs
the same validated read-processing, alignment, and variant-calling scientific
path used by the CLI.

Filename suffix is not part of the scientific input contract; `.ab1` is a common sequencing filename, while ABIF validity is determined from file content and required tags.

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

Variant positions are explicitly one-based. Variants are
configuration-eligible, anchored primary-sequence differences produced by the
current scientific path. Their caller representation preserves the
alignment-selected topology.

## CalledVariantSet

`VariantAnalysisResult::called_variants()` projects one analysis result into
the cross-modality called-variant boundary:

```rust
pub struct CalledVariantSet {
    pub reference: ReferenceIdentity,
    pub variants: Vec<Variant>,
}
```

A `CalledVariantSet` means that the variants are evidence-backed calls against
one identified reference. It does **not** mean that the variants have been
right-aligned, left-aligned, VCF-normalized, HGVS-formatted, or passed through
target nomenclature.

## Variant Normalization

Optional post-calling normalization is exposed as:

```rust
pub fn normalize(
    reference_path: &Path,
    called: &CalledVariantSet,
    policy: NormalizationPolicy,
) -> dna::error::Result<VariantNormalizationResult>
```

through:

```rust
dna::variant_normalization::normalize(...)
```

The initial policy is:

```rust
pub enum NormalizationPolicy {
    MtDnaRightAligned,
}
```

The result preserves both source and selected representations:

```rust
pub struct VariantNormalizationResult {
    pub reference: ReferenceIdentity,
    pub source_variants: Vec<Variant>,
    pub alternate_sequence: String,
    pub normalized_variants: Vec<Variant>,
}
```

`MtDnaRightAligned` selects sequence-equivalent 3'/right-most insertion and
deletion representations without rotating across the FASTA/rCRS coordinate
seam. Every accepted movement must reconstruct exactly the same complete
alternate sequence as the source calls.

Normalization is optional and does not apply mtDNA special-region nomenclature,
sample reconciliation, VCF/HGVS formatting, genotype interpretation, or
clinical interpretation.

## Variant Nomenclature

The immutable nomenclature input seam remains:

```rust
pub struct NomenclatureInput<'a> {
    pub reference: &'a ReferenceIdentity,
    pub source_variants: &'a [Variant],
    pub alternate_sequence: &'a str,
    pub normalized_variants: &'a [Variant],
}

pub fn from_normalization(
    normalized: &VariantNormalizationResult,
) -> NomenclatureInput<'_>
```

The first implemented target rule is the human-mtDNA HVS-II 309/315 poly-C
representation:

```rust
pub fn apply_hv2_polyc(
    reference_path: &Path,
    input: NomenclatureInput<'_>,
) -> dna::error::Result<VariantNomenclatureResult>
```

through:

```rust
dna::variant_nomenclature::mtdna::apply_hv2_polyc(...)
```

The result preserves every prior representation and adds the selected target
representation:

```rust
pub struct VariantNomenclatureResult {
    pub reference: ReferenceIdentity,
    pub source_variants: Vec<Variant>,
    pub alternate_sequence: String,
    pub normalized_variants: Vec<Variant>,
    pub represented_variants: Vec<Variant>,
}
```

The HVS-II rule recognizes the validated rCRS 303-315 poly-C window around T310.
Sequence-equivalent anchor movement is represented as C-run length change at the
309 and 315 boundaries while the complete alternate haplotype remains identical.
Decimal strings such as `309.1C` or `315.1C` remain an outer notation concern.

HVS-III 513-524, HVS-I 16189/16193, Sanger artifact interpretation, sample
reconciliation, and NGS-specific behavior are not implemented by this API.

## Boundary rules

The Rust API returns stable typed data for the current Variant Analysis capability. The versioned JSON documents under
this reference directory remain separate serialization/publication contracts;
`dna.analysis/v7` is not the Rust API result model.

The Sanger adapter name is source-specific by design. Future input modalities
may provide additional adapters while converging on compatible
`CalledVariantSet` semantics. Optional normalization and target nomenclature
remain separate capabilities rather than hidden stages of the source adapter.
