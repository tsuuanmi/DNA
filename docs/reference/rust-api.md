# Rust Public API

This document defines the current stable Rust library boundary exposed by the
`dna` crate.

## Crate surface

The public modules are `cli`, `error`, `profile`, `variant`, `variant_analysis`,
`variant_normalization`, and `variant_nomenclature`, plus the CLI dispatcher
`dna::run(cli::Cli) -> dna::error::Result<()>`. All other modules, including
configuration, internal models, scientific stages, and reporting, are private
implementation detail.

The workspace plugin crates (`dna-kernel`, `dna-core`, `dna-sanger`,
`dna-post`; [ADR-0070](../decisions/adr/0070-workspace-split-by-plugin-family.md))
are unpublished implementation crates, not part of this boundary. The `dna`
facade re-exports exactly these public items from them:
- `error` is `dna_kernel::error`;
- `variant` is `dna_kernel::variant`;
- `profile` exposes `Profile` and `ProfileIdentity`;
- `variant_normalization` exposes `normalize`, `NormalizationPolicy`, and
  `VariantNormalizationResult`;
- `variant_nomenclature` exposes `apply`, `from_normalization`,
  `NomenclatureInput`, and `VariantNomenclatureResult`.

Types follow one evolution rule:

- results that only the library produces (`VariantAnalysisResult`,
  `VariantNormalizationResult`, `VariantNomenclatureResult`) and the enums
  `VariantKind`, `NormalizationPolicy`, `cli::Command`, and `error::Error` are
  `#[non_exhaustive]`, so fields and variants can be added without a breaking
  change; callers read their fields and match enums with a wildcard arm;
- input types that callers construct (`CalledVariantSet`, `Variant`,
  `ReferenceIdentity`, `ReferenceSegment`, `NomenclatureInput`) and the
  identity record `ProfileIdentity` keep exhaustive public fields;
- `Profile` is opaque: it is built only by `Profile::load`, which validates it.

## Target profiles

```rust
impl Profile {
    pub fn load(path: &Path) -> dna::error::Result<Profile>;
    pub fn identity(&self) -> &ProfileIdentity;
}

pub struct ProfileIdentity {
    pub id: String,
    pub sha256: String,
}
```

`dna::profile::Profile::load` reads and validates one [target profile](profiles.md);
`identity()` returns the declared `id` and the SHA-256 of the file bytes. The
profile's contents are not public fields, so the file format can grow without
an API break.

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
one single-record FASTA reference, and one explicit DNA configuration with the
target profile it names, then runs
the same validated read-processing, alignment, and variant-calling scientific
path used by the CLI.

Filename suffix is not part of the scientific input contract; `.ab1` is a common sequencing filename, while ABIF validity is determined from file content and required tags.

It does not derive or validate a CLI result path, create `results/`, open an
operational log, serialize JSON, or publish a file.

### VariantAnalysisResult

```rust
#[non_exhaustive]
pub struct VariantAnalysisResult {
    pub input_sha256: String,
    pub reference: ReferenceIdentity,
    pub configuration_sha256: String,
    pub profile: ProfileIdentity,
    pub reference_segments: Vec<ReferenceSegment>,
    pub variants: Vec<Variant>,
}
```

`input_sha256`, `reference.sha256`, `configuration_sha256`, and `profile`
identify the exact source artifact, reference sequence, validated configuration
content, and target profile used for the result.

The canonical contract types `ReferenceIdentity`, `Variant`, `VariantKind`, and
`CalledVariantSet` live in `dna::variant`, re-exported from `dna_kernel::variant` (ADR-0070), so
normalization, nomenclature, and future modalities share them without depending
on the Sanger capability. `VariantAnalysisResult` and `ReferenceSegment` live in
`dna::variant_analysis`.

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

#[non_exhaustive]
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
#[non_exhaustive]
pub enum NormalizationPolicy {
    RightAligned,
}
```

The result preserves both source and selected representations:

```rust
#[non_exhaustive]
pub struct VariantNormalizationResult {
    pub reference: ReferenceIdentity,
    pub source_variants: Vec<Variant>,
    pub alternate_sequence: String,
    pub normalized_variants: Vec<Variant>,
}
```

`RightAligned` selects sequence-equivalent 3'/right-most insertion and
deletion representations without rotating across the FASTA coordinate seam
(a profile's `indel_placement = "right"`). Every accepted movement must reconstruct exactly the same complete
alternate sequence as the source calls.

Normalization is optional and does not apply nomenclature windows,
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

Nomenclature applies the windows of a target profile:

```rust
pub fn apply(
    reference_path: &Path,
    profile: &Profile,
    input: NomenclatureInput<'_>,
) -> dna::error::Result<VariantNomenclatureResult>
```

through:

```rust
dna::variant_nomenclature::apply(...)
```

The reference is loaded with the profile topology and must be the profile's
reference when the profile pins one (`Error::Profile(ProfileError::ReferenceMismatch)`).

The result preserves every prior representation and adds the selected target
representation:

```rust
#[non_exhaustive]
pub struct VariantNomenclatureResult {
    pub reference: ReferenceIdentity,
    pub source_variants: Vec<Variant>,
    pub alternate_sequence: String,
    pub normalized_variants: Vec<Variant>,
    pub represented_variants: Vec<Variant>,
}
```

Within each profile window, the first rule whose candidate reconstructs the
identical window haplotype is used (SRS-NOM-004 to SRS-NOM-015). The shipped
human-mtDNA profile declares:

| Window | Forms |
| --- | --- |
| HVS-II 303-315 poly-C around T310 | C-run lengths at 309/315 (`309.1C`, `315.1C`), `311T 315.1C`, `310C 315DEL`, EMPOP terminal `315.1C` |
| HVS-III 513-524 AC repeat | one-motif loss as `523DEL 524DEL` plus substitutions, e.g. `513A` |
| HVS-I 16181-16193 poly-C around T16189 | `16183C 16184A 16189C`, `16189C 16193DEL` |

Variants outside the windows and window haplotypes no rule represents keep
their normalized form; an edit crossing a window boundary fails explicitly.
Decimal strings such as `309.1C` remain an outer notation concern. Sanger
artifact interpretation, sample reconciliation, and NGS-specific behavior are
not implemented by this API.

## Errors

Every fallible API returns `dna::error::Result<T>`, an alias for
`Result<T, dna::error::Error>`. `Error` names the failing stage, and each stage
variant wraps that stage's own `#[non_exhaustive]` failure enum, re-exported from
`dna::error`:

| `Error` variant | Wrapped failure | Display prefix |
| --- | --- | --- |
| `Config` | `ConfigError` | `invalid configuration value:` |
| `Profile` | `ProfileError` | `invalid target profile:` |
| `Abif` | `AbifError` | `invalid ABIF input:` |
| `Fasta` | `FastaError` | `invalid reference FASTA:` |
| `Sequence` | `SequenceError` | `invalid sequence FASTA:` |
| `Variants` | `VariantsError` | `invalid variants document:` |
| `Basecalling` | `BasecallingError` | `base re-calling failed:` |
| `Signal` | `SignalError` | `signal processing failed:` |
| `Callability` | `CallabilityError` | `read callability failed:` |
| `QualityControl` | `QualityControlError` | `quality control failed:` |
| `Evidence` | `EvidenceError` | `read evidence failed:` |
| `Alignment` | `AlignmentError` | `alignment failed:` |
| `Variant` | `VariantError` | `variant calling failed:` |
| `VariantNormalization` | `NormalizationError` | `variant normalization failed:` |
| `VariantNomenclature` | `NomenclatureError` | `variant nomenclature failed:` |
| `Sample` | `SampleError` | `sample evidence failed:` |
| `Report` | `ReportError` | `failed to assemble analysis report:` |

Failures that depend on observed values carry them as structured fields, for
example `AlignmentError::LowIdentity { identity, minimum }`,
`CallabilityError::TooFewCallableCalls { callable, minimum }`, or
`FastaError::UnsupportedBase { base }`. Fixed rules are distinct variants, some
of which carry the rule's static text (for example `ConfigError::Constraint` or
`AbifError::PeakLocations`). Arithmetic overflow and internal-consistency
guards use `Overflow(&'static str)` or `Inconsistent(&'static str)`; they
indicate a defect rather than bad input. `NormalizationError` and
`NomenclatureError` wrap the shared `RepresentationError` for allele, edit, and
anchoring failures.

`Path` reports a rejected path with a static reason. `Read`, `Log`, and
`Output` keep the `std::io::Error` as their `source`. Third-party parser and serializer
failures (`ConfigParse`, `ProfileParse`, `VariantsParse`, `Serialize`) are erased to `ForeignError`, so
dependency types never appear in the public API. Stage failures render inline,
`"<prefix> <failure>"`, and are reached by matching rather than through
`std::error::Error::source`.

## Instrumentation

Library capabilities emit `tracing` instrumentation under the targets of the
crates that run each stage (`dna`, `dna_core`, `dna_sanger`) and
never install a subscriber, write files, or print. `analyze_sanger` opens one
`info` span per scientific stage, in order and not nested in one another
(`basecalling`, `signal_processing`, `callability`, `quality_control`,
`read_evidence`, `alignment`, `variant_calling`), and emits one structured completion event per stage plus
`warn` events for removed variant candidates and warning summaries. Span names
are stable; event fields are operational detail and may grow.

## Boundary rules

The Rust API returns stable typed data for the current Variant Analysis capability. The versioned JSON documents under
this reference directory remain separate serialization/publication contracts;
`dna.analysis/v9` is not the Rust API result model.

The Sanger adapter name is source-specific by design. Future input modalities
may provide additional adapters while converging on compatible
`CalledVariantSet` semantics. Optional normalization and target nomenclature
remain separate capabilities rather than hidden stages of the source adapter.
