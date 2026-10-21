---
implements: [SRS-API, SRS-CALL, SRS-NOM, SRS-SAMPLE, SRS-VAR, SRS-VN]
decided-by: [ADR-0070, ADR-0076, ADR-0077, ADR-0078, ADR-0079, ADR-0087]
covers:
  - crates/dna-kernel/src/variant.rs
  - crates/dna-normalize/src/
  - src/lib.rs
  - src/model/review_result.rs
  - src/pipeline/review.rs
  - src/pipeline/sample_call.rs
  - src/report/review.rs
  - src/recorded.rs
  - src/sample_analysis/
  - src/variant_analysis/
---

# Rust Public API

This document defines the current stable Rust library boundary exposed by the
`dna` crate.

Sample-level analysis (`SampleAnalysis`) is in [[reference/rust-api-sample]].

## Crate surface

The public modules are `cli`, `error`, `profile`, `recorded`,
`sample_analysis`, `variant`, `variant_analysis`, `variant_normalization`, and
`variant_nomenclature`, plus
the CLI dispatcher `dna::run(cli::Cli) -> dna::error::Result<()>`. All other
modules, including configuration, internal models, scientific stages, and
reporting, are private implementation detail.

The workspace plugin crates (`dna-kernel`; the core's `dna-align`,
`dna-call`, `dna-sample`, `dna-consensus`; `dna-sanger`, `dna-normalize`,
`dna-empop`;
[[ADR-0070]],
[[ADR-0077]],
[[ADR-0079]])
are unpublished implementation crates, not part of this boundary. The `dna`
facade re-exports exactly these public items from them:
- `error` re-exports the failure vocabularies of `dna_kernel::error` (its
  `Error` as `KernelError`) and of `dna_sanger::error`, beside the facade's
  own [`Error`](#errors);
- `variant` is `dna_kernel::variant`;
- `profile` exposes `Profile` and `ProfileIdentity`;
- `recorded` makes a trace's `dna.read_evidence/v1` document in memory and
  analyses a sample from such documents, returning `EvidenceOutcome`
  ([[reference/rust-api-sample|sample analysis]]). `placed_read_document`
  makes a read's `dna.placed_read/v1` document, and `analyze_placed_outcome`
  analyses a sample from `ReadDocuments` (evidence, with a placed read when
  given) without searching the given placements again ([[ADR-0087]],
  proposed);
- `sample_analysis` exposes `SampleEvidence` (`dna-sample`) and `Consensus` (`dna-consensus`) and
  `SangerAttachment` and `SangerRejection` (`dna-sanger`), the records of a
  [sample analysis](#sample-analysis);
- `variant_normalization` exposes `normalize`, `NormalizationPolicy`, and
  `VariantNormalizationResult`;
- `variant_nomenclature` exposes `apply`, `from_normalization`,
  `NomenclatureInput`, and `VariantNomenclatureResult`.

Types follow one evolution rule:

- results that only the library produces (`VariantAnalysisResult`,
  `SampleAnalysis`, `RejectedRead`, `ConsensusCall`, `NotationCall`,
  `VariantNormalizationResult`, `VariantNomenclatureResult`) and the enums
  `VariantKind`, `RejectionCause`, `NormalizationPolicy`, `cli::Command`,
  `error::Error`, `error::KernelError`, and `error::SangerError` are
  `#[non_exhaustive]`, so fields and variants can be added without a breaking
  change; callers read their fields and match enums with a wildcard arm;
- the plugin-crate records re-exported by `sample_analysis` are the crates'
  own evidence records with exhaustive public fields; they evolve with the
  crates, and their field types are named through the crates;
- input types that callers construct (`CalledVariantSet`, `Variant`,
  `ReferenceIdentity`, `ReferenceSegment`, `NomenclatureInput`) and the
  identity record `ProfileIdentity` keep exhaustive public fields;
- `Profile` is opaque: it is built only by `Profile::load`, which validates it.

## Target profiles

```rust
impl Profile {
    pub fn load(path: &Path) -> Result<Profile, dna::error::KernelError>;
    pub fn identity(&self) -> &ProfileIdentity;
}

pub struct ProfileIdentity {
    pub id: String,
    pub sha256: String,
}
```

`dna::profile::Profile::load` reads and validates one [[reference/profiles|target profile]];
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
) -> Result<VariantNormalizationResult, dna::error::KernelError>
```

through:

```rust
dna::variant_normalization::normalize(...)
```

The initial policy is:

```rust
#[non_exhaustive]
pub enum NormalizationPolicy {
    RightAligned { max_trade: usize },
}

impl NormalizationPolicy {
    pub fn for_notation(notation: &Notation) -> Self;
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
(a profile's `indel_placement = "right"`). A deletion that cannot shift
further may also move 3' across substitutions it can replace (the deletion
trade, SRS-VN-008: `16293DEL 16298C` becomes `16293C 16298DEL`), which
rewrites those substitutions, by at most `max_trade` bases
(`for_notation` takes the profile's `normalization.max_trade`). `normalize`
bounds the trade by no window or read; the
DNA workflows bound it by the read's callable bases and the profile's
nomenclature windows. Every accepted movement must reconstruct exactly the
same complete alternate sequence as the source calls.

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
) -> Result<VariantNomenclatureResult, dna::error::KernelError>
```

through:

```rust
dna::variant_nomenclature::apply(...)
```

The reference is loaded with the profile topology and must be the profile's
reference when the profile pins one (`KernelError::Profile(ProfileError::ReferenceMismatch)`).

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
| HVS-II 303-315 poly-C around T310 | C-run lengths at 309/315 (`309.1C`, `315.1C`), the rest left-anchored (`312T 315.1C 315.2C 315.3C`), `311T 315.1C`, `310C 315DEL`, terminal `315.1C` |
| HVS-III 512-524 AC repeat | whole-copy loss or gain at the 3' end plus substitutions (`523DEL 524DEL`, `513A 523DEL 524DEL`, `524.1A 524.2C`) |
| HVS-I 16180-16193 A-tract and poly-C around T16189 | `16183C 16184A 16189C`, `16183DEL`, `16188.1C`, `16182C 16183C 16189C`, `16189C 16190T 16193.1C`, `16189C 16193DEL` |

Variants outside the windows, window haplotypes no rule represents, and a
window an edit crosses keep their normalized form. The DNA workflows apply the
windows per read with `dna_normalize::variant_nomenclature::apply_read`, which
the facade does not re-export: it also keeps the normalized form of a window
whose representation would place a call outside the read's callable
reference intervals, or would rest on a decisive base the read does not read
(SRS-NOM-010).
`dna_normalize::length_convention::ShortestDominant` is the optional
"shortest dominant" length convention (ISFG 2014 3.3.2): from a length
distribution it reads the dominant length, or the next-shorter one when its
share is at least the configured ratio of the dominant's (`[normalization]`,
at most `MAX_RATIO`, 0.85). The review applies it as a labelled suggestion
([[reference/review-lengths]]); the facade does not re-export it.
Decimal strings such as `309.1C` remain an outer notation concern. Sanger
artifact interpretation, sample reconciliation, and NGS-specific behavior are
not implemented by this API.

## Errors

Every fallible facade API returns `dna::error::Result<T>`, an alias for
`Result<T, dna::error::Error>`. `Error` says which crate owns the failure:

| `Error` variant | Wrapped failure |
| --- | --- |
| `Kernel` | `KernelError`, a failure of the shared kernel or the core |
| `Sanger` | `SangerError`, a failure of the Sanger modality |
| `OperationAndLog { operation, logging }` | an operation failure whose error log could not be persisted |

`Kernel` and `Sanger` render transparently. The capabilities re-exported from
the plugin crates (`Profile::load`, `variant_normalization::normalize`,
`variant_nomenclature::apply`) return `Result<T, KernelError>`, which converts
into `Error` with `?`. Each kernel and Sanger variant names a stage and wraps
that stage's own `#[non_exhaustive]` failure enum, re-exported from
`dna::error`:

| Variant | Wrapped failure | Display prefix |
| --- | --- | --- |
| `KernelError::Config` | `ConfigError` | `invalid configuration value:` |
| `KernelError::Profile` | `ProfileError` | `invalid target profile:` |
| `KernelError::Fasta` | `FastaError` | `invalid reference FASTA:` |
| `KernelError::Sequence` | `SequenceError` | `invalid sequence FASTA:` |
| `KernelError::Variants` | `VariantsError` | `invalid variants document:` |
| `KernelError::ConsensusDocument` | `ConsensusDocumentError` | `invalid consensus document:` |
| `KernelError::RunContext` | `RunContextError` | `invalid run context:` |
| `KernelError::Evidence` | `EvidenceError` | `read evidence failed:` |
| `KernelError::Alignment` | `AlignmentError` | `alignment failed:` |
| `KernelError::Variant` | `VariantError` | `variant calling failed:` |
| `KernelError::VariantNormalization` | `NormalizationError` | `variant normalization failed:` |
| `KernelError::VariantNomenclature` | `NomenclatureError` | `variant nomenclature failed:` |
| `KernelError::Sample` | `SampleError` | `sample evidence failed:` |
| `KernelError::Report` | `ReportError` | `failed to assemble analysis report:` |
| `SangerError::Abif` | `AbifError` | `invalid ABIF input:` |
| `SangerError::Basecalling` | `BasecallingError` | `base re-calling failed:` |
| `SangerError::Signal` | `SignalError` | `signal processing failed:` |
| `SangerError::Callability` | `CallabilityError` | `read callability failed:` |
| `SangerError::QualityControl` | `QualityControlError` | `quality control failed:` |
| `SangerError::RunSignal` | `RunSignalError` | `invalid run signal:` |
| `SangerError::ReadEvidence` | `ReadEvidenceError` | `invalid read evidence:` |
| `SangerError::Products` | `ProductsError` | `alignment failed:` or `variant calling failed:` |

`ProductsError` reports one read's Sanger products disagreeing while its core
evidence is assembled (call counts, trim, a missing call or channel peaks),
under the stage prefix the core gave it before the modality owned it, so
command-line and log text is the same whichever crate fails. A Sanger
configuration or evidence-validation failure is the kernel's
(`SangerError::Kernel`).

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

`KernelError::Path` reports a rejected path with a static reason. `Read`, `Log`, and
`Output` keep the `std::io::Error` as their `source`. Third-party parser and serializer
failures (`ConfigParse`, `ProfileParse`, `VariantsParse`, `ConsensusParse`, `RunContextParse`, `Serialize`) are erased to `ForeignError`, so
dependency types never appear in the public API. Stage failures render inline,
`"<prefix> <failure>"`, and are reached by matching rather than through
`std::error::Error::source`.

## Instrumentation

Library capabilities emit `tracing` instrumentation under the targets of the
crates that run each stage (`dna`, `dna_call`, `dna_sanger`) and
never install a subscriber, write files, or print. `analyze_sanger` opens one
`info` span per scientific stage, in order and not nested in one another
(`basecalling`, `signal_processing`, `callability`, `quality_control`,
`read_evidence`, `alignment`, `variant_calling`), and emits one structured completion event per stage plus
`warn` events for removed variant candidates and warning summaries.
`SampleAnalysis::call` emits the stage events and the `nomenclature` span of
the `call` command, without its `input_loading` and `reporting` spans. Span names
are stable; event fields are operational detail and may grow.

`sample_analysis::analyze_sanger` runs each read's stages, with the same spans,
under the dispatcher and inside the span that were current on the calling
thread, so the stage spans of every read have the caller's span as parent on
any worker. Each read's spans and events keep their order; with more than one
worker those of different reads may interleave.

## Boundary rules

The Rust API returns stable typed data for the current Variant Analysis and
Sample Analysis capabilities. The versioned JSON documents under
this reference directory remain separate serialization/publication contracts;
`dna.analysis/v10`, `dna.sample_evidence/v10`, and `dna.consensus/v1` are not
Rust API result models.

The Sanger adapter name is source-specific by design. Future input modalities
may provide additional adapters while converging on compatible
`CalledVariantSet` semantics. Optional normalization and target nomenclature
remain separate capabilities rather than hidden stages of the source adapter.
