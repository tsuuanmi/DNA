# Rust Public API

This document defines the current stable Rust library boundary exposed by the
`dna` crate.

## Crate surface

The public modules are `cli`, `error`, `profile`, `sample_analysis`, `variant`,
`variant_analysis`, `variant_normalization`, and `variant_nomenclature`, plus
the CLI dispatcher `dna::run(cli::Cli) -> dna::error::Result<()>`. All other
modules, including configuration, internal models, scientific stages, and
reporting, are private implementation detail.

The workspace plugin crates (`dna-kernel`, `dna-core`, `dna-sanger`,
`dna-post`; [ADR-0070](../decisions/adr/0070-workspace-split-by-plugin-family.md))
are unpublished implementation crates, not part of this boundary. The `dna`
facade re-exports exactly these public items from them:
- `error` is `dna_kernel::error`;
- `variant` is `dna_kernel::variant`;
- `profile` exposes `Profile` and `ProfileIdentity`;
- `sample_analysis` exposes `SampleEvidence` and `Consensus` (`dna-core`) and
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
  `VariantKind`, `RejectionCause`, `NormalizationPolicy`, `cli::Command`, and
  `error::Error` are
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

## Sample Analysis

One placement pass over a sample's Sanger traces yields both the sample
evidence and the consensus:

```rust
pub fn analyze_sanger(
    traces: &[PathBuf],
    reference: &Path,
    config: &Path,
    workers: NonZeroUsize,
) -> dna::error::Result<SampleAnalysis>
```

through:

```rust
dna::sample_analysis::analyze_sanger(...)
```

It loads the traces, reference, configuration, and target profile as the
`sample` and `consensus` commands do, processes and places every trace once,
sets aside rejected reads (SRS-SAMPLE-027, SRS-SAMPLE-029, SRS-SAMPLE-030), aggregates the
admitted reads, and builds their consensus. The two commands run this same
pass, so the result describes exactly what they publish. A caller that needs
both no longer places every trace twice.

Reads are processed on up to `workers` scoped threads; with one worker they run
one after another on the calling thread. A caller that already runs many
samples in parallel passes `NonZeroUsize::MIN`; the CLI passes the available
parallelism. The result does not depend on the number of workers or on
scheduling. When a read fails other than by rejection, the failure of the
first such read in trace order is reported, as a sequential pass would: its
error is returned, or its panic (a defect) resumed. Workers stop taking reads
after the lowest failed one. Like `analyze_sanger` for one read, it opens no
log, derives no output path, serializes nothing, and publishes nothing.

```rust
#[non_exhaustive]
pub struct SampleAnalysis {
    pub reference: ReferenceIdentity,
    pub profile: ProfileIdentity,
    pub evidence: SampleEvidence,
    pub consensus: Consensus,
    pub reads: BTreeMap<String, SangerAttachment>,
    pub rejected: BTreeMap<String, RejectedRead>,
}

impl SampleAnalysis {
    pub fn regions(&self) -> &[[usize; 2]];
    pub fn loci(
        &self,
        first_1based: usize,
        last_1based: usize,
    ) -> dna::error::Result<Vec<SampleLocusEvidence>>;
}

#[non_exhaustive]
pub struct RejectedRead {
    pub input_name: String,
    pub input_sha256: String,
    pub cause: RejectionCause,
}

#[non_exhaustive]
pub enum RejectionCause {
    TooShort(ShortRead),
    Callability { sanger: SangerRejection, rejection: ReadRejection },
    Placement { sanger: SangerRejection, rejection: PlacementRejection },
}

impl RejectionCause {
    pub const fn sanger(&self) -> Option<&SangerRejection>;
}
```

- `evidence` is the core's sample evidence (`dna_core`): the read registry in
  content-SHA-256 order, the rejected reads, coverage, overlaps, locus
  differences, and normalized variants with support and opposition. It carries
  the reference and configuration SHA-256.
- `consensus` is the core's consensus (`dna_core`): segments, sites, and
  summary; each site names its reads in registry order.
- `reads` holds each admitted read's Sanger evidence (`dna_sanger`): calls,
  signal analysis, callability, and quality control.
- `regions()` returns the reportable regions of the target profile the
  analysis ran with (inclusive, 1-based), so a caller reads them from the
  profile rather than repeating them.
- `loci(first_1based, last_1based)` returns, for each position of the
  inclusive interval that an admitted read covers, a `SampleLocusEvidence`
  (`dna_core`) built exactly as `evidence.locus_differences` builds one, but
  whether or not a read differs there: every covering read's state,
  reference-strand base, source call index, and reference-oriented A/C/G/T
  profile, with read indices into `evidence.reads` (SRS-API-013). It is
  computed on demand from the placed reads the analysis keeps and changes
  nothing; uncovered positions have no record. The sample evidence keeps only
  differential loci; a caller that needs the minor channel at a position where
  every read calls the reference reads it here.
- `rejected` holds each rejected read's name, identity, and cause: fewer
  calls than one signal window, too few callable calls, or a placement the
  core rejected. The last two carry the read's integrity and callability
  (`SangerRejection`); a read too short to process has none.

`reads` and `rejected` are keyed by read content SHA-256, the identity the
core records use. Every coordinate field names its base (`_0based`,
`_1based`, `_0based_exclusive`). The published projections of these records
are `dna.sample_evidence/v10` and `dna.consensus/v1`. Notation is not part of
this result; the consensus call below derives it.

### Sample outcome

`analyze_sanger` fails with `SampleError::NoAdmittedReads` when every read is
rejected, so the caller learns how many but not which. A caller that must name
those traces, such as the evaluation harness asking for a rerun, uses the same
pass through:

```rust
pub fn analyze_sanger_outcome(
    traces: &[PathBuf],
    reference: &Path,
    config: &Path,
    workers: NonZeroUsize,
) -> dna::error::Result<SampleOutcome>

pub enum SampleOutcome {
    Analyzed(Box<SampleAnalysis>),
    Rejected(BTreeMap<String, RejectedRead>),
}
```

It returns `Analyzed` with exactly the analysis `analyze_sanger` returns, or
`Rejected` with every read, its cause, and its Sanger evidence, keyed by read
content SHA-256, when none is admitted. Every other failure is the error
`analyze_sanger` returns; `analyze_sanger` is this call with `Rejected` turned
into that error. The two cases are complete, so the enum is not marked
non-exhaustive.

### Consensus call

`SampleAnalysis::call` completes the recommended per-sample workflow,
`consensus` then `call` ([workflow](consensus.md#workflow)), in process:

```rust
impl SampleAnalysis {
    pub fn call(&self, sample_id: &str) -> dna::error::Result<ConsensusCall>;
}
```

It validates the sample identifier, renders the consensus FASTA `dna
consensus` would publish for that sample, and calls its records as `dna call`
reads them from that file, with the configuration, profile, and reference of
the analysis. It fails where `call` would fail, for example when the consensus
has no segment. Like `analyze_sanger`, it opens no log, derives no output path,
serializes nothing, and publishes nothing. A caller that runs many samples in
parallel calls `analyze_sanger` with `NonZeroUsize::MIN` workers, then `call`.

```rust
#[non_exhaustive]
pub struct ConsensusCall {
    pub callable_segments: Vec<ReferenceSegment>,
    pub notation: Option<Vec<NotationCall>>,
}

#[non_exhaustive]
pub struct NotationCall {
    pub position: usize,
    pub ordinal: usize,
    pub reference: Option<char>,
    pub alternate: Option<char>,
    pub text: String,
}
```

- `callable_segments` holds the reference intervals the consensus sequences
  cover with callable bases (zero-based, half-open), in consensus read order:
  the `callable_reference_segments` of the `dna.variants/v1` reads.
  `ReferenceSegment` is the type of `variant_analysis`, re-exported from
  `sample_analysis`.
- `notation` holds the sample's calls in the profile notation, distinct and in
  reference order: the `notation.calls[].call` values of the `dna.variants/v1`
  document `dna call` would publish. It is `None` when the profile declares no
  notation. Each call gives the 1-based `position` of the changed base (or of
  the base an insertion follows), the insertion `ordinal` (`0` for a
  substitution or deletion), the `reference` base (`None` for an inserted
  base), the called `alternate` base (`None` for a deleted base), and the
  rendered `text` (`16189C`, `309.1C`, `249DEL`).

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
`dna.analysis/v9`, `dna.sample_evidence/v10`, and `dna.consensus/v1` are not
Rust API result models.

The Sanger adapter name is source-specific by design. Future input modalities
may provide additional adapters while converging on compatible
`CalledVariantSet` semantics. Optional normalization and target nomenclature
remain separate capabilities rather than hidden stages of the source adapter.
