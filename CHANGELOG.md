# Changelog

All notable user-visible and production-readiness changes are recorded here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
DNA uses semantic versions for tagged releases, while scientific JSON contracts
remain independently versioned and are never silently changed in place.

## [Unreleased]

### Added

- Target profiles (ADR-0063, SRS-PRF-001 to SRS-PRF-007): target knowledge —
  reference identity and topology, reportable regions, indel placement,
  nomenclature windows and their rules, notation style — moves out of code and
  configuration into strict versioned TOML profiles. The shipped profile is
  `config/profiles/human-mtdna-rcrs.toml`; a profile that pins a reference fails
  closed on another sequence. Scientific output for human mtDNA is unchanged.
  **Breaking:** configuration schema 6 gains a root `profile` key and loses
  `[reference]` and `variant_calling.regions`; `dna.analysis/v8` and
  `dna.sample_evidence/v9` record `provenance.profile` (`id`, file `sha256`);
  sample `notation.policy` becomes `notation.style = "per_base_decimal"`.
  **Breaking (Rust API):** new `dna::profile::{Profile, ProfileIdentity}`;
  `variant_nomenclature::mtdna::apply_control_region` is replaced by
  `variant_nomenclature::apply(reference, &profile, input)`;
  `NormalizationPolicy::MtDnaRightAligned` is renamed `RightAligned`;
  `VariantAnalysisResult` gains `profile`; `Error` gains `Profile` and
  `ProfileParse`; `ConfigError::RegionOutOfBounds` moves to `ProfileError`;
  `NomenclatureError` window names are `String`s.
- Read callability (ADR-0062, SRS-VAR-013): variants with a call near either
  read end or right after a long homopolymer are ineligible with reasons
  `read_end` / `post_homopolymer`. **Breaking:** configuration schema 6 adds
  `variant_calling.read_end_margin`, `homopolymer_min_length`, and
  `post_homopolymer_window`.
- `dna.sample_evidence/v9`: when the target profile declares notation, an
  optional `notation` view publishes each read's eligible calls after the
  profile's right alignment and nomenclature windows (HVS-II, HVS-III, HVS-I for
  human mtDNA), rendered per base (`73G`, `249DEL`, `309.1C`) with supporting
  reads (SRS-NOM-010 to SRS-NOM-012). All other v8 fields are unchanged.
- Production-oriented CI/CD and supply-chain verification, including strict cargo-shear dependency/source hygiene and an explicit Ubuntu 24.04 runner baseline.
- Declare the crate proprietary (`license = "LicenseRef-Proprietary"`) so SBOMs
  and audits record its ownership.
- Crate-wide Rust lint policy in `Cargo.toml` (`missing_docs`, `unreachable_pub`,
  Clippy `pedantic`) and a tuned release profile (thin LTO, one codegen unit).
- Rust/Python source-boundary enforcement.
- Dependency policy, dependency review, CodeQL, OpenSSF Scorecard, fuzzing,
  self-contained explicit-target auditable release bundles with post-strip metadata verification, release SBOMs, and artifact attestations.

### Changed

- **Breaking (Rust API):** `variant_nomenclature::mtdna::apply_hv2_polyc` is
  replaced by `apply_control_region`, which adds HVS-II `311T 315.1C`,
  `310C 315DEL` and EMPOP `315.1C`, HVS-III `513A 523DEL 524DEL`, and HVS-I
  `16183C 16184A 16189C` / `16189C 16193DEL` (SRS-NOM-013 to SRS-NOM-015).
  `NomenclatureError` window failures name their window.
- **Breaking (Rust API):** the empty public `dna::config` and `dna::model`
  modules are now private; library-produced results (`VariantAnalysisResult`,
  `VariantNormalizationResult`, `VariantNomenclatureResult`) and the enums
  `VariantKind`, `NormalizationPolicy`, `cli::Command`, and `error::Error` are
  `#[non_exhaustive]`.
- Operational logging uses `tracing`: scientific stages emit structured events
  and per-stage spans, and a DNA layer renders them in the unchanged
  per-operation record format (ADR-0061). Library callers can observe stage
  progress through their own `tracing` subscriber.
- **Breaking (Rust API):** `dna::error::Error` stage variants now wrap typed,
  `#[non_exhaustive]` per-stage failure enums (`AbifError`, `FastaError`,
  `ConfigError`, `AlignmentError`, and so on) instead of `String` messages;
  `DNAProcessing` is renamed `Signal`; `Path::reason` is `&'static str`; and
  `ConfigParse`/`Serialize` no longer expose `toml`/`serde_json` types. CLI and
  log error text is unchanged, except two internal-invariant messages that
  valid input cannot reach: `VariantError::NoSupportingCalls` uses the
  uppercase variant label, and a non-ASCII rendered allele reports
  `RepresentationError::InvalidAllele`.
- CI documentation now matches the workflows: CodeQL, fuzzing, MSRV, release
  build, and RustSec audit run on `main` and on schedules, not on pull requests.
- Corrected the minimum supported Rust version to match language features used
  by the codebase.

### Fixed

- An indel flanking call whose primary base is unresolved (`N`) is omitted from
  public calls (SRS-VAR-006). Previously a tied flank aborted `analyze` or a
  whole `sample` operation, and a mixed-signal flank was published with the
  schema-invalid base `N`.
- PLOC loci one sample apart, which SRS-IN-003 accepts, no longer fail
  basecalling with an empty locus window.

## [0.1.0] - Unreleased

Initial pre-production development baseline.
