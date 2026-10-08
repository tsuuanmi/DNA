# Changelog

All notable user-visible and production-readiness changes are recorded here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
DNA uses semantic versions for tagged releases, while scientific JSON contracts
remain independently versioned and are never silently changed in place.

## [Unreleased]

### Added

- Production-oriented CI/CD and supply-chain verification, including strict cargo-shear dependency/source hygiene and an explicit Ubuntu 24.04 runner baseline.
- Crate-wide Rust lint policy in `Cargo.toml` (`missing_docs`, `unreachable_pub`,
  Clippy `pedantic`) and a tuned release profile (thin LTO, one codegen unit).
- Rust/Python source-boundary enforcement.
- Dependency policy, dependency review, CodeQL, OpenSSF Scorecard, fuzzing,
  self-contained explicit-target auditable release bundles with post-strip metadata verification, release SBOMs, and artifact attestations.

### Changed

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

- An indel whose flanking call is unresolved (its strongest channels tie, so it
  has no primary event) no longer aborts `analyze` or a whole `sample`
  operation; the unresolved flank is omitted from public calls, as SRS-VAR-006
  now states.

## [0.1.0] - Unreleased

Initial pre-production development baseline.
