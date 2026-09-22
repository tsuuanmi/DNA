# Changelog

All notable user-visible and production-readiness changes are recorded here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
DNA uses semantic versions for tagged releases, while scientific JSON contracts
remain independently versioned and are never silently changed in place.

## [Unreleased]

### Added

- Production-oriented CI/CD and supply-chain verification, including strict cargo-shear dependency/source hygiene.
- Rust/Python source-boundary enforcement.
- Dependency policy, dependency review, CodeQL, OpenSSF Scorecard, fuzzing,
  auditable release binaries with post-strip metadata verification, release SBOMs, and artifact attestations.

### Changed

- Corrected the minimum supported Rust version to match language features used
  by the codebase.

## [0.1.0] - Unreleased

Initial pre-production development baseline.
