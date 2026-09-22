# Release Operations

This is the operational companion to ADR-0018.

A release record must identify the exact:

- DNA semantic version and source revision;
- release tag and supported target triple;
- Rust/Cargo toolchain;
- `Cargo.lock` identity;
- artifact SHA-256;
- SPDX SBOM;
- GitHub artifact provenance and SBOM attestations;
- dependency-policy/audit status;
- CodeQL and CI status;
- fuzz/adversarial status;
- synthetic regression status;
- approved real-AB1 validation status;
- documented runtime and peak-memory measurement.

Use [the release evidence template](release-evidence-template.md) for any release candidate intended to carry the production-ready label.

## Release states

Use explicit language:

- **build verified** — engineering, dependency, security, provenance, and artifact gates passed;
- **scientifically validated for the documented corpus/domain** — approved real-trace evidence passed;
- **production-ready release** — the full ADR-0018 release contract is satisfied.

Do not collapse these into one status.

## Automated delivery

Tags matching `v*` trigger `.github/workflows/release.yml`.

The workflow refuses a tag that does not match `Cargo.toml` or whose commit is not reachable from `main`. It reruns the Rust and dependency gates, builds using the pinned release toolchain, locked dependency graph, and pinned `cargo-auditable`, explicitly preserves `.dep-v0` during stripping, verifies the packaged binary with RustSec, generates an SPDX SBOM from that auditable binary, packages the supported Linux artifact, writes checksums, generates cryptographic GitHub/Sigstore attestations, and publishes the archive/SBOM/checksums as release assets.

The release job does not restore shared CI build caches. This keeps the artifact-producing path independent of cache state written by ordinary CI.

Consumers can verify a downloaded artifact with GitHub CLI attestation verification in addition to checking `SHA256SUMS`.

## Artifact behavior

Existing versioned scientific JSON contracts must not gain unversioned build metadata. Software/build provenance belongs in release metadata/attestations or in a separately versioned scientific output contract.

## Production decision

Automated delivery establishes a **build verified** artifact. It does not by itself establish biological validity. A release must not be described as production-ready until the exact revision also has the scientific, adversarial, performance, and approved real-corpus evidence required by ADR-0018.


## Release immutability

Repository settings must enable GitHub immutable releases before production delivery. The workflow is create-only and does not replace existing release assets. Once published, GitHub must prevent the associated release tag and assets from being modified.


## Runtime and memory evidence

The Linux validation tooling includes `tools/python/scripts/measure_command.py`.
Use it around the exact approved validation command to record wall time, CPU time,
peak RSS, exit status, platform identity, and command arguments as
`dna.performance/v1` JSON. Retain that file (or its checksum and controlled
location) with the release evidence record.


## Target identity

The release target is explicit: `x86_64-unknown-linux-gnu`. The workflow installs that Rust target and builds with `--target`; it does not infer the supported artifact identity from whatever host target the current runner happens to use.


## Self-contained bundle

The Linux release archive contains the executable plus the authoritative default
`config/dna.toml` and `references/rCRS.fasta`. Both runtime assets carry
SHA-256 sidecars inside the archive, and the workflow verifies them before
publishing. A user extracting the archive therefore does not need a source
checkout merely to satisfy DNA's default configuration path or use the bundled
rCRS reference.


## Pre-tag packaging smoke

Pull requests execute the same package-building path before a tag is created. This catches target, strip/audit metadata, runtime-asset, SBOM, archive, and checksum defects without creating attestations or publishing a release.

The tag workflow then adds two privilege-separated stages: a source-free attestation job and a source-free publication job.
