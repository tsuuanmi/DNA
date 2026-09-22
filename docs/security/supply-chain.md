# Supply-Chain and Vulnerability Management

Dependency, CI, and release integrity are part of the production trust boundary.

## Dependencies

- Rust and Python tooling use committed lockfiles.
- `cargo-shear --deny-warnings` rejects unused/misplaced dependencies and unlinked Rust source.
- `cargo-deny` enforces license, source, ban, yank, wildcard, and advisory policy.
- pinned `cargo-audit` independently checks RustSec.
- GitHub dependency review blocks moderate-or-higher vulnerable additions on pull requests.
- Dependabot maintains Cargo, uv-tooling, GitHub Actions, and Rust-toolchain update PRs.
- dependency/advisory findings are reviewed; they are not suppressed merely to make CI green.

## GitHub Actions

Third-party Actions are pinned to full 40-character commit SHAs. Repository policy rejects mutable pins and privileged untrusted triggers; actionlint and zizmor analyze workflow syntax and security.

Linux verification/release jobs pin Ubuntu 24.04 rather than following the moving `ubuntu-latest` label.

## Release provenance

The supported release target is explicit: `x86_64-unknown-linux-gnu`.

Pull requests exercise the real release-package path. Tagged delivery separates three trust zones:

1. package/build with read-only repository access;
2. source-free OIDC/SBOM/provenance attestation;
3. source-free publication with release-write permission.

The package path:

- builds with pinned `cargo-auditable` and locked dependencies;
- strips while preserving `.dep-v0`, then verifies and audits the packaged binary;
- bundles authoritative `config/dna.toml` and `references/rCRS.fasta` with checksums;
- emits an SPDX JSON SBOM and SHA-256 checksums;
- verifies the archive before staging it for attestation/publication.

Release builds do not restore shared CI build caches.

## Vulnerabilities

A newly discovered vulnerability is triaged by:

1. affected dependency/component and reachable code path;
2. impact on confidentiality, integrity, availability, or scientific correctness;
3. fixed-version or mitigation availability;
4. required source/config/contract changes;
5. regression/security validation before release.

If a mitigation changes architecture, public behavior, or trust boundaries, update the corresponding decision and canonical docs rather than maintaining a hidden workaround.

See [dependency policy](../engineering/dependencies.md), [release engineering](../engineering/release.md), [repository governance](../governance/repository.md), and [ADR-0018](../decisions/adr/0018-production-readiness-release-contract.md).
