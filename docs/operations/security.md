# Security and Trust Boundaries

DNA is a local scientific CLI, but it processes untrusted binary/text input and sensitive biological data.

## Untrusted inputs

Treat as untrusted:

- ABIF/AB1 bytes;
- FASTA bytes;
- TOML configuration;
- filenames and filesystem metadata;
- external manifests/scripts inputs.

Required protections include bounded file sizes, checked arithmetic, checked slices, validated cardinality, allocation caps, and typed failures.

## Rust boundary

First-party production code forbids unsafe Rust. A future unsafe/FFI boundary requires an explicit ADR and targeted safety validation.

## Resource exhaustion

Input-controlled dimensions must be bounded before large allocation. Alignment-cell limits, reference limits, AB1 size caps, and indel limits are correctness/security controls, not just performance settings.

## Filesystem

- result publication is atomic and no-overwrite;
- cleanup scripts validate targets before deletion;
- symlink/path-escape behavior must be explicit;
- scientific failures must not leave a partial result presented as successful.

## Data privacy

Sequences, chromatograms, variant evidence, and sample-linked metadata may be identifying biological data. Logs and test fixtures must follow [data policy](../data.md) and must not leak complete sequences or dense signal payloads unless the contract explicitly requires them.

## Dependencies

Release dependencies are reviewed for known advisories, source policy, licenses, and explicit exceptions as defined by ADR-0018.


## Static and adversarial analysis

- CodeQL runs the Rust `security-extended` query suite on pull requests, `main`, and weekly.
- The ABIF parser has a dedicated `cargo-fuzz` harness with pull-request smoke campaigns and longer scheduled campaigns.
- A fuzz-discovered defect is incomplete until a minimized regression input is retained when practical.

## Supply chain

- GitHub Actions are pinned by immutable full commit SHA.
- Dependabot proposes reviewed updates for Cargo, uv tooling, Actions, and the Rust toolchain.
- `cargo-deny` enforces trusted dependency sources, license policy, bans, yanked/advisory policy, and wildcard rules.
- `cargo-audit` independently checks RustSec.
- dependency review blocks vulnerable dependency additions before merge.
- OpenSSF Scorecard provides a scheduled external posture check.

## Release integrity

Release artifacts include SHA-256 checksums and an SPDX SBOM. GitHub artifact attestations bind the published archive and SBOM to the repository, workflow, source revision, and OIDC-backed signing identity.

The repository-level disclosure process is defined in [SECURITY.md](../../SECURITY.md), and the settings required to prevent CI bypass are defined in [repository governance](repository-governance.md).
