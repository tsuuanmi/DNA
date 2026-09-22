# Platform Support

Platform support is evidence-based.

## Release-supported

| Target | Status | Evidence |
|---|---|---|
| `x86_64-unknown-linux-gnu` | release-supported | Explicit `--target` self-contained release bundle plus CI verification, runtime-asset checksums, SBOM, checksums, and attestations |

A release-supported target is one that the release process builds and verifies for
the exact release revision. The current Linux build/verification runner baseline is Ubuntu 24.04 and changes to that baseline are reviewed explicitly.

## Buildable but not release-supported

Other targets may compile or run successfully, but DNA does not claim production
support for them until they are included in automated build/test/release evidence.

Adding a supported target requires:

- CI compilation and tests on the target or an equivalent justified runner;
- target-specific filesystem/atomic-publication validation;
- documented runtime and resource evidence;
- tagged release artifacts and checksums;
- provenance/SBOM attestation;
- scientific regression evidence showing no platform-specific semantic drift.

Do not infer support from Rust's target list or from a successful local build.
