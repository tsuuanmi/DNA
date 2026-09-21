# Versioning Policy

DNA versions the software and its scientific output contracts separately.

## Software versions

Tagged binary releases use semantic versioning.

While the software version is below 1.0, minor releases may contain intentional
breaking software-interface changes, but those changes must still be documented
in `CHANGELOG.md` and must not silently mutate an existing scientific JSON
contract.

A release tag is exactly `v<version from Cargo.toml>`.

## Scientific contracts

JSON schema identifiers are explicit independent contracts, for example:

- `dna.basecalls/v2`;
- `dna.analysis/v7`;
- `dna.sample_evidence/v8`.

An incompatible shape or semantic change requires a new contract version. A
software release may therefore change without changing a JSON contract, and a
contract version bump does not imply a particular software major version.

## Production support

Only releases whose exact revision satisfies ADR-0018 and has a completed release
evidence record may be described as production-ready.

The latest source on `main` is a development line, not a production release.
Security support is defined separately in root `SECURITY.md`.
