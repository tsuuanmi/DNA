# Versioning

DNA versions externally meaningful contracts explicitly.

## Public JSON contracts

Closed schemas such as `dna.analysis/v7`, `dna.basecalls/v2`, and `dna.sample_evidence/v8` are immutable by version. Incompatible shape or semantic changes require a new contract version.

## Configuration

Configuration is strict and versioned by its documented schema/version field. Unknown or unsupported values fail rather than silently falling back.

## Software releases

Tagged binary releases use semantic versioning. A release tag is exactly `v<version from Cargo.toml>`.

While the software version is below 1.0, minor releases may contain intentional breaking software-interface changes, but those changes must still be documented in `CHANGELOG.md` and must not silently mutate an existing scientific JSON contract.

Software versions identify a source revision and reproducible build/release evidence. A software-version change does not silently redefine an existing JSON schema version.

Only releases whose exact revision satisfies ADR-0018 and has a completed release-evidence record may be described as production-ready. The latest source on `main` is a development line, not a production release.

## ADRs and proposals

ADR/proposal identifiers are stable historical identifiers and are never reused. A changed decision is represented through supersession; an evolved proposal updates status/implementation metadata without pretending prior review history did not occur.

## Documentation

Living current-state documents do not need artificial version suffixes. Git history records previous revisions; obsolete current-state copies are deleted rather than kept as `old`, `legacy`, or `v2-final` documentation.
