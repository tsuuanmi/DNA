---
decided-by: [ADR-0018, ADR-0069]
---

# Versioning

DNA versions externally meaningful contracts explicitly.

## Public JSON contracts

Closed schemas such as `dna.analysis/v10`, `dna.basecalls/v3`, and `dna.sample_evidence/v10` are immutable by version. Incompatible shape or semantic changes require a new contract version.

A contract version introduced in the current unreleased cycle may be revised in place until the first tagged release that emits it; every such revision is recorded as **Breaking** under `[Unreleased]` in `CHANGELOG.md`.

## Evaluation contracts

The evaluation harness's output contract, `dna.eval_run/v1`, is versioned
like a result contract, but no file of a run carries the identifier: a run's
`run.json` records the harness revision, whose schema describes the run. The
harness is never released, and runs from any revision are read back as
baselines and by `dna-eval power`, so the unreleased-cycle exception does not
apply: a change to a file's columns, formats, or meaning, or to which files a
run writes, is a new version.

## Configuration

Configuration is strict and versioned by its documented schema/version field. Unknown or unsupported values fail rather than silently falling back. The same unreleased-cycle rule applies: a configuration schema version introduced in the current unreleased cycle may be revised in place before the first tagged release that accepts it, each revision recorded as Breaking.

## Plugins

Each plugin of the static registry
([[ADR-0069]])
carries a method version, recorded in result provenance. A change that can
alter a plugin's output for the same inputs and configuration increases that
plugin's version in the same change, and the CHANGELOG records it. A pure
refactor with byte-identical output leaves it unchanged. Plugin versions start
at 1 and are independent of the software version and of contract versions.

## Software releases

Tagged binary releases use semantic versioning. A release tag is exactly `v<version from Cargo.toml>`.

While the software version is below 1.0, minor releases may contain intentional breaking software-interface changes, but those changes must still be documented in `CHANGELOG.md` and must not silently mutate an existing scientific JSON contract.

Software versions identify a source revision and reproducible build/release evidence. A software-version change does not silently redefine an existing JSON schema version.

Only releases whose exact revision satisfies ADR-0018 and has a completed release-evidence record may be described as production-ready. The latest source on `main` is a development line, not a production release.

## ADRs and proposals

ADR/proposal identifiers are stable historical identifiers and are never reused. A changed decision is represented through supersession; an evolved proposal updates status/implementation metadata without pretending prior review history did not occur.

## Documentation

Living current-state documents do not need artificial version suffixes. Git history records previous revisions; obsolete current-state copies are deleted rather than kept as `old`, `legacy`, or `v2-final` documentation.
