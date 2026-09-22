# Versioning

DNA versions externally meaningful contracts explicitly.

## Public JSON contracts

Closed schemas such as `dna.analysis/v7`, `dna.basecalls/v2`, and `dna.sample_evidence/v8` are immutable by version. Incompatible shape or semantic changes require a new contract version.

## Configuration

Configuration is strict and versioned by its documented schema/version field. Unknown or unsupported values fail rather than silently falling back.

## Software releases

Software versions identify a source revision and reproducible build/release evidence. A software-version change does not silently redefine an existing JSON schema version.

## ADRs and proposals

ADR/proposal identifiers are stable historical identifiers and are never reused. A changed decision is represented through supersession; an evolved proposal updates status/implementation metadata without pretending prior review history did not occur.

## Documentation

Living current-state documents do not need artificial version suffixes. Git history records previous revisions; obsolete current-state copies are deleted rather than kept as `old`, `legacy`, or `v2-final` documentation.
