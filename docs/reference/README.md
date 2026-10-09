# Reference and Production Contracts

This directory is the canonical reference for public, machine-visible, configuration, coordinate, schema, and result-contract semantics.

## Human-readable contracts

- [Configuration](configuration.md): strict TOML and environment behavior.
- [Target profiles](profiles.md): target knowledge referenced by the configuration.
- [Basecall result](basecalls.md): `dna.basecalls/v3`.
- [Analysis result](analysis/README.md): `dna.analysis/v9`.
- [Sample evidence result](sample-evidence/README.md): `dna.sample_evidence/v10`.
- [Coordinate conventions](coordinates.md): shared coordinate domains and interval semantics.
- [Rust public API](rust-api.md): capability-oriented typed library contract.

## Machine-readable contracts

- [Schema index](schemas/README.md)
  - [Analysis JSON Schema](schemas/analysis-v9.schema.json)
- [Basecall JSON Schema](schemas/basecalls-v3.schema.json)
- [Sample evidence JSON Schema](schemas/sample-evidence-v10.schema.json)
- [Synthetic examples](examples/README.md)

A schema is authoritative for the exact serialized shape of its named version. Human contract documentation defines semantics and interpretation boundaries that JSON Schema cannot express alone.

Detailed mechanisms live in [design](../design/README.md); normative behavioral requirements live in [requirements](../requirements/README.md). Reference documentation explains stable interfaces and shapes rather than implementation internals.
