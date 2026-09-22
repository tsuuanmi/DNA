# Production Contracts

This directory owns public, machine-visible, and configuration contracts.

## Human-readable contracts

- [Configuration](configuration.md): strict TOML and environment behavior.
- [Basecall result](basecalls.md): `dna.basecalls/v2`.
- [Analysis result](analysis/README.md): `dna.analysis/v7`.
- [Sample evidence result](sample-evidence/README.md): `dna.sample_evidence/v8`.
- [Coordinate conventions](coordinates.md): shared coordinate domains and interval semantics.

## Machine-readable contracts

- [Schema index](schemas/README.md)
  - [Analysis JSON Schema](schemas/analysis-v7.schema.json)
- [Basecall JSON Schema](schemas/basecalls-v2.schema.json)
- [Sample evidence JSON Schema](schemas/sample-evidence-v8.schema.json)
- [Synthetic examples](examples/README.md)

A schema is authoritative for the exact serialized shape of its named version. Human contract documentation defines semantics and interpretation boundaries that JSON Schema cannot express alone.
