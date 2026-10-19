---
implements: [SRS-OUT]
covers: [src/operation_log.rs, src/pipeline/, src/report/]
---

# Reference and Production Contracts

This directory is the canonical reference for public, machine-visible, configuration, coordinate, schema, and result-contract semantics.

## Human-readable contracts

- [Configuration](configuration.md): strict TOML and environment behavior.
- [Target profiles](profiles.md): target knowledge referenced by the configuration.
- [Basecall result](basecalls.md): `dna.basecalls/v3`.
- [Read evidence](read-evidence.md): `dna.read_evidence/v1`, the modality → core boundary, also the input of `dna analyze --evidence` and of `dna sample`, `dna consensus`, and `dna review` with `--evidence`.
- [Analysis result](analysis/README.md): `dna.analysis/v10`.
- [Sample evidence result](sample-evidence/README.md): `dna.sample_evidence/v10`.
- [Variants result](variants.md): `dna.variants/v1`.
- [Notation result](notation.md): `dna.notation/v1`.
- [Consensus result](consensus.md): `dna.consensus/v1` and its FASTA.
- [Review list](review.md): `dna.review/v1`.
- [Run signal](signal.md): `dna.signal/v1`.
- [Coordinate conventions](coordinates.md): shared coordinate domains and interval semantics.
- [Rust public API](rust-api.md): capability-oriented typed library contract.
  - [Sample analysis API](rust-api-sample.md)
- [Glossary](glossary.md): terms used across the contracts.

## Machine-readable contracts

- [Schema index](schemas/README.md)
  - [Analysis JSON Schema](schemas/analysis-v10.schema.json)
- [Basecall JSON Schema](schemas/basecalls-v3.schema.json)
- [Read evidence JSON Schema](schemas/read-evidence-v1.schema.json)
- [Sample evidence JSON Schema](schemas/sample-evidence-v10.schema.json)
- [Variants JSON Schema](schemas/variants-v1.schema.json)
- [Notation JSON Schema](schemas/notation-v1.schema.json)
- [Consensus JSON Schema](schemas/consensus-v1.schema.json)
- [Review JSON Schema](schemas/review-v1.schema.json)
- [Run signal JSON Schema](schemas/signal-v1.schema.json)
- [Evaluation run JSON Schema](schemas/eval-run-v1.schema.json): `dna.eval_run/v1`, the evaluation harness output, described in [validation](../validation/evaluation-harness.md#outputs)
- [Synthetic examples](examples/README.md)

A schema is authoritative for the exact serialized shape of its named version. Human contract documentation defines semantics and interpretation boundaries that JSON Schema cannot express alone.

Detailed mechanisms live in [design](../design/README.md); normative behavioral requirements live in [requirements](../requirements/README.md). Reference documentation explains stable interfaces and shapes rather than implementation internals.
