# `src/report/mod.rs`

## Purpose

Defines typed result projection, shared deterministic serialization, and atomic
no-overwrite publication.

## Responsibilities

- Re-export the analysis, basecall, and sample completed-input wrappers plus their
  focused builders.
- Share deterministic JSON serialization across
  `dna.analysis/v7`, `dna.basecalls/v2`, and
  `dna.sample_evidence/v7`.
- Keep analysis projection, basecall projection, sample projection, shared
  trace-integrity/signal projection, variant projection, and publication in
  separate focused modules.

## Non-responsibilities

No basecalling, signal feature calculation, quality scoring, alignment, variant
calling, sample aggregation, or operational logging.

## Key types and functions

- `build_analysis(completed)`: projects one `ReadObservation` to analysis v7.
- `build_basecall(completed)`: projects reference-free read stages to basecalls
  v2.
- `build_sample(completed)`: projects `SampleEvidence` to sample-evidence v7.
- `serialize(result) -> Result<Vec<u8>>`: deterministic pretty JSON plus one
  trailing newline.
- `publish(path, bytes) -> Result<()>`: atomic no-overwrite publication.
- Child modules: `json`, `basecall`, `sample`, `signal`, `variant`, and
  `atomic`.

## Invariants and errors

Each builder emits exactly its named versioned contract. Reporting performs
projection/consistency validation only; scientific decisions remain upstream.

## Dependencies

Completed model/stage outputs, `error`, and serialization/filesystem support.

## Tests

Integration tests and schema validation cover contract projection; atomic
publication has focused unit coverage.

## Status

Implemented.
