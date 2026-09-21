# `src/report/basecall.rs`

## Purpose

Builds the typed `dna.basecalls/v2` document without filesystem effects.

## Responsibilities

- Consume one completed config, trace, base-call, signal, and quality-control set.
- Derive the ambiguity sequence from ordered call records.
- Validate sequence lengths, trim bounds, and retained-primary slice consistency.
- Project provenance, read sequences, trace-integrity evidence, merged noisy regions, and public warnings.

## Non-responsibilities

No scientific processing, logging, serialization, publication, reference,
alignment, variant, or compatibility output.

## Dependencies

`config`, `error`, basecall/shared result models, scientific read models, and the
shared `report::signal` projection.

## Tests

`tests/basecall.rs` and `scripts/validate_result_schemas.py` exercise this builder.
