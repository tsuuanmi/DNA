# `src/pipeline/sample.rs`

## Purpose

Runs one multi-read sample-evidence operation.

## Responsibilities

- Load one validated sample identifier, one shared reference/configuration, and one or more traces.
- Process every trace through shared `pipeline::sample_reads`, which delegates each trace to the authoritative `pipeline::observation` path and is reused by validation export.
- Aggregate completed observations through `sample::aggregate` using the strict sample-reconciliation overlap policy.
- Delegate `sample_aggregation_completed` metric extraction to `pipeline::sample_metrics`, including profile availability, contribution states, support mass, mean-profile availability, threshold-free profile heterogeneity, directional profile distance, signal-context counts, and locus/variant topology. These remain evidence summaries rather than consensus or confidence verdicts.
- Build and atomically publish one `dna.sample_evidence/v7` document at `results/<sample-id>.sample.json`.
- Write exactly one persistent `$DNA_LOG_DIR/<sample-id>.log`; all nested per-trace read/basecall/signal/QC/alignment/variant records use that same logger between `sample_read_started` and `sample_read_completed`, with no per-trace log opened by the sample command.

## Non-responsibilities

No F/R pairing, metadata-derived placement, consensus calling, sample-level variant adjudication, or per-read compatibility publication.

## Status

Implemented.
