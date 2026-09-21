# `src/pipeline/sample_metrics.rs`

## Purpose

Derives operational summary metrics from completed internal sample evidence.

## Responsibilities

- Centralize all `sample_aggregation_completed` counters and numeric summaries.
- Consume profile availability, contribution state, support topology, profile geometry, signal-context counts, and variant-call summaries.
- Keep operational metric extraction out of sample orchestration.

## Non-responsibilities

No scientific evidence mutation, filtering, thresholds, report projection, or consensus logic.

## Status

Implemented.
