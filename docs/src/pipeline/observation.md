# `src/pipeline/observation.rs`

## Purpose

Owns the one authoritative reference-guided scientific path from a decoded trace to `ReadObservation`.

## Responsibilities

Run shared read processing, pass post-trim `EvidenceProfile` observations into profile-aware reference alignment, run configured normalized primary-sequence variant calling, stage logging, exclusion diagnostics, and warning aggregation including PLOC/vendor mismatch and exact clipping counts. Both single-read analysis and sample analysis consume this path.

## Non-responsibilities

No input loading, output publication, JSON projection, sample aggregation, or sample interpretation.

## Architectural role

This module removes duplicated read science between command pipelines. It is also the dependency boundary where immutable signal evidence becomes input to reference placement without allowing alignment to mutate upstream calls/profiles. `pipeline::analyze` and `pipeline::sample` differ only in orchestration after one-read observation creation.

## Status

Implemented.
