# `src/report/dna.rs`

## Purpose

Owns shared projection from internal signal analysis to compact trace-integrity
and merged-region result records.

## Responsibilities

Consume `DNAAnalysis`, project immutable `TraceIntegrity`, omit rolling
windows, and map merged call/sample intervals plus minimum primary SNR into
`DNAQualityResult`. `project_integrity` is reused by sample read summaries.

## Non-responsibilities

No feature calculation, region merging, policy decisions, serialization, or
publication.

## Dependencies

`model::signal` and shared interval/signal result records.

## Tests

Analysis and basecall integration/schema tests exercise both projection callers.
