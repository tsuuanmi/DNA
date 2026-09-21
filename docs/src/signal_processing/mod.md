# `src/signal_processing/mod.rs`

## Purpose

Provides the observation-only signal-processing stage.

## Responsibilities

- Accept the immutable decoded chromatogram, signal-derived calls, and strict signal-processing configuration.
- Sequence rolling feature calculation, basecall-independent locus-evidence extraction, whole-trace integrity assessment, and candidate-noisy region merging.
- Return one typed `DNAAnalysis` without mutating calls or channels.

## Non-responsibilities

No ABIF decoding, channel mutation or smoothing, waveform baseline correction, peak
selection, quality scoring, alignment, variant filtering, reporting, logging, or
file I/O.

## Key function

- `analyze(trace, calls, config) -> Result<DNAAnalysis>`.

## Dependencies

- `locus_evidence`, `integrity`, `features`, `regions`, and `statistics` for algorithms.
- `config`, `error`, and `model` for boundary types.

## Status

Implemented.
