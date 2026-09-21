# `src/model/mod.rs`

## Purpose

Owns the validated domain types shared across scientific stages.

## Responsibilities

- Declare the model submodules: `alignment`, `basecall_result`, `basecalls`, `coordinate`, `locus_evidence`, `nucleotide`, `quality`, `read_observation`, `reference`, `result`, `sample_evidence`, `sample_result`, `signal`, `trace`, and `variant`.
- Provide validated domain vocabulary with no I/O, CLI, configuration loading, or
  algorithm dependencies.

## Non-responsibilities

No filesystem access, argument parsing, logging, or algorithm execution.

## Key types and functions

- `coordinate`: zero-based to one-based conversion.
- `nucleotide`: canonical bases and reverse complement.
- `trace`: four-channel `Chromatogram` plus vendor evidence.
- `basecalls`: primary/ambiguous calls with trace positions.
- `locus_evidence`: basecall-independent per-locus A/C/G/T evidence and normalized profile.
- `signal`: ordered locus evidence, rolling SNR windows, and merged candidate-noisy regions.
- `quality`: quality vector, trim bounds, and QC result.
- `read_observation`: complete one-read scientific products after
  evidence-derived placement.
- `reference`: name, sequence, topology, and checksums.
- `alignment`: selected orientation, score, reference segments, metrics, and
  per-column coordinates without final gapped-row duplication.
- `variant`: normalized `Variant`, `VariantKind`, alleles, mappings, and
  exclusion diagnostics without report-only labels.
- `result`: compact `AnalysisResult` plus result records shared by both
  contracts.
- `basecall_result`: reference-free `BasecallResult` matching basecalls schema v2.
- `sample_evidence`: compact internal read-registry, differential-locus, and normalized-variant evidence.
- `sample_result`: public `dna.sample_evidence/v8` records.

## Invariants and errors

Constructors and loaders enforce aligned vector lengths, valid nucleotide
symbols, explicit coordinate systems, bounded qualities, and consistent internal
records. Model files have no I/O, CLI, configuration loading, or algorithm
dependencies.

## Dependencies

- `serde` for serialization where records are emitted.
- `error` for typed failures in the few fallible conversions.

## Apollo mapping

Replaces loosely coupled `Trace`, `BaseCalls`, `ReferenceSlice`, alignment
arrays, and `Variant` structs with validated Rust types.

## Requirements and decisions

ADR-0002, ADR-0028; `SRS-IN-003`, `SRS-CFG-005`, `SRS-BC-005`, `SRS-SIG-008`, `SRS-SIG-009`, `SRS-ALN-005`,
`SRS-VAR-005`, and `SRS-OUT-002`.

## Tests

Focused invariant tests live beside each child implementation; the integration
tests exercise the full model through the pipeline.

## Status

Implemented.
