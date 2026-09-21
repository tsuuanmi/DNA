# `src/pipeline/input.rs`

## Purpose

Validates command-specific path cardinality and filesystem types, then loads each
declared input exactly once.

## Responsibilities

- Require every declared AB1 and reference input to be a non-empty regular file.
- Load the strict configuration and require its source file to exist.
- Derive and validate output targets for publishing commands; sample scientific input loading remains independent from the production result target so validation tooling can reuse the same loaded traces/reference/config.
- Load one chromatogram for `analyze`/`basecall`, or one or more chromatograms
  for `sample`; load the reference only for reference-guided commands.
- Validate the sample identifier used only for deterministic sample result/log
  naming.

## Non-responsibilities

No directory scanning, globs, manifests, scientific placement, reconciliation, or
algorithm execution.

## Key types and functions

- `AnalysisInputs`: loaded config, chromatogram, reference, and analysis target.
- `BasecallInputs`: loaded config, chromatogram, and reference-free target.
- `SampleInputs`: loaded config, one or more chromatograms, and shared reference.
- `sample_output(sample_id)`: derives and validates the production sample-evidence target.
- `load_analysis(args)`, `load_basecall(args)`, and `load_sample(args)`:
  command-specific entry points.
- `trace_stem(trace) -> Result<&str>`: validates the UTF-8 stem used by
  single-read result/log paths.
- `validate_sample_id(sample_id)`: accepts 1..=128 ASCII characters beginning
  with an alphanumeric and otherwise limited to alphanumeric, `_`, `.`, or
  `-`.

## Invariants and errors

- Required files must be non-empty regular files.
- Output targets must not already exist.
- A missing output parent is allowed and created only at publication time.
- Sample IDs affect naming/provenance only; they never influence placement.
- No manifest/list compatibility input is accepted by the core CLI.

## Dependencies

- `cli` for `AnalyzeArgs`, `BasecallArgs`, and `SampleArgs`.
- `config`, `reference`, `trace`, validated model types, and `error`.

## Biological semantics

None; this module validates and loads declared inputs only.

## Tests

CLI/integration tests exercise command cardinality, invalid paths, invalid sample
IDs, and no-overwrite behavior.

## Status

Implemented.
