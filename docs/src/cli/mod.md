# `src/cli/mod.rs`

## Purpose

Defines the stable command-line syntax and typed arguments for single-read and sample operations.

## Responsibilities

- Expose `signal analyze <trace.ab1> --reference <fasta>`, reference-free `signal basecall <trace.ab1>`, and `signal sample <sample-id> <trace.ab1>... --reference <fasta>`.
- Derive outputs as `results/<trace-stem>.json`, `results/<trace-stem>.basecalls.json`, or `results/<sample-id>.sample.json`; no output-path option is accepted.

## Inputs and outputs

Converts operating-system arguments into `Cli`, `Command`, and `AnalyzeArgs`.
Path existence and file content are validated later by input modules.

## Key types and functions

- `Cli`: the top-level parser.
- `Command::Analyze(AnalyzeArgs)`: reference-guided analysis.
- `Command::Basecall(BasecallArgs)`: reference-free read processing.
- `AnalyzeArgs`: one positional AB1 path and one required `--reference` FASTA path.
- `BasecallArgs`: one positional AB1 path and no reference or output option.
- `SampleArgs`: one sample ID, one or more positional AB1 paths, and one required `--reference` FASTA path.

## Invariants and errors

Argument syntax failures are handled by Clap. The CLI does not run algorithms.
`analyze` and `basecall` accept exactly one AB1. `sample` accepts one or more AB1 paths. Directories, manifests, globs, and list-file compatibility inputs are rejected. `DNA_CONFIG` is process configuration, not a CLI trace input.

## Dependencies

- `clap` for `Parser`, `Args`, and `Subcommand`.

## Apollo mapping

DNA uses one focused command rather than Apollo's six-command surface.

## Requirements and decisions

ADR-0001, ADR-0003, ADR-0007; `SRS-IN-001`, `SRS-IN-010` through
`SRS-IN-012`, and `SRS-CFG-004`.

## Tests

`tests/cli.rs`.

## Status

Implemented.
