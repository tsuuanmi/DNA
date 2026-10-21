# CLI

Owns command-line syntax and typed arguments for `analyze`, `basecall`, `evidence`, `signal`,
`sample`, `consensus`, `review`, `run-context`, `call`, and `notation`.
`run-context` takes a batch's read list as exactly one of `--traces` and
`--evidence`; `consensus` and `review` take `--run-context`.

`analyze` takes a trace or, with `--evidence`, the read's evidence document,
exactly one of them (a `clap` argument group); `sample`, `consensus`, and
`review` take the sample's traces or, with `--evidence`, their evidence
documents, likewise.

The CLI is an outer adapter. It parses arguments only; `lib.rs` translates the
parsed values into operation inputs before calling the pipeline. Pipeline and
scientific modules must not depend on `clap` argument structs.

Filesystem/content validation, orchestration, and scientific algorithms belong
downstream.

See [input requirements](../../docs/requirements/input.md),
[configuration contract](../../docs/reference/configuration.md), and
[root usage](../../README.md).
