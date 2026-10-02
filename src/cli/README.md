# CLI

Owns command-line syntax and typed arguments for `analyze`, `basecall`, and
`sample`.

The CLI is an outer adapter. It parses arguments only; `lib.rs` translates the
parsed values into operation inputs before calling the pipeline. Pipeline and
scientific modules must not depend on `clap` argument structs.

Filesystem/content validation, orchestration, and scientific algorithms belong
downstream.

See [input requirements](../../docs/requirements/input.md),
[configuration contract](../../docs/reference/configuration.md), and
[root usage](../../README.md).
