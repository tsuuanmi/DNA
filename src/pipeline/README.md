# Pipeline

Owns production operation orchestration for `analyze`, `basecall`, and
`sample`.

Pipeline entry points accept operation values such as trace/reference paths,
sample identifiers, and an explicit configuration path rather than CLI/`clap`
argument structs or process environment. This keeps orchestration independent of
the command-line frontend and prepares the same boundary for non-CLI callers.

Key children separate input loading, shared read processing, reference-guided
observation processing, sample-read processing, metrics, and result
publication.

Pipeline code sequences stages and preserves typed failures; algorithm internals
remain in their owning modules.

See [pipeline method](../../docs/design/pipeline.md),
[pipeline invariants](../../docs/architecture/invariants/pipeline.md), and
[output requirements](../../docs/requirements/output.md).
