# Pipeline

Owns production command orchestration for `analyze`, `basecall`, and `sample`.

Key children separate input loading, shared read processing, reference-guided
observation processing, sample-read processing, metrics, and command publication.

Pipeline code sequences stages and preserves typed failures; algorithm internals
remain in their owning modules.

See [pipeline method](../../docs/methods/pipeline.md),
[pipeline invariants](../../docs/architecture/invariants/pipeline.md), and
[output requirements](../../docs/srs/output.md).
