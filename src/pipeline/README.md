# Pipeline

Owns production operation orchestration for `analyze`, `basecall`, and
`sample`.

Pipeline entry points accept operation values such as trace/reference paths,
sample identifiers, and an explicit configuration path rather than CLI/`clap`
argument structs or process environment. This keeps orchestration independent of
the command-line frontend and prepares the same boundary for non-CLI callers.

Key children own CLI/application filesystem naming, overwrite protection,
logging/publication, sample-read orchestration, and sample metrics. Scientific
Sanger source loading is owned by `input::sanger`; shared reference-free read
processing lives in `read_processing.rs`; reference-guided read observation is
owned by `variant_analysis`.

Scientific input adapters are independent of deterministic CLI publication
targets. Pipeline validates output naming and overwrite protection separately,
so `results/*` state is not part of the scientific input contract.

Pipeline passes its file-backed logger through the internal `StageLog`
capability to reusable scientific modules. It owns log construction, terminal
error persistence, synchronization, JSON projection, and publication, but not
the shared scientific implementation.

Pipeline code composes capabilities and preserves typed failures; algorithm
internals remain in their owning modules.

See [pipeline method](../../docs/design/pipeline.md),
[pipeline invariants](../../docs/architecture/invariants/pipeline.md), and
[output requirements](../../docs/requirements/output.md).
