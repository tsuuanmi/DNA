# Pipeline

Owns production operation orchestration for `analyze`, `basecall`, `sample`,
`call`, and `notation`.

Pipeline entry points accept operation values such as trace/reference paths,
sample identifiers, and an explicit configuration path rather than CLI/`clap`
argument structs or process environment. This keeps orchestration independent of
the command-line frontend and prepares the same boundary for non-CLI callers.

Key children own CLI/application filesystem naming, overwrite protection,
logging/publication, sample-read orchestration, sample metrics, and per-read
profile representation for the notation view of sample, call, and notation
documents (`represent.rs`, which composes `variant_normalization` and
`variant_nomenclature` when the target profile declares notation). Scientific
Sanger source loading is owned by `input::sanger`; shared reference-free read
processing lives in `read_processing` (`prepare` through callability, then
`finish`); reference-guided read observation is owned by `variant_analysis`.
`sample_reads.rs` processes each trace in those two steps so that a read with
too few callable calls is recorded as rejected while the others continue.
`notation.rs` runs the post-calling plugins over one variants document:
representation, conformance, and `dna.notation/v1`. `call.rs` runs the core
alone over reviewed consensus sequences: the sequence
adapter, `read_call`, optional notation, and `dna.variants/v1`.
`sample.rs` splits every read into its `CalledRead`, which the neutral `sample`
aggregation consumes, and its Sanger attachment, which `sample_metrics.rs` and
the sample report join by read identity (ADR-0069).

Scientific input adapters are independent of deterministic CLI publication
targets. Pipeline validates output naming and overwrite protection separately,
so `results/*` state is not part of the scientific input contract.

Each command runs as one `Operation`: it opens the operation's `OperationLog`,
records the start event (failing fast if that write fails), runs the stages in
the log's `tracing` scope, synchronizes the log before publishing, and on failure
records the terminal event with the last-entered stage span. Pipeline owns log
construction, terminal error persistence, synchronization, JSON projection, and
publication, but not the shared scientific implementation.

Pipeline code composes capabilities and preserves typed failures; algorithm
internals remain in their owning modules.

See [pipeline method](../../docs/design/pipeline.md),
[pipeline invariants](../../docs/architecture/invariants/pipeline.md), and
[output requirements](../../docs/requirements/output.md).
