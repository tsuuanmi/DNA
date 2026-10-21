# Pipeline

Owns production operation orchestration for `analyze`, `basecall`, `evidence`, `signal`, `sample`,
`consensus`, `review`, `run-context`, `call`, and `notation`, and the in-process sample call.
`run_context.rs` reads a batch's `sample_id<TAB>path` list, analyses each
sample as `consensus` does, sums its cross-strand pairs into the batch's run
context (PROP-0014), and publishes `dna.run_context/v1`; `consensus` and
`review` load one with `--run-context` (`load_run_context`) into the
sample's inputs.

Pipeline entry points accept operation values such as trace/reference paths,
sample identifiers, and an explicit configuration path rather than CLI/`clap`
argument structs or process environment. This keeps orchestration independent of
the command-line frontend and prepares the same boundary for non-CLI callers.

Key children own CLI/application filesystem naming, overwrite protection,
logging/publication, sample metrics, and per-read
profile representation for the notation view of sample, call, and notation
documents (`represent.rs`, which composes `variant_normalization` and
`variant_nomenclature` when the target profile declares notation). Scientific
Sanger source loading is owned by `input::sanger`; shared reference-free read
processing lives in `read_processing` (`prepare` through callability, then
`finish`); reference-guided read observation is owned by `variant_analysis`.
`sample.rs` and `consensus.rs` run the sample placement pass of
`sample_analysis`, which records a read shorter than one signal window
(`ShortRead::of` reads it from the failure of `read_processing::prepare`), one
with too few callable calls, or one the core cannot place
(`dna_call::read_call::PlacementRejection`) as rejected while the others
continue; the operation log implements its `ReadScope`, so records
of concurrently processed reads land in trace order. They run the reads on as
many workers as the available parallelism, which `mod.rs` reads once per
process.
`consensus.rs` then builds the consensus and publishes `dna.consensus/v1` with
its FASTA.
`review.rs` runs the same pass through `sample_analysis::analyze` under the
operation log, calls the consensus with `SampleAnalysis::call` when a read is
admitted, and publishes `dna.review/v1` for every sample: called, every read
rejected, or the call failed; the review projects each read's phase
mixtures (`report::review_mixtures`), given the profile notation of a
component's edits, which `represent.rs` writes. `plugins.rs` composes its
workflow by the post-calling plugins the profile declares.
`review.rs` also implements the public `SampleAnalysis::review_list`
(SRS-API-015): the check items and the noted length tracts of the same review
list, without writing it.
`sample_call.rs` implements the public `SampleAnalysis::call`: it names an
analysis's consensus segments as the consensus FASTA does and runs the
scientific step of `call.rs` for a consensus document (`call_consensus`, which
the `call` command shares: the consensus's decided edits, SRS-VAR-021) with the
analysis's configuration, profile, and reference, without paths, logs, or
publication. It lives here, not in `sample_analysis`,
because it composes `report` and the `call` pipeline, which depend on
`sample_analysis`.
`notation.rs` runs the post-calling plugins over one variants document:
representation, conformance, the EMPOP query, and `dna.notation/v1`; `sample_call.rs`
builds the same query in process. `call.rs` runs `call`: over a consensus document, the consensus's decided
edits (`call_consensus`); over reviewed consensus sequences, the core
alone: the sequence
adapter, `read_call`, optional notation, and `dna.variants/v1`.
`analyze.rs` runs `analyze` from a trace or, with `--evidence`, from the
read's `dna.read_evidence/v1` document (`input::evidence`, SRS-API-018): both
paths place and call the read with `variant_analysis::observation::call` and
build the same document from the core's `CalledRead` and the read's
`SangerRead`, under the trace's stem (`path::evidence_stem` drops the
document's `.evidence`).
`sample.rs`, `consensus.rs`, and `review.rs` take a `SampleInput`: the
sample's traces, or with `--evidence` the reads' `dna.read_evidence/v1`
documents (`input::evidence`, SRS-API-019), which `mod.rs` loads into the
reads of one placement pass. `sample.rs` splits every read into its
`CalledRead`, which the neutral `sample` aggregation consumes, and its
`SangerRead`, which `sample_metrics.rs` and the sample report join by read
identity (ADR-0069).
`evidence.rs` admits the read as the sample commands do
(`read_processing::admit`): a read set aside before the core gets a document
with its `Rejection` instead of its calls.

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
