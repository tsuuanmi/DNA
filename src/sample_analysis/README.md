# Sample Analysis

Owns the sample placement pass and the public capability built on it.

`reads.rs` runs every trace of a sample once through the one-read path of
`variant_analysis`: an admitted read keeps its `CalledRead` and Sanger
attachment; a read shorter than one signal window (checked before any stage),
one with too few callable calls, or one the core cannot place becomes a
`RejectedRead`, an outcome like an admitted read rather than a failure. Reads run on up to the caller's number of scoped
threads; each read's error or panic is kept as its outcome, workers stop after
the lowest failed read, and outcomes are taken in trace order up to the first
failure, as a sequential pass would. A `ReadScope` decides where each read's
`tracing` instrumentation goes: the CLI operation log captures it per read and
appends it in trace order, and the public call sends it to the caller's
dispatcher inside the caller's span.

`mod.rs` aggregates the placed reads into the core sample evidence, builds the
consensus in the evidence's read order, and exposes `analyze_sanger`, which
returns both as a typed `SampleAnalysis` from one pass. The `sample` and
`consensus` pipelines call the same pass, aggregation, and consensus, then add
their own logging, notation, projection, and publication. `mod.rs` also keeps
the loaded configuration, profile, and reference in the result, exposes the
profile's reportable regions (`SampleAnalysis::regions`), keeps the placed
reads so that `SampleAnalysis::observations` gives every read's observation of
each covered position of an interval as the core records it at a differential
locus (`dna_sample::observations`),
returns the rejected
reads instead of an error when none is admitted (`analyze_sanger_outcome`, which
`analyze_sanger` wraps; its body over loaded inputs, `analyze`, also serves the
`review` operation with the operation log as the read scope), and defines the
public `ConsensusCall` and `NotationCall` types of `SampleAnalysis::call`, which
`pipeline` implements because it composes the `consensus` report and the `call`
pipeline.

It does not open logs, derive output paths, serialize JSON, or publish files.

See [Rust API contract](../../docs/reference/rust-api.md#sample-analysis),
[pipeline method](../../docs/design/pipeline.md#sample-placement-pass), and
[sample evidence requirements](../../docs/requirements/sample-evidence/README.md).
