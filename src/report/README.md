# Report

Owns projection of completed scientific state into versioned public contracts,
deterministic JSON serialization, and atomic no-overwrite publication. It also
holds the review's labelling rules: the `review` post-calling plugin
(`review.rs`, with `length.rs` for run lengths) is the one part of this module
that judges evidence rather than projecting it.

Output contract: every public JSON contract
([schemas](../../docs/reference/schemas/README.md)); the
[step contracts](../../docs/architecture/data-flow.md#step-contracts) map each
to its step and command. The review list,
[`dna.review/v1`](../../docs/reference/schemas/review-v1.schema.json), is
derived here.

Key children separate analysis (`analysis.rs`, with `aligned_calls.rs`, each aligned call's reference column and the core's trust in it, `dna.analysis/v10`), basecall, read evidence (`evidence.rs`, `dna.read_evidence/v1`), sample, run signal (`run_signal.rs`, `dna.signal/v1`), run context (`run_context.rs`, `dna.run_context/v1`), variant,
`per_base_decimal` notation rendering (`notation.rs`), what every document
shares (`json.rs`: the plugin and profile provenance projections and
serialization), and atomic publication concerns. Every document's provenance records the workflow's
plugins; analysis, sample, variants, notation, consensus, and review provenance
also record the target-profile identity.

Core records carry no Sanger evidence
([ADR-0069](../../docs/decisions/adr/0069-plugin-first-modality-core-post-calling.md)).
`call.rs` projects core-only calls into `dna.variants/v1`, `consensus.rs`
projects a sample consensus into `dna.consensus/v1` and its FASTA (also
rendered from the consensus alone, with the same record names),
`review.rs` derives a sample's review items and labels from its analysis
and call (the `review` plugin it declares, whose version is the labelling
rules') and projects `dna.review/v1`; `length.rs` judges the runs of the
profile's anchored repeat windows from the reads' tract decodes and in-phase
readings for its `length_unconfirmed`, `length_not_recorded`, and
`tract_length` labels, a profile length tract's run noted when its dominant
length is determined (SRS-API-020),
and `notation.rs` renders the per-read notation once (`calls`, with each call's
reference and called base and the reads containing it) and projects it for sample,
call, and notation documents and builds `dna.notation/v1`. Every document reads a read's Sanger evidence through `model::sanger_read`
(`SangerRead` and `Rejection`, projected there from the Sanger products, so
the commands from read evidence report the same values), and
`sanger_call.rs` joins a call's reference-strand peaks and
relative quality, and
`sample.rs` and `consensus.rs` join each read's `SangerRead` (integrity,
callability, quality, noisy regions) and each rejected read's `Rejection` by
read identity and call index. `evidence.rs` writes a read's evidence, or the
rejection of a read the Sanger plugin set aside.

Calls and decisions remain upstream: the review rules label evidence the
analysis and the call already hold, read the `SangerRead` only for each
call's scan position, relative quality, second peak with its cause, and
window peaks (which the run signal counts), and change no call.

See [contracts](../../docs/reference/README.md),
[output requirements](../../docs/requirements/output.md), and
[output/operations invariants](../../docs/architecture/invariants/output-operations.md).
