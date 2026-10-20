# dna-sanger

Owns the Sanger modality plugin
([ADR-0069](../../../docs/decisions/adr/0069-plugin-first-modality-core-post-calling.md)):
everything from ABIF bytes to one read's `ReadEvidence`. It depends on
`dna-kernel` and on [dna-trace](../../dna-trace/src/README.md), the trace
signal library it calls for phase mixtures and tract decodes (ADR-0084).

Output contract: [`dna.basecalls/v3`](../../../docs/reference/schemas/basecalls-v3.schema.json) (`dna basecall`),
[`dna.read_evidence/v1`](../../../docs/reference/schemas/read-evidence-v1.schema.json) (`dna evidence`), and
[`dna.signal/v1`](../../../docs/reference/schemas/signal-v1.schema.json) (`dna signal`), serialized by the facade
([step contracts](../../../docs/architecture/data-flow.md#step-contracts)).

- [abif](abif/README.md) — ABIF decoding into a chromatogram; the `fuzzing` feature exposes the container parser.
- `locus.rs` — Sanger locus window geometry.
- [basecalling](basecalling/README.md), [signal_processing](signal_processing/README.md), [callability](callability/README.md), and [quality_control](quality_control/README.md) — the reference-free stages.
- [read_processing](read_processing/README.md) — stage orchestration, `SangerConfig`, the Sanger plugin descriptor, and the `ReadEvidence` adapter.
- [run_signal](run_signal/README.md) — the second-channel signal a run's reads share (PROP-0009).
- [second_cause](second_cause/README.md) — why a call shows a second peak (SRS-CALL-021); labels only.
- [tract](tract/README.md) — the read's homopolymer tract decodes (SRS-CALL-025) and the `[trace]` section; evidence only.
- [model](model/README.md) — canonical Sanger evidence.
- [error](error/README.md) — the Sanger stage failures, `SangerError`.

See [crates](../../README.md).
