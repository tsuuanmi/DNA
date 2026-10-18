# dna-sanger

Owns the Sanger modality plugin
([ADR-0069](../../../docs/decisions/adr/0069-plugin-first-modality-core-post-calling.md)):
everything from ABIF bytes to one read's `ReadEvidence`. It depends only on
`dna-kernel`.

- [abif](abif/README.md) — ABIF decoding into a chromatogram; the `fuzzing` feature exposes the container parser.
- `locus.rs` — Sanger locus window geometry.
- [basecalling](basecalling/README.md), [signal_processing](signal_processing/README.md), [callability](callability/README.md), and [quality_control](quality_control/README.md) — the reference-free stages.
- [read_processing](read_processing/README.md) — stage orchestration, `SangerConfig`, the Sanger plugin descriptor, and the `ReadEvidence` adapter.
- [run_signal](run_signal/README.md) — the second-channel signal a run's reads share (PROP-0009).
- [model](model/README.md) — canonical Sanger evidence.
- [error](error/README.md) — the Sanger stage failures, `SangerError`.

See [crates](../../README.md).
