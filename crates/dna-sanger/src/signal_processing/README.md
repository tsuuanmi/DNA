# Signal Processing

Owns observation-only signal analysis over immutable decoded chromatogram
channels.

It coordinates rolling features, basecall-independent locus evidence,
Sanger-integrity observations, and candidate-noisy region merging.

This module does not mutate channels, perform reference interpretation, or decide
variant eligibility.

`statistics.rs` owns the baseline, noise, and SNR primitives and `SpanEstimates`,
which estimates the channel statistics of each sample span once per read: the
rolling feature windows (`features.rs`) and the locus contexts
(`locus_evidence.rs`) span the same runs of locus windows and share them.

`config.rs` owns the `[signal_processing]` configuration section: its raw record, validated
record, and rules (ADR-0069).

See [signal requirements](../../../../docs/requirements/signal-processing.md),
[signal methods](../../../../docs/design/signal-processing/README.md), and
[evidence invariants](../../../../docs/architecture/invariants/evidence.md).
