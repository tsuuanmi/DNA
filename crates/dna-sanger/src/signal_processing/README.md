# Signal Processing

Owns observation-only signal analysis over immutable decoded chromatogram
channels.

It coordinates rolling features, basecall-independent locus evidence,
Sanger-integrity observations, and candidate-noisy region merging.

`mod.rs` owns the one length rule (SRS-SIG-007): `check_length` fails a read
with fewer calls than one window. Read processing applies it before
basecalling, so the feature and locus-context code assumes a full window and
keeps no length guard of its own.

This module does not mutate channels, perform reference interpretation, or decide
variant eligibility.

`statistics.rs` owns the baseline, noise, and SNR primitives and `SpanEstimates`,
which holds one trace and estimates the channel statistics of each of its
sample spans once: the rolling feature windows (`features.rs`) and the locus
contexts (`locus_evidence.rs`) span the same runs of locus windows, share the
estimates, and read the trace only through that record.

`config.rs` owns the `[signal_processing]` configuration section: its raw record, validated
record, and rules (ADR-0069).

See [signal requirements](../../../../docs/requirements/signal-processing.md),
[signal methods](../../../../docs/design/signal-processing/README.md), and
[evidence invariants](../../../../docs/architecture/invariants/evidence.md).
