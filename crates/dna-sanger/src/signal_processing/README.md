# Signal Processing

Owns observation-only signal analysis over immutable decoded chromatogram
channels.

It coordinates rolling features, basecall-independent locus evidence,
Sanger-integrity observations, candidate-noisy region merging, and
electrophoretic spikes.

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

`spike.rs` finds the read's electrophoretic spikes (SRS-SIG-012): four
channels peaking within a call spacing or so, each far above the local
primary peak of the read's own calls, and the calls whose primary peak lies
inside. It is the one signal feature that withholds calls, and only through
the Sanger evidence (SRS-CALL-027, ADR-0088); it reads the basecalls' selected
peaks, never a reference.

`config.rs` owns the `[signal_processing]` configuration section: its raw record, validated
record, and rules (ADR-0069).

See [signal requirements](../../../../docs/requirements/signal-processing.md),
[signal methods](../../../../docs/design/signal-processing/README.md), and
[evidence invariants](../../../../docs/architecture/invariants/evidence.md).
