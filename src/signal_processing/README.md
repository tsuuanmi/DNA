# Signal Processing

Owns observation-only signal analysis over immutable decoded chromatogram
channels.

It coordinates rolling features, basecall-independent locus evidence,
trace-integrity observations, and candidate-noisy region merging.

This module does not mutate channels, perform reference interpretation, or decide
variant eligibility.

See [signal requirements](../../docs/srs/signal-processing.md),
[signal methods](../../docs/methods/signal-processing/README.md), and
[evidence invariants](../../docs/architecture/invariants/evidence.md).
