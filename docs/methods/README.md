# Current Scientific and Algorithmic Methods

These documents describe current production algorithms and scientific behavior, not public serialization shape.

- [Pipeline](pipeline.md): canonical orchestration map and stage boundaries.
- [ABIF decoding](abif-decoding.md): bounded container decode and canonical scientific-tag extraction.
- [Basecalling](basecalling.md): PLOC-window peak selection and primary/ambiguity calls.
- [Signal processing](signal-processing.md): rolling SNR, locus evidence, trace-integrity observations, and interpretation limits.
- [Quality control](quality-control.md): penalty, relative quality, best-section selection, and end trimming.
- [Alignment](alignment.md): profile-aware Gotoh placement, traceback, orientation, and circular mapping.
- [Variant calling](variant-calling.md): difference extraction, normalization, eligibility, and deterministic ordering.
- [Sample evidence aggregation](sample-evidence.md): aggregation after independent read placement.

Public configuration and result semantics live under [contracts](../contracts/README.md). Exploratory methods belong under [research](../research/README.md) until explicitly promoted into the production SRS/ADR/method/contract system.
