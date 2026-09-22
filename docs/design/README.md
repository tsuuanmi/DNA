# Design and Methods

These documents describe **current mechanisms**: how DNA implements its scientific and technical behavior. They are living design documentation, not public serialization contracts and not historical decision records.

- [Pipeline](pipeline.md): orchestration and stage boundaries.
- [ABIF decoding](abif-decoding.md): bounded container decode and scientific-tag extraction.
- [Basecalling](basecalling.md): PLOC-window peak selection and primary/ambiguity calls.
- [Signal processing](signal-processing/README.md): rolling SNR, locus evidence, trace-integrity observations, and interpretation limits.
- [Quality control](quality-control.md): penalty, relative quality, best-section selection, and end trimming.
- [Alignment](alignment.md): profile-aware Gotoh placement, traceback, orientation, and circular mapping.
- [Variant calling](variant-calling.md): difference extraction, normalization, eligibility, and deterministic ordering.
- [Sample evidence](sample-evidence/README.md): aggregation after independent read placement.

Normative behavior lives in [requirements](../requirements/README.md), rationale in [decisions](../decisions/README.md), exact public/configuration shapes in [reference](../reference/README.md), and exploratory work in [research](../research/README.md).
