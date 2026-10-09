# Read Processing

Owns shared reference-free Sanger read processing and the Sanger adapter to the
modality-neutral core.

- `mod.rs` runs basecalling, signal processing, callability, and quality
  control for one chromatogram. It decides whether the read can be analyzed at
  all (rejected reads) and returns the read products with their warning counts.
  Basecall and Variant Analysis both reuse it.
- `evidence.rs` is the Sanger modality adapter
  ([ADR-0069](../../docs/decisions/adr/0069-plugin-first-modality-core-post-calling.md)).
  It builds the per-read `ReadEvidence` that alignment and variant calling
  consume:
  - the trim interval;
  - per call: the primary base, the locus evidence profile, the mask, and the
    support vetoes;
  - the Sanger support vetoes: `peak_below_minimum`,
    `relative_quality_not_above_threshold`, and `mixed_supporting_dna`
    (SNVs only);
  - the Sanger mask labels: `post_homopolymer`, or else the phase-state label.
    Dephased calls become anchoring masks; every other masked call is
    unresolved.

  This file is the only place where Sanger peaks, relative quality,
  qualifying channels, and phase states turn into core evidence.

`config.rs` owns the `[sanger_evidence]` section and `SangerConfig`, which
validates every Sanger-owned section together, including the cross-section
check of the relative-quality threshold. `mod.rs` also declares the Sanger
plugin descriptor. The trim's context margin comes from the composer, which
passes the core's `read_end_margin`.

This module does not use reference context and does not decide variant
eligibility: it supplies labels, and the core reports them.

See [read evidence](../read_evidence.rs),
[variant-calling method](../../docs/design/variant-calling.md), and
[callability method](../../docs/design/callability.md).
