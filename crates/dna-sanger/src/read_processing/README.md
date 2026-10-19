# Read Processing

Owns shared reference-free Sanger read processing and the Sanger adapter to the
modality-neutral core.

- `mod.rs` runs basecalling, signal processing, callability, and quality
  control for one chromatogram. It decides whether the read can be analyzed at
  all (rejected reads: the signal-window length check before any stage, too few callable calls
  after callability) and returns the read products with their warning counts,
  each call's second-peak label among them ([second_cause](../second_cause/README.md))
  and its phase mixtures ([decompose](../decompose/README.md)).
  Basecall and Variant Analysis both reuse it.
- `evidence.rs` is the Sanger modality adapter
  ([ADR-0069](../../../../docs/decisions/adr/0069-plugin-first-modality-core-post-calling.md)).
  It builds the per-read `ReadEvidence` that alignment and variant calling
  consume:
  - the trim interval;
  - per call: the primary base, the locus evidence profile, the mask, and the
    support vetoes;
  - the Sanger support vetoes: `peak_below_minimum`,
    `relative_quality_not_above_threshold`, and `mixed_supporting_dna`
    (SNVs only);
  - the second-peak assertions (SRS-CALL-022, ADR-0080): whether each
    call's second base stands out, whether its second peak is unexplained,
    and whether a call is `reference_only`: an isolated mixed call whose
    peak is a neighbour's tail or a dye blob, or a single-channel call
    withheld only for its quality or at a benign signal edge (the untrimmed
    start, a weak mask clear of repeat runs);
  - the read's phase mixtures (SRS-CALL-023, ADR-0082), as `decompose`
    asserted them;
  - the Sanger mask labels: `post_homopolymer`, or else the phase-state label.
    Dephased calls become anchoring masks; every other masked call is
    unresolved.

  - the signal-edge verdict of every call (SRS-CALL-013), the leading calls
    vouched for by their base only (SRS-CALL-014; start noise, never a call
    masked after a repeat run), the edge calls whose own
    signal takes their base off the edge (SRS-CALL-017, and at the read's
    start against the start's own norms, SRS-CALL-018; using the peak widths
    signal processing measures), and every call's second base (SRS-CALL-015,
    computed by [second_cause](../second_cause/README.md)'s `second_base`).

  This file is the only place where Sanger peaks, relative quality,
  qualifying channels, and phase states turn into core evidence.

`config.rs` owns the `[sanger_evidence]` section and `SangerConfig`, which
validates every Sanger-owned section together, including the cross-section
check of the relative-quality threshold. `mod.rs` also declares the Sanger
plugin descriptor. Every read-edge parameter, the trim's context margin
included, is the plugin's own configuration; no core key is read.

This module does not use reference context and does not decide variant
eligibility: it supplies labels, and the core reports them.

See [read evidence](../../../dna-kernel/src/read_evidence.rs),
[variant-calling method](../../../../docs/design/variant-calling.md), and
[callability method](../../../../docs/design/callability.md).
