# Read Processing

Owns shared reference-free Sanger read processing and the Sanger adapter to the
modality-neutral core.

- `mod.rs` runs basecalling, signal processing, callability, and quality
  control for one chromatogram. It decides whether the read can be analyzed at
  all (rejected reads: the signal-window length check before any stage, too few callable calls
  after callability; `admit` returns either as an `Admission` with the evidence
  its stages produced, as the sample commands, `dna evidence`, and the run
  signal take it) and returns the read products with their warning counts,
  each call's second-peak label among them ([second_cause](../second_cause/README.md))
  and its phase mixtures, decomposed by
  [dna-trace](../../../dna-trace/src/README.md) (ADR-0084) from the read's
  calls, locus profiles, callability, and trim start, and its tract decodes
  ([tract](../tract/README.md), SRS-CALL-025).
  Basecall and Variant Analysis both reuse it.
- `evidence.rs` is the Sanger modality adapter
  ([ADR-0069](../../../../docs/decisions/adr/0069-plugin-first-modality-core-post-calling.md)).
  It builds the per-read `ReadEvidence` that alignment and variant calling
  consume:
  - the trim interval;
  - per call: the primary base, the locus evidence profile, the mask, and the
    support vetoes;
  - the Sanger support vetoes: `peak_below_minimum`,
    `relative_quality_not_above_threshold`, `mixed_supporting_dna`
    (SNVs only), and `signal_spike`;
  - the second-peak assertions (SRS-CALL-022, ADR-0080): whether each
    call's second base stands out, whether its second peak is unexplained,
    and whether a call is `reference_only`: a mixed call whose peak is a
    neighbour's tail or a dye blob, between neighbours with one qualifying
    channel or a tail or blob second peak of their own, withheld otherwise at most for
    its quality, off any signal edge but a benign one; or a single-channel
    call without an unexplained second peak, withheld only for its quality
    or at a benign signal edge (the untrimmed start, a weak mask clear of
    repeat runs; never the read's physical end). `reference_only_calls`
    holds the one benign-edge test both use;
  - the read's phase mixtures (SRS-CALL-023, ADR-0082), as `dna-trace`
    asserted them, and its tract decodes (SRS-CALL-025), as
    [tract](../tract/README.md) turned them into `TractLengths`;
  - the calls beside a junction whose spacing marks a missing peak
    (callability's `gaps`, SRS-CALL-026), vouched for by their bases only:
    they place no indel and show no run's end;
  - the calls of every electrophoretic spike signal processing finds
    (SRS-SIG-012), withheld for their base and for indels with
    `signal_spike` after every vouching and edge rule (SRS-CALL-027): no
    run's end, no reference confirmation;
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
  qualifying channels, and phase states turn into core evidence. `reason`
  maps a published label back onto the plugin's own reason, so a recorded
  `dna.read_evidence/v1` document is read back without new labels
  (PROP-0011).

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
