# Callability

Owns signal-derived per-read callability: per-position features, phase-state
segments, the typed mask, and the callable span (`dna.read_callability/v1`).

Key children:

- `features.rs` — dominance, secondary ratio, slippage shift offset, spacing
  deviation, and weakness per position;
- `runs.rs` — homopolymer and dinucleotide repeat runs in the read's own
  primary calls, used only as priors for the window after each run;
- `phase.rs` — rolling defect fractions, repeat priors, and the hysteresis
  state machine that decides where the read is callable;
- `classify.rs` — the phase state of every masked segment and its repeat
  attribution;
- `mask.rs` — the per-position mask, the callable span, and their invariants;
- `sanger.rs` — the Sanger adapter that builds the core's evidence records from
  the chromatogram loci, locus evidence, and primary calls.

Dependency rule: `features`, `runs`, `phase`, `classify`, and `mask` form a
modality-generic core that depends only on the plain records in
`model::callability` (plus the shared metric rounding of `signal_processing`);
they never import `Chromatogram`, `BaseCalls`, `LocusEvidence`, or any other
Sanger type. `sanger.rs` is the only Sanger-aware file and owns no thresholds.
A future modality supplies its own `PositionEvidence` records and reuses the
core unchanged; no trait is introduced until a second adapter exists.

The stage runs after signal processing and before quality control, reference-,
profile-, and filename-free, in trace order for both strands. It never modifies
channels, loci, calls, selected peaks, or locus evidence; a masked position keeps
its original call and evidence. While the mask is observation-only it alters no
trim bound, alignment, warning total, or variant eligibility.

See [callability requirements](../../docs/requirements/callability.md),
[callability method](../../docs/design/callability.md), ADR-0067, and the
[evidence invariants](../../docs/architecture/invariants/evidence.md).
