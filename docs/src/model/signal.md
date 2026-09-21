# `src/model/dna.rs`

## Purpose

Defines the internal observation-only signal-quality result.

## Responsibilities

- Represent one rolling call/sample window with minimum primary SNR, maximum
  secondary SNR, and the candidate-noisy decision. Maximum secondary SNR remains
  internal for the pipeline's operational aggregate and is omitted from compact
  v5.
- Own one immutable `TraceIntegrity` summary containing PLOC/vendor cardinality, PLOC spacing, exact signed-16-bit clipping count, and optional event-signal scale ratio.
- Own the ordered collection of basecall-independent `LocusEvidence` records alongside rolling windows and merged noisy regions.
- Mark candidate-noisy windows without changing any scientific call.
- Represent the deterministic union of overlapping or adjacent candidate-noisy windows.
- Provide aggregate noisy-window and noisy-call counts for operational logging.

## Non-responsibilities

No feature calculation, artifact classification, threshold validation, smoothing, base calling, trimming, variant filtering, serialization, or file I/O.

## Key types

- `DNAWindow`: 0-based half-open call/sample intervals, minimum primary SNR,
  internal maximum secondary SNR, and candidate-noisy flag.
- `LocusEvidence` is defined in `model/locus_evidence.rs` and retained here as the authoritative per-locus signal-evidence layer.
- `NoisyRegion`: merged call/sample intervals and minimum primary SNR.
- `TraceIntegrity`: observation-only whole-trace structural/signal-scale evidence plus vendor-length mismatch count helper.
- `DNAAnalysis`: trace integrity, ordered locus evidence, windows, and merged regions, with
  aggregate count helpers. Isolated candidate windows can remain in `windows`
  without appearing in `noisy_regions`.

## Invariants

All intervals are ordered in original trace orientation. A/C/G/T arrays follow
`Nucleotide::ALL`; SNR values are finite and non-negative; corrected amplitudes
are non-negative. Regions do not overlap and do not bridge clean gaps.

## Dependencies

None outside the standard library.

## Status

Implemented.
