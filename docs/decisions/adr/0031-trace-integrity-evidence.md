# ADR-0031: Preserve PLOC and signal-integrity evidence without artifact reclassification

- **Status:** Accepted
- **Date:** 2026-09-19

## Context

DNA's current scientific method is explicitly anchored at ABI `PLOC.2`
positions. Tracy's public issue history exposes two evidence-foundation failure
modes that should be visible before adding more sample interpretation:

1. some AB1 files contain vendor base/quality series whose cardinality does not
   match the usable PLOC series;
2. high-amplitude or clipped channel artifacts can dominate a chromatogram even
   when downstream profile methods appear numerically confident.

Before this decision DNA rejected any PBAS/PCON length mismatch with PLOC.
That behavior conflated optional vendor metadata consistency with the validity of
the PLOC-defined scientific locus series. DNA also retained rich per-locus
corrected amplitudes internally without exposing concise whole-trace integrity
evidence.

## Decision

DNA adds one observation-only `TraceIntegrity` record to
`DNAAnalysis`. It is carried through `ReadObservation`, reference-free
basecall results, single-read analysis results, and each sample read summary.

### PLOC remains mandatory

`PLOC.2` remains a hard dependency of the current re-calling method. Decode
still requires a non-empty, strictly increasing, in-range PLOC series.

DNA does **not** add a PLOC-independent event detector or fallback.

### Optional vendor cardinality mismatch becomes evidence

`PBAS.2` and `PCON.2` remain optional vendor evidence. When present, their
decoded lengths no longer have to equal the PLOC count.

DNA processes exactly the PLOC-defined loci. Downstream vendor comparisons use
bounded indexing, so a shorter optional vendor series simply provides no vendor
observation beyond its available prefix; a longer vendor series does not create
additional DNA calls.

The integrity record retains:

- `ploc_count`;
- optional `vendor_primary_count`;
- optional `vendor_quality_count`.

Each present vendor series whose count differs from `ploc_count` contributes
one `ploc_vendor_length_mismatches` warning.

### PLOC spacing observations

When at least two PLOC positions exist, DNA records the minimum, median, and
maximum adjacent PLOC spacing in trace-sample units. A one-locus trace has no
spacing summary.

These are observations only. No new spacing threshold or automatic rejection is
introduced.

### Exact clipping observation

ABIF analyzed channels are decoded from signed 16-bit DATA values. DNA counts
channel samples exactly equal to `i16::MIN` or `i16::MAX` as
`clipped_channel_samples`.

This is a factual representation-limit observation, not a claim that every
clipped sample is a dye blob or that every artifact clips.

A nonzero count contributes to the public warning summary and operational warning
total. It does not alter calls, quality, alignment, or variants.

### Relative event-signal scale

For every PLOC locus, DNA already refines an event sample and derives
non-negative corrected A/C/G/T amplitudes. Define:

~~~text
event_signal(locus) = sum(corrected_amplitudes[A,C,G,T])
~~~

Across loci with positive event signal:

~~~text
maximum_to_median_event_signal_ratio
    = max(event_signal) / median(event_signal)
~~~

The emitted ratio is rounded with the existing six-decimal signal metric rule.
If no locus has positive event signal, the ratio is absent.

The ratio is intentionally **not** thresholded into an artifact flag in this ADR.
A large ratio is validation/review evidence until an artifact classifier is
separately specified and validated.

### Public contracts

The change is intentionally incompatible:

- `dna.analysis/v6` -> `dna.analysis/v7`;
- `dna.basecalls/v1` -> `dna.basecalls/v2`;
- `dna.sample_evidence/v4` -> `dna.sample_evidence/v5`.

Single-read contracts expose integrity under `signal_quality.integrity`.
Sample evidence stores the same concise integrity summary once per read.

Earlier schema/example files are removed; no compatibility result or alias is
retained.

## Consequences

- DNA can analyze a valid PLOC-defined prefix even when optional vendor series
  are longer or shorter.
- The number of scientific calls remains exactly the PLOC count.
- Trace-integrity evidence survives into sample reconciliation for future
  artifact-aware admission/consensus work.
- Exact clipping is visible without pretending to solve all high-amplitude
  artifact detection.
- High-amplitude event imbalance is preserved quantitatively without an
  unvalidated biological or artifact label.
- Future PLOC-independent event detection, dye-blob classification, baseline
  shift detection, and neighbor-interference classification remain separate
  methods.
