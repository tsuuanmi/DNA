# `src/model/locus_evidence.rs`

## Purpose

Defines immutable basecall-independent evidence at one PLOC-defined locus.

## Responsibilities

- Represent raw co-located A/C/G/T channel values, local baseline/noise, corrected amplitudes, and SNR.
- Preserve the source PLOC window and the independently refined event sample.
- Represent an optional normalized `EvidenceProfile` whose A/C/G/T weights come directly from corrected signal mass.
- Provide channel complementation for reverse-orientation alignment without changing profile mass.

## Non-responsibilities

No primary/IUPAC calling, threshold membership, quality scoring, reference placement, variant interpretation, genotype, heteroplasmy, or reporting. The model only provides the immutable evidence and strand-complement operation consumed by alignment.

## Invariants

A/C/G/T arrays follow `Nucleotide::ALL`. Corrected amplitudes and SNR are non-negative. A profile exists only when corrected signal mass is positive and then sums to one apart from floating-point roundoff. Complementation maps `[A,C,G,T]` to `[T,G,C,A]` exactly.

## Status

Implemented.
