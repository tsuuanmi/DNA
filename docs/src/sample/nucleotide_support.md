# `src/sample/nucleotide_support.rs`

## Purpose

Aggregates eligible reference-oriented nucleotide profiles at one retained
differential locus without making a consensus decision.

## Responsibilities

- Consume only observations classified as `NucleotideContribution::Eligible`.
- Require every eligible observation to retain a real `EvidenceProfile`.
- Add each normalized profile with unit read mass.
- Preserve total, forward-only, and reverse-only A/C/G/T support vectors plus
  contributor counts.
- Derive total support from the forward/reverse partitions.
- Derive an arithmetic mean `EvidenceProfile` for each non-empty total/forward/reverse partition by dividing support by contributor count.
- Accumulate per-read profile impurity and delegate threshold-free heterogeneity decomposition plus forward/reverse Total Variation to `profile_geometry`.

## Non-responsibilities

No relative-quality, SNR, amplitude, noisy-region, or alignment-score weighting.
No normalization into a consensus profile, base call, confidence class, variant
verdict, genotype, or heteroplasmy inference.

## Invariants

Each eligible read contributes exactly one normalized profile at the locus.
Missing-profile and deletion observations contribute no nucleotide mass.
Forward plus reverse contributor counts equal total contributors, and total
support is the channel-wise sum of forward and reverse support.

## Status

Implemented.
