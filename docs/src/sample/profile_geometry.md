# `src/sample/profile_geometry.rs`

## Purpose

Provides threshold-free mathematical geometry over normalized A/C/G/T evidence
profiles.

## Responsibilities

- Compute Gini impurity for one normalized profile.
- Decompose mean-profile heterogeneity into within-read impurity and between-read dispersion.
- Compute Total Variation distance between forward/reverse mean profiles.
- Keep all functions independent of quality, SNR, noisy-region state, alignment score, configuration, and consensus policy.

## Non-responsibilities

No thresholds, labels, confidence, consensus, genotype, heteroplasmy, or indel interpretation.

## Status

Implemented.
