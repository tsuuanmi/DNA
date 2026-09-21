# ADR-0042: Decompose nucleotide profile heterogeneity

- **Status:** Accepted
- **Date:** 2026-09-19

## Context

A mean nucleotide profile cannot distinguish reproducible mixed signal within each
chromatogram from disagreement between otherwise clean chromatograms. For
example, two pure A/G-disagreeing reads and two identical 50/50 A/G-mixed reads
have the same arithmetic mean profile.

DNA needs threshold-free evidence that preserves this distinction before any
consensus, confidence, or heteroplasmy interpretation.

## Decision

For every non-empty eligible-profile partition, DNA derives normalized profile
heterogeneity using Gini impurity.

For one normalized profile p:

~~~text
impurity(p) = 1 - sum(p_i^2)
~~~

For n contributors with arithmetic mean profile mu:

~~~text
within_profile_impurity = mean(impurity(p_r))
total_profile_heterogeneity = impurity(mu)
between_profile_dispersion = total_profile_heterogeneity - within_profile_impurity
~~~

Mathematically:

~~~text
total_profile_heterogeneity = within_profile_impurity + between_profile_dispersion
~~~

The decomposition is retained independently for total, forward, and reverse
eligible-profile partitions.

These are geometry summaries only. They do not define mixed-signal thresholds,
biological mixture, heteroplasmy, agreement classes, confidence, or consensus.

## Consequences

- Reproducible within-read mixture can be separated from inter-read disagreement
  even when the mean profile is identical.
- Forward and reverse contributor groups retain their own internal geometry.
- No new quality/SNR/noise weighting or public schema field is introduced.

## Non-goals

No thresholds, categorical states, genotype/heteroplasmy inference, consensus
base, calibrated confidence, or indel/gap interpretation.
