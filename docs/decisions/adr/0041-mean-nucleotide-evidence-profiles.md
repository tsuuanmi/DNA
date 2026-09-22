# ADR-0041: Derive mean nucleotide evidence profiles

- **Status:** Accepted
- **Date:** 2026-09-19

## Context

ADR-0040 accumulates every structurally eligible reference-oriented
`EvidenceProfile` with unit read mass and retains total plus forward/reverse
A/C/G/T support vectors.

Each contributing `EvidenceProfile` is already normalized to unit positive
signal mass. Therefore dividing an accumulated support vector by its contributor
count yields the arithmetic mean profile for that contributor set without
introducing any additional weighting or scientific threshold.

This normalized summary is useful for later directional comparison and consensus
research, but it must not be confused with a consensus call or calibrated
confidence.

## Decision

For each retained differential locus, derive:

~~~text
mean_profile
forward_mean_profile
reverse_mean_profile
~~~

from the existing unit-mass support vectors:

~~~text
mean_profile[channel] =
    support[channel] / contributors

forward_mean_profile[channel] =
    forward_support[channel] / forward_contributors

reverse_mean_profile[channel] =
    reverse_support[channel] / reverse_contributors
~~~

A mean profile is absent when its contributor count is zero.

The arithmetic mean preserves A/C/G/T evidence shape only. DNA does not select
the largest channel, assign a consensus base, threshold mixed evidence, or map
the profile to confidence.

No relative-quality, SNR, corrected-amplitude, noisy-region, or alignment-score
weight is introduced.

The current public `dna.sample_evidence/v7` contract remains unchanged.
Production sample logging consumes mean-profile presence when counting loci with
nucleotide support and bidirectional support.

## Consequences

- Downstream interpretation can compare total and orientation-specific evidence
  on the same normalized A/C/G/T scale.
- Contributor count remains explicit and is not erased by normalization.
- Missing-profile and deletion observations still contribute neither support mass
  nor a mean nucleotide profile.
- Mean profiles are mathematical evidence summaries, not consensus, confidence,
  genotype, or heteroplasmy estimates.

## Non-goals

This decision does not define contribution weighting, directional-discordance
thresholds, winning-nucleotide rules, consensus states, gap/deletion consensus,
sample variants, or calibrated confidence.
