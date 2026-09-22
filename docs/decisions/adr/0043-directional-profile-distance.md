# ADR-0043: Measure directional nucleotide-profile distance

- **Status:** Accepted
- **Date:** 2026-09-19

## Context

When both selected read orientations contribute eligible nucleotide profiles,
DNA already retains independent forward and reverse mean profiles.

A threshold-free orientation comparison is useful before any directional
agreement/discordance policy is defined.

## Decision

When both orientation-specific mean profiles exist, DNA derives Total
Variation distance:

~~~text
directional_profile_distance = 0.5 * sum(abs(forward_mean_i - reverse_mean_i))
~~~

The value is bounded to [0, 1] for normalized A/C/G/T profiles:

- 0 means identical orientation-specific profile distributions.
- 1 means disjoint profile mass.

The metric is absent unless both orientations have eligible profile evidence.

The field is named directional_profile_distance, not conflict, discordance, or
confidence. No threshold is attached.

## Consequences

- Direction-specific disagreement is measurable without collapsing profiles to
  primary basecalls.
- Same-direction biological or technical replicates remain distinct from
  forward/reverse comparison.
- The metric does not imply independent biological replication.

## Non-goals

No discordance threshold, strand-bias test, winning nucleotide, confidence,
consensus, heteroplasmy, or quality weighting.
