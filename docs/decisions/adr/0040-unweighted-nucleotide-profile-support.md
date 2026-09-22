# ADR-0040: Accumulate eligible profiles with unit read mass

- **Status:** Accepted
- **Date:** 2026-09-19

## Context

ADR-0038 made profile availability explicit and ADR-0039 introduced a separate
nucleotide-contribution eligibility policy.

The next sample-consensus foundation needs a deterministic way to preserve the
combined nucleotide evidence at one retained locus before any quality-aware
weighting or consensus threshold is defined.

Raw corrected amplitude and per-channel SNR are not suitable default weights:
their scale is not calibrated across traces, instruments, runs, or local signal
contexts. Reusing those values now would silently convert observation features
into confidence.

Tracy demonstrates that normalized trace profiles can be combined by direct
vector addition. DNA can retain that useful intermediate evidence without
adopting it as a final consensus rule.

## Decision

For each retained differential locus, DNA derives one internal
`LocusNucleotideSupport` from observations whose
`NucleotideContribution == Eligible`.

Each eligible observation contributes its normalized reference-oriented
`EvidenceProfile` with unit read mass:

~~~text
support_A = sum(profile_A)
support_C = sum(profile_C)
support_G = sum(profile_G)
support_T = sum(profile_T)
~~~

The same accumulation is retained independently for forward and reverse selected
read orientations.

The object records:

~~~text
contributors
forward_contributors
reverse_contributors

support[A,C,G,T]
forward_support[A,C,G,T]
reverse_support[A,C,G,T]
~~~

with:

~~~text
contributors = forward_contributors + reverse_contributors
support = forward_support + reverse_support
~~~

Missing-profile observations and deletion events contribute no nucleotide mass.

This is explicitly an **unweighted evidence accumulator**, not a consensus
profile. DNA does not normalize the accumulated vector into a call and does
not select a winning nucleotide.

The current public `dna.sample_evidence/v7` contract remains unchanged.
Production logging consumes aggregate contributor topology and total unweighted
profile mass without logging nucleotide-specific composition.

## Consequences

- Cross-read mixed signal survives as continuous A/C/G/T evidence instead of
  collapsing to primary-character majority vote.
- Forward and reverse support remain inspectable separately.
- Unit read mass gives a deterministic baseline while avoiding unvalidated
  cross-trace amplitude/SNR weighting.
- A future validated weighting policy can operate at the contribution boundary
  without redefining profile availability or eligibility.
- No consensus base, confidence, sample variant, genotype, or heteroplasmy
  semantics are introduced.

## Non-goals

This decision does not claim equal reads are biologically independent or equally
reliable. It does not define quality/SNR/artifact weights, normalization,
dominance thresholds, conflict classes, gap/deletion weighting, or final
consensus.
