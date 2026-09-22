# ADR-0039: Define structural nucleotide contribution eligibility

- **Status:** Accepted
- **Date:** 2026-09-19

## Context

ADR-0037 unified call-backed sample signal evidence and ADR-0038 separated
profile availability from future contribution policy.

The next boundary must answer whether a retained locus observation can
participate in future nucleotide aggregation without prematurely turning
uncalibrated relative quality, per-channel SNR, or candidate-noisy-region
membership into a rejection rule.

Profile availability alone is necessary but is not the long-term definition of
eligibility. DNA therefore needs an explicit policy state now, before future
quality/artifact-aware refinements.

## Decision

Every retained differential-locus observation has one internal
`NucleotideContribution` state:

~~~text
Eligible
MissingProfile
DeletionEvent
~~~

The first policy is intentionally structural:

~~~text
deletion
    -> DeletionEvent

call-backed observation + CallDNAEvidence.profile
    -> Eligible

call-backed observation + missing profile
    -> MissingProfile
~~~

An unresolved call with a real basecall-independent profile remains
`Eligible`. Its collapsed character does not erase the underlying nucleotide
evidence.

Candidate-noisy-region membership, corrected amplitudes, per-channel SNR, and
uncalibrated relative quality remain preserved observations but do not change
this first eligibility state.

Profile availability remains the independent denominator introduced by
ADR-0038. The current policy therefore admits every profile-bearing observation,
but availability and eligibility remain separate concepts so future validated
policy can evolve without redefining evidence availability.

Deletion observations remain on the separate gap/event path and are never
converted into nucleotide contributors.

The current public `dna.sample_evidence/v7` contract remains unchanged.
Production logging consumes all three contribution states.

## Consequences

- Future nucleotide consensus has an explicit contributor gate rather than
  implicitly equating all coverage with support.
- Zero-signal loci do not gain fabricated nucleotide evidence.
- Mixed/unresolved profile evidence survives into future aggregation.
- Deletions remain distinct from nucleotide support.
- No new threshold, weight, vote, confidence, genotype, or heteroplasmy
  semantics are introduced.

## Non-goals

This decision does not define quality-aware weighting, noisy-tail suppression,
artifact-aware rejection, per-channel SNR thresholds, gap confidence, consensus
states, or sample-level variant calling.
