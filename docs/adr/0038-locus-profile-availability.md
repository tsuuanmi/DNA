# ADR-0038: Expose differential-locus profile availability before consensus

- **Status:** Accepted
- **Date:** 2026-09-19

## Context

DNA now preserves one reference-oriented `CallDNAEvidence` for every
call-backed differential-locus observation. The optional `EvidenceProfile`
inside that object is absent only when the source locus has no positive corrected
signal mass.

A future profile-aware sample consensus cannot consume a nucleotide profile where
no profile exists. Before defining any quality/noise weighting or contributor
eligibility policy, the sample model therefore needs an explicit denominator of
which retained locus observations actually have profile evidence.

Tracy pairwise consensus directly adds two trace-profile vectors at overlapping
columns before normalizing. That demonstrates profile combination, but it gives
the two traces equal profile mass and does not provide a calibrated local
weighting model that DNA should copy as final consensus policy.

## Decision

Extend internal differential-locus support topology with:

~~~text
profile_reads
profile_forward_reads
profile_reverse_reads
~~~

A read counts as profile-available exactly when its retained call-backed
observation contains `CallDNAEvidence.profile = Some(...)`.

The topology MUST satisfy:

~~~text
profile_reads = profile_forward_reads + profile_reverse_reads
profile_reads <= reads - deletion_reads
~~~

Unresolved primary calls may still count when their basecall-independent profile
exists. A valid zero-signal call with no profile does not count. Deletions never
count as nucleotide-profile observations.

Profile-availability counts were initially internal. ADR-0053 later includes those same counts in the public differential-locus support topology in `dna.sample_evidence/v8` without changing contributor eligibility.

## Consequences

- Future sample interpretation has an explicit local nucleotide-profile denominator.
- Forward/reverse profile availability is visible without implying bidirectional agreement.
- Missing profile evidence remains distinct from deletion and from noisy-region context.
- Production logging consumes the topology directly instead of rescanning observations.
- No profile vectors are summed and no quality/noise/read weight is assigned.

## Non-goals

This decision does not define local contribution eligibility, equal-read profile
voting, profile summation, consensus thresholds, confidence, artifact
classification, sample-level variant adjudication, genotype, or heteroplasmy.
