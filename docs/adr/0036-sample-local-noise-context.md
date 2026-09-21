# ADR-0036: Preserve local noisy-region context in sample evidence

- **Status:** Accepted
- **Date:** 2026-09-19

## Context

DNA already derives merged candidate-noisy call intervals under
`dna.windowed_snr/v1`. Those intervals are observation-only: a merged region
can contain calls that are not individually weak, so current scientific contracts
explicitly forbid using the region as automatic variant eligibility or artifact
classification.

After sample aggregation, call-backed differential-locus and normalized-variant
evidence retained quality and nucleotide profiles but not whether the original
call fell inside one of those existing noisy regions. Future local contribution
research would otherwise have to reach backward into `ReadObservation` or lose
that context entirely.

## Decision

Sample aggregation preserves existing noisy-region membership without changing
its interpretation.

For every call-backed differential-locus observation:

~~~text
in_noisy_region = Some(true | false)
~~~

For a deletion, where no source nucleotide call exists:

~~~text
in_noisy_region = None
~~~

For every variant-associated source call, sample evidence preserves the same
boolean membership.

Membership is computed only from the existing half-open
`NoisyRegion.call_start_0based..call_end_0based_exclusive` intervals. Sample
code does not recalculate SNR, merge windows, or introduce a threshold.

This context was initially internal. ADR-0053 later projects call-backed differential-locus noisy-region membership into `dna.sample_evidence/v8` without changing its observational semantics.
Production sample logging consumes aggregate noisy locus-observation and
variant-call counts.

## Consequences

- Future contribution-policy research can inspect local noise context without
  reconstructing upstream signal state.
- Deletions do not receive fabricated nucleotide-call noise evidence.
- Existing candidate-noisy semantics remain observation-only.
- This decision changes no call, alignment, read admission, variant eligibility, consensus weight, or public result; ADR-0053 separately promotes the already-observed locus context to v8.

## Non-goals

This decision does not define per-call error probability, artifact labels, local
contribution eligibility, consensus weighting, sample confidence, genotype, or
heteroplasmy.
