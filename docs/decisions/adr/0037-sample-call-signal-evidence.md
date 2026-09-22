# ADR-0037: Unify reference-oriented call signal evidence at sample scope

- **Status:** Accepted
- **Date:** 2026-09-19

## Context

Sample aggregation already preserves two call-backed signal projections:

- basecall-independent `EvidenceProfile`;
- membership in existing merged candidate-noisy regions.

Those values were resolved through separate sample helpers even though both refer
to the same original call and authoritative `LocusEvidence`. Future local
contribution policy also needs direct quantitative signal context, especially the
baseline-corrected A/C/G/T amplitudes and per-channel SNR already computed by the
single-read signal-processing stage.

Adding another independent lookup path would duplicate index validation and make
orientation projection easier to apply inconsistently.

## Supersedes in part

This decision supersedes the separate sample projection boundaries introduced by ADR-0034 and ADR-0036. It preserves their evidence semantics while replacing their split lookup/projection structure with one authoritative call-signal projection.

## Decision

Replace the split profile/noise projection helpers with one authoritative
`CallDNAEvidence` projection for every source call retained at sample scope.

The object preserves:

~~~text
corrected_amplitudes[A,C,G,T]
snrs[A,C,G,T]
profile: Option<EvidenceProfile>
in_noisy_region
~~~

All channel-valued evidence is projected to the selected reference orientation.
Forward reads keep A/C/G/T order. Reverse reads reorder trace-strand T/G/C/A into
reference-oriented A/C/G/T; profile complementation follows the same A↔T and C↔G
mapping.

The call index must resolve to the matching `LocusEvidence.call_index_0based`.
Missing or misindexed locus evidence remains a typed sample error.

Differential-locus deletions have no source nucleotide call and therefore retain
`signal = None`. Variant-associated calls and non-deletion differential-locus
observations retain one `CallDNAEvidence`.

The full `CallDNAEvidence` remains internal. ADR-0053 later projects only its optional normalized profile and noisy-region membership at differential loci into `dna.sample_evidence/v8`; corrected amplitudes and per-channel SNR remain internal.

## Production observability

Sample aggregation logs continue to consume profile/noisy-region state and also
record aggregate counts of positive corrected-amplitude and positive-SNR channels
for retained locus observations and variant-associated calls. These are
observation-only operational summaries, not calibrated quality metrics.

## Consequences

- One call index has one authoritative sample signal projection.
- Profile, corrected amplitude, SNR, and noisy-region context share one reference-oriented channel frame.
- Future local-contribution research can consume direct quantitative evidence without reaching backward into `ReadObservation`.
- Obsolete split `sample::profile` and `sample::noise` modules are removed.
- No read admission, variant eligibility, consensus weight, artifact classifier,
  error probability, genotype, or heteroplasmy behavior is introduced.

## Non-goals

This decision does not specify local contributor eligibility, consensus scoring,
sample-level variant adjudication, calibrated confidence, or public bulk signal
output.
