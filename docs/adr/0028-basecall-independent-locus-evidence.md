# ADR-0028: Locus evidence and evidence profiles are basecall-independent

- **Status:** Accepted
- **Date:** 2026-09-19

## Context

DNA already retained four-channel peak evidence and local signal statistics, but
the previous per-call signal metrics were anchored to `BaseCall.primary_peak_evidence`.
That made richer downstream profile work depend on whether the current caller had
already selected a unique primary event.

This coupling is inappropriate for profile-aware alignment and persistent
mixed-signal analysis. A continuous evidence profile must describe measured A/C/G/T
signal before threshold membership, IUPAC classification, or a primary-base verdict.

## Decision

DNA introduces one authoritative internal `LocusEvidence` layer at each
validated PLOC-defined locus.

PLOC remains the current locus anchor. Shared neighboring-midpoint geometry defines
the bounded locus window. Within that window, signal processing refines one event
sample directly from the analyzed A/C/G/T channels:

1. derive a deterministic fixed-width local context directly from neighboring PLOC windows and estimate the existing per-channel baseline and noise there;
2. for every sample in the locus window, sum the non-negative baseline-corrected
   A/C/G/T amplitudes;
3. select the sample with the greatest total corrected signal;
4. break equal totals by nearest PLOC, then lower sample coordinate.

At the selected event sample, `LocusEvidence` retains raw A/C/G/T values, local
baseline/noise, corrected amplitudes, and SNR.

`EvidenceProfile` is the corrected A/C/G/T amplitude vector normalized by total
positive corrected dna. If the total is zero, the profile is absent. DNA does
not substitute a uniform vector or another synthetic fallback.

The locus-evidence calculation does not consume `BaseCalls` or rolling `DNAWindow` records. The profile calculation also does not consume:

- primary base;
- ambiguity/IUPAC code;
- selected basecall peak identity;
- qualifying-channel membership or `secondary_peak_ratio`.

Primary calling, quality control, current Gotoh alignment, variant calling, and
public JSON remain unchanged in this decision.

The previous `CallDNAMetrics` / `PrimaryEventDNAMetrics` path is removed
rather than retained as a compatibility layer.

## Consequences

- Unresolved or exact-tie basecalls can still have complete nucleotide evidence.
- Sub-threshold channel mass remains available for future profile-aware scoring.
- Future Gotoh profile scoring can consume a stable evidence boundary without
  redefining basecalling.
- Future persistent mixed-indel detection can inspect the same immutable evidence
  sequence.
- Event refinement is deterministic and reference-independent.
- No public schema version changes in this ADR.

## Non-goals

This decision does not change primary basecalls, add prominence/artifact
classification, align evidence profiles to the reference, detect phase shifts,
decompose alleles, estimate heteroplasmy, or calibrate probabilities.
