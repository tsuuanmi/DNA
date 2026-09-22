# ADR-0034: Preserve reference-oriented evidence profiles in sample evidence

- **Status:** Superseded in part by ADR-0037
- **Date:** 2026-09-19

## Context

ADR-0028 introduced one basecall-independent `EvidenceProfile` per PLOC-defined
locus, and ADR-0029 uses those profiles for evidence-aware read placement.

After placement, sample reconciliation still reduced called differential-locus
observations to reference-oriented base plus uncalibrated relative quality.
Normalized variant call evidence retained co-located raw peak heights and quality,
but not the basecall-independent profile.

That creates an evidence-loss boundary immediately before future sample
interpretation: mixed/sub-threshold nucleotide mass can influence placement, then
disappear from the sample science model.

## Decision

Sample aggregation preserves the existing `EvidenceProfile` for every
call-backed differential-locus observation and every call mapped to normalized
variant evidence.

### Authoritative lookup

The profile is resolved by the alignment column or variant mapping's
`original_call_index_0based`.

The corresponding `LocusEvidence.call_index_0based` MUST match. A missing or
misindexed locus record is a typed sample-aggregation error rather than a silent
fallback.

### Reference orientation

`LocusEvidence` is stored in original trace-channel orientation.

For a selected forward read, sample evidence keeps the profile unchanged.

For a selected reverse read, DNA complements A/C/G/T weights:

~~~text
A <-> T
C <-> G
~~~

The original read evidence is not mutated.

### Missing profile semantics

A valid PLOC locus can have no `EvidenceProfile` when its total positive
baseline-corrected signal is zero.

Sample aggregation MUST preserve that absence as `None`.

It MUST NOT substitute:

- the primary/called base;
- an IUPAC projection;
- a one-hot vector;
- a uniform A/C/G/T vector;
- a reference-derived profile.

Deletion observations have no source call and therefore no nucleotide profile.

### Scope

The profile is retained in the internal `SampleEvidence` science model for:

- called `LocusDifferenceObservation` records;
- `VariantCallEvidence` records.

At the time of this decision the public sample contract remained unchanged. ADR-0053 later promotes each differential-locus observation's normalized profile into `dna.sample_evidence/v8`; corrected amplitudes, SNR, and variant-associated call-signal context remain internal so the reviewer-facing projection stays compact.

## Relationship to ADR-0037

ADR-0037 supersedes the separate profile-projection implementation boundary with one authoritative call-signal projection. The requirements to preserve real profiles, preserve absence, and project A/C/G/T evidence to reference orientation remain in force.

## Consequences

- Basecall-independent nucleotide evidence survives the read-to-sample boundary.
- Reverse reads can be compared with forward reads in one reference-oriented
  A/C/G/T coordinate system.
- Unresolved calls can retain useful mixed nucleotide evidence without being
  promoted to a canonical call.
- Future sample interpretation can consume profile evidence without re-reading
  chromatogram internals or reconstructing evidence from thresholded calls.
- No new weighting, confidence, consensus, variant-eligibility, genotype, or
  heteroplasmy policy is introduced.

## Non-goals

This decision does not define a consensus score, contributor threshold,
artifact classifier, gap profile, sample-level variant verdict, public profile
schema, genotype, or heteroplasmy estimate.
