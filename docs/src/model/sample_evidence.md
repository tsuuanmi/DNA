# `src/model/sample_evidence.rs`

## Purpose

Defines compact internal sample-level evidence after independent read analysis.

## Responsibilities

- Retain each contributing read once with source basename, stable SHA-256, trace-integrity evidence, and concise evidence-derived post-trim alignment summary.
- Represent maximal reference intervals with constant total/forward/reverse read depth.
- Represent pairwise reference-coordinate overlap/admission evidence without pair-first merging.
- Represent retained sample reference loci with every covering read and factor support by total reads, selected orientation, reference/alternate/unresolved/deletion state, and basecall-independent profile availability. Production retains only differential loci; validation may reuse the same internal type for all covered loci.
- Keep retained locus observations as deterministic read index, optional original call index, state/base/quality, one optional reference-oriented `CallDNAEvidence` containing corrected A/C/G/T amplitudes, per-channel SNR, optional basecall-independent `EvidenceProfile`, existing merged noisy-region membership, and one explicit structural nucleotide-contribution state. Deletions have no source call index.
- Retain one internal unweighted nucleotide-support aggregate per retained locus with eligible contributor counts, total/forward/reverse A/C/G/T profile sums, optional arithmetic mean profiles, threshold-free within/between/total heterogeneity geometry for each non-empty partition, and optional forward/reverse Total Variation distance.
- Represent normalized variant observations with internal deterministic read indexes, eligibility/exclusion reasons, factorized read/eligibility/orientation support topology, and call evidence retaining role, reference-oriented base, co-located four-channel heights, quality, and one required reference-oriented `CallDNAEvidence`. Public v8 projects only differential-locus profile/noisy context; corrected amplitudes, per-channel SNR, aggregate profile geometry, and variant-associated `CallDNAEvidence` remain internal.
- Bind sample evidence to one reference and one scientific configuration identity.

## Non-responsibilities

No filename-driven placement, F/R pairing, consensus, genotype/heteroplasmy inference, JSON read-name projection, or filesystem I/O.

## Invariants

Coverage depth derives only from selected read segments and counts every independently placed read regardless of pairwise eligibility. Internal read indexes address the SHA-sorted registry only inside the domain layer. Pairwise overlap edges use those indexes only after deterministic SHA ordering; non-overlapping reads have no edge. Nucleotide agreement counts canonical base/base observations only, leaving unresolved/gap evidence outside the denominator. Missing coverage is not reference support. Routine all-reference positions are not materialized in production `SampleEvidence.locus_differences`; the shared validation path may materialize them separately without changing the public sample contract. At a retained locus, explicit observations include reference supporters as well as alternate, unresolved, or deleted reads. Differential-locus support topology is derived only from those explicit observations; total reads must equal both the forward/reverse partition and the reference/alternate/unresolved/deletion partition. Profile-bearing reads form a separate subset whose forward/reverse counts sum exactly to profile availability and cannot include deletions. Variant support topology is derived only from per-read observations of that exact normalized variant and cannot absorb covering reference/unresolved/competing-event reads. Call-backed sample signal must resolve from the matching original call's `LocusEvidence`; reverse reads reorder corrected amplitudes and SNRs and complement profile A/T and C/G channels in the same reference-oriented frame. Missing profiles remain missing, noisy-region membership remains observation-only, and deletions have no call signal object.

## Status

Implemented.
