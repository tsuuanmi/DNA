# `src/sample/contribution.rs`

## Purpose

Defines structural per-locus nucleotide-contribution eligibility after profile
availability has already been preserved.

## Responsibilities

- Classify one retained differential-locus observation as `Eligible`,
  `MissingProfile`, or `DeletionEvent`.
- Consume the authoritative reference-oriented `CallDNAEvidence`; do not
  re-resolve `LocusEvidence` or recompute signal features.
- Keep unresolved calls eligible when a real basecall-independent profile exists.
- Keep zero-signal/missing-profile calls out of nucleotide aggregation.
- Keep deletions on the separate gap/event evidence path.

## Non-responsibilities

No relative-quality threshold, noisy-region rejection, SNR threshold, profile
weighting, consensus base/state, artifact classifier, variant verdict, genotype,
or heteroplasmy inference.

## Invariants

Eligibility is distinct from profile availability even though the current first
policy admits every profile-bearing call. Candidate-noisy membership and
uncalibrated relative quality remain preserved context and cannot alter this
structural eligibility.

## Status

Implemented.
