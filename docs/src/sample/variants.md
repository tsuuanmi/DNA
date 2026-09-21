# `src/sample/variants.rs`

## Purpose

Aggregates normalized read-level variant observations into deterministic sample variant evidence.

## Responsibilities

- Group normalized observed variants by `(position, reference, alternate, kind)`.
- Preserve internal deterministic read index, eligibility, and exclusion reasons.
- Derive exact observed/eligible and forward/reverse/eligible-forward/eligible-reverse read counts from support records plus selected read orientation.
- Reject duplicate biological variant identity from one read.
- Resolve each associated call to the original call/quality records and one matching basecall-independent `LocusEvidence` projection.
- Require primary-event evidence and project base plus co-located A/C/G/T channel heights to reference orientation.
- Retain one internal reference-oriented `CallDNAEvidence` containing corrected amplitudes, SNRs, optional profile, and noisy-region membership while the reviewer-facing v7 projection remains role/base/peaks/quality only.

## Non-responsibilities

No read placement, locus classification, coverage/reference-opposition inference, JSON read-name projection, confidence scoring, consensus, or genotype inference.

## Status

Implemented.
