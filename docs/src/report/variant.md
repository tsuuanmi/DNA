# `src/report/variant.rs`

## Purpose

Projects normalized variants into reviewer-facing call evidence for `dna.analysis/v6`.

## Responsibilities

- Resolve each internal variant call mapping to its original base call and quality record.
- Reject missing or mismatched call/quality records.
- Require primary-event evidence for every public variant-associated call.
- Project trace-strand called base and co-located A/C/G/T channel heights into reference orientation.
- Emit only `role`, reference-oriented `base`, four-channel `peaks`, and uncalibrated `quality`.

## Non-responsibilities

No variant extraction, normalization, filtering, genotype inference, or filesystem publication.

## Review semantics

For reverse reads the public base and channel labels are reverse-complement projected, so a reviewer can compare the normalized alternate allele directly with the strongest displayed channel.

Original call index, PLOC, mapped call position, trace-strand ambiguity, selected peak positions/sources, and vendor evidence remain internal.

## Status

Implemented.
