# `src/sample/loci.rs`

## Purpose

Builds reference-coordinate sample-locus evidence for either sparse production differences or dense validation coverage.

## Responsibilities

- Select differential loci for production or every covered reference locus for validation.
- Retain every covering read at selected loci, including reference, alternate, unresolved, and deletion observations.
- Resolve and retain each call-backed observation's original call index, quality, and reference-oriented `CallDNAEvidence` once from the authoritative read observation; deletions retain no source call index.
- Classify structural nucleotide contribution and derive support topology plus nucleotide-profile support/geometry.
- Reject insertion columns, duplicate read contributions to one reference coordinate, inconsistent reference bases, and missing call/quality evidence.

## Selection semantics

`LocusSelection::Differential` preserves the public sparse sample behavior: only positions where at least one read is non-reference are retained.

`LocusSelection::AllCovered` is validation-only and retains clean reference-matching positions needed to estimate null geometry distributions.

## Non-responsibilities

No thresholds, consensus, truth labels, genotype/heteroplasmy inference, or publication.

## Status

Implemented.
