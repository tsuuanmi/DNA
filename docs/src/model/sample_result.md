# `src/model/sample_result.rs`

## Purpose

Defines serializable `dna.sample_evidence/v8` result records.

## Responsibilities

Represent sample identity, shared provenance, the read registry with per-read trace-integrity evidence, run-length coverage/orientation topology, pairwise overlap/admission evidence, sparse differential loci with factorized support topology plus concise per-read A/C/G/T profile/noisy context, and normalized variant support with factorized read/eligibility/orientation topology.

Public overlap, locus, and variant records refer to reads by unique human-readable filename stem rather than numeric registry index. Differential-locus call observations may expose one normalized reference-oriented A/C/G/T profile and existing noisy-region membership; deletions expose neither. Variant call evidence contains only role, reference-oriented base, four-channel peaks, and quality.

## Coordinates

Locus/variant `position` is 1-based. Read reference segments are 0-based half-open. Internal call/PLOC coordinates are not part of this public contract.

## Status

Implemented.
