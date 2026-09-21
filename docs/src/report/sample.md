# `src/report/sample.rs`

## Purpose

Projects internal `SampleEvidence` into `dna.sample_evidence/v8`.

## Responsibilities

- Validate evidence/reference identity.
- Convert each source basename to a reviewer-facing filename stem.
- Require those stems to be unique within one sample result.
- Keep SHA-256 as the stable scientific content identity.
- Resolve internal read indexes to human-readable read names in overlap, sparse locus, and variant evidence.
- Project per-read trace-integrity and post-trim alignment summaries, run-length coverage topology, pairwise overlap/admission metrics, differential-locus support topology plus state/base/quality/profile/noisy context, normalized-variant support topology, and variant role/base/peaks/quality.
- Reuse shared `AlignmentResult` and `PeakHeightsResult` types.

## Non-responsibilities

No scientific placement, aggregation, signal recomputation, consensus, or filesystem publication.

## Status

Implemented.
