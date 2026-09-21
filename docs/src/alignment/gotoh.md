# `src/alignment/gotoh.rs`

## Purpose

Runs bounded semi-global affine-gap Gotoh dynamic programming over a retained
primary sequence plus one aligned evidence-profile record per query row.

## Responsibilities

- Require query/profile cardinality equality before allocating DP state.
- Build match/insertion/deletion matrices using fixed-point profile substitution scores.
- Preserve affine gap semantics as `open + k × extension` in the same 1024 score scale.
- Enforce the compiled `MAX_ALIGNMENT_CELLS` cap.
- Decode each maximum-score traceback and pass it through `alignment::canonical` before placement deduplication.
- Recover up to two genuinely distinct equally scoring placements after repeat-equivalent gap placements have been canonicalized.
- Preserve the one-reference-length bound and wrapped-start deduplication for circular references.

## Non-responsibilities

No profile construction, QC trimming, forward/reverse selection, public metric
policy, variant extraction, or reporting.

## Key function

- `align(query, profiles, reference, config, modulo_length) -> Result<Vec<RawAlignment>>`.

The query string remains necessary for traceback characters and downstream
primary-sequence metrics. Substitution scores are driven by the profile at the
same query row.

## Invariants

- Query and reference are non-empty.
- `profiles.len() == query.len()`.
- DP comparisons are integer-only after profile quantization.
- Existing predecessor tie order remains Match > Deletion > Insertion.
- Existing gap-extension-vs-open equality behavior is unchanged.
- Circular traceback may consume at most one original reference length.

## Tests

Tests cover clean one-hot compatibility, free reference flanks, circular one-span placement, affine gap scoring, unresolved primary character with usable profile evidence, profile-cardinality rejection, right-most homopolymer insertion/deletion integration, preserved genuine repeated-placement ambiguity, and scores beyond `i32` range.

## Traceability

ADR-0029 and ADR-0047; `SRS-ALN-001` through `SRS-ALN-004`, `SRS-ALN-006`, `SRS-ALN-010`, and `SRS-ALN-012`.

## Status

Implemented.
