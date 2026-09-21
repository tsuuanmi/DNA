# `src/alignment/canonical.rs`

## Purpose

Canonicalizes repeat-equivalent maximum-score gap placements to one deterministic 3'/right-most mtDNA alignment.

## Responsibilities

- Inspect decoded `RawAlignment` gap runs after Gotoh traceback.
- Shift one gap run right by one reference event only when primary-sequence repeat equivalence permits the shift.
- Recompute the complete fixed-point alignment score for every candidate shift and accept the shift only when it remains exactly equal to the DP optimum.
- Prioritize the right-most gap run on each iteration so multiple equivalent runs converge deterministically.
- Recompute alignment metrics after accepted topology changes.
- Stop circular canonicalization at the rCRS origin seam.

Deletion shifts require the first deleted reference base to equal the following aligned reference base. Insertion shifts require the first inserted primary base to equal the following aligned reference base; because profile evidence can differ between otherwise identical primary bases, insertion shifts additionally survive only when the full profile-aware score is unchanged.

## Non-responsibilities

No DP scoring search, orientation selection, basecalling, profile construction, variant serialization, phylogenetic/EMPOP special-region realignment, or lower-score rightward preference.

## Invariants

- `RawAlignment.score` remains unchanged.
- Canonicalization never changes ungapped query or reference sequence order.
- A lower-scoring rightward candidate is rejected.
- Non-repeat edits are never shifted.
- Circular shifts cannot cross the rCRS 16569|1 seam.
- Variant code must consume the resulting canonical columns rather than normalize position independently.

## Tests

Unit tests cover homopolymer deletion/insertion, tandem-repeat deletion, profile-sensitive insertion rejection, non-repeat deletion, circular seam protection, and metric preservation. Gotoh integration tests verify canonical placement before placement dedup and preserve genuine repeated-placement ambiguity without an indel.

## Traceability

ADR-0047; `SRS-ALN-006`, `SRS-ALN-012`; `INV-ALN-001` through `INV-ALN-004`.

## Status

Implemented.
