# Alignment canonicality Invariants

These invariants are part of the canonical [system invariant set](../invariants.md).

- **INV-ALN-001:** Scientific alignment score optimality precedes canonicalization. A right-most gap cannot replace a higher-scoring traceback.
- **INV-ALN-002:** Repeat-equivalent optimal indel placements have one canonical reference-oriented topology: preserve equivalent gap content contiguously where possible, then place it furthest 3' on the rCRS light strand (right-most in ordinary increasing rCRS coordinates).
- **INV-ALN-003:** Canonicalization cannot collapse genuinely different edit explanations or rotate an indel across the rCRS 16569|1 seam.
- **INV-ALN-004:** Alignment topology and serialized variant normalization are separate contracts. A representation-layer convention such as VCF left-normalization cannot silently redefine the authoritative alignment columns.
