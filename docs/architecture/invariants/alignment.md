# Alignment canonicality Invariants

These invariants are part of the canonical [system invariant set](README.md).

- **INV-ALN-001:** Scientific alignment score optimality precedes canonicalization. A right-most gap cannot replace a higher-scoring traceback, except through INV-ALN-007.
- **INV-ALN-002:** Repeat-equivalent optimal indel placements have one canonical reference-oriented topology: preserve equivalent gap content contiguously where possible, then place it furthest 3' on the rCRS light strand (right-most in ordinary increasing rCRS coordinates).
- **INV-ALN-003:** Canonicalization cannot collapse genuinely different edit explanations or rotate an indel across the rCRS 16569|1 seam, except through INV-ALN-007, which never crosses the seam.
- **INV-ALN-004:** Alignment topology and serialized variant normalization are separate contracts. A representation-layer convention such as VCF left-normalization cannot silently redefine the authoritative alignment columns.

- **INV-ALN-005:** An alignment acceleration may bypass dynamic programming only with a proof that preserves the authoritative profile-score optimum and the same orientation/placement ambiguity outcomes. Failure to establish that proof must fall back to full Gotoh rather than narrow the scientific search space heuristically.

- **INV-ALN-006:** Score-bounded seeded pruning is complete with respect to its proven threshold: every placement that can tie or exceed that threshold must be represented by at least one certified seed window. A candidate locator may provide a lower bound, but it cannot by itself exclude any placement, orientation, repeat-equivalent optimum, or circular-seam case.

- **INV-ALN-007:** Run-structure re-expression (ADR-0072) runs only after placement selection and admission. It keeps the selected placement, its span and score, and every column outside a rewritten stretch. It rewrites only a stretch that holds a masked call and whose read calls spell the reference's runs with every changed run longer, or every changed run shorter; a read without masked calls is never changed.
