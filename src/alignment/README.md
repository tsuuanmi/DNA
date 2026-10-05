# Alignment

Owns bounded, deterministic affine-gap alignment, strand selection, traceback,
and canonical repeat-equivalent gap placement.

Entry point: `align_best` from `mod.rs`.

Key children: `scoring.rs`, `exact.rs`, `gotoh.rs`, `traceback.rs`, `canonical.rs`,
and `orient.rs`.

This module does not extract variants or mutate upstream signal evidence.

The evidence-profile Gotoh implementation is intentionally first-party after reuse evaluation rather
than by default. `exact.rs` owns proof-based acceleration: direct upper-bound placement for exact profile-optimal reads and score-bounded q-gram pruning for low-edit reads. Candidate seeds can provide a lower bound, but only a completeness proof may exclude reference space or an orientation. Both tiers reuse `memchr` for maintained SIMD substring localization and fall back to full Gotoh whenever proof, circular-seam, or cost conditions are not satisfied. Current `rust-bio` pairwise alignment does not directly satisfy
DNA's required evidence-profile scoring and `open + k * extension` gap contract;
DNA also requires deterministic repeat-equivalent placement, explicit placement
ambiguity handling, and circular-reference span semantics. The owning method
document records the full rationale. A future maintained implementation may
replace this code only if it satisfies the same scientific contract and is
validated independently.

See [alignment requirements](../../docs/requirements/alignment.md),
[alignment method](../../docs/design/alignment.md), and
[alignment invariants](../../docs/architecture/invariants/alignment.md).
