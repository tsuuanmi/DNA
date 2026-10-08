# Exact alignment acceleration

Owns proof-based shortcuts that return exactly the placements full-reference
Gotoh would return, or report that no proof applies so the caller runs Gotoh.

- `mod.rs` — the shared per-locus profile bound (optimal sequence, score upper
  bound, minimum edit loss) and the provability preconditions.
- `upper_bound.rs` — direct placement of reads that attain the profile-score
  upper bound, located by exact substring search.
- `seeded.rs` — score-bounded q-gram pruning: a seed-derived lower bound sets
  an edit budget, and Gotoh runs only on reference windows certified to contain
  every placement that can reach it.
- `tests.rs` — equivalence of both tiers with full-reference Gotoh.

Seeds never exclude reference space or an orientation on their own; only a
completeness proof may. Circular-seam crossings and uneconomical window sets
fall back to Gotoh.

See [alignment method](../../../docs/design/alignment.md) and
[alignment invariants](../../../docs/architecture/invariants/alignment.md).
