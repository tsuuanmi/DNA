# Alignment Method

Part of the canonical [DNA pipeline](pipeline.md).

Aligns the retained basecall-independent evidence-profile sequence to the reference with affine-gap Gotoh dynamic programming. The retained primary sequence stays attached to traceback coordinates for admission metrics and downstream primary-sequence variant extraction. Alignment is semi-global: the retained query is fully consumed while unaligned reference flanks are allowed.

## Implementation sourcing

The current Gotoh implementation is retained deliberately under ADR-0059 after
evaluating maintained pairwise-alignment machinery, including `rust-bio`.

The current `rust-bio` pairwise contract is not a drop-in implementation of
DNA's method:

- DNA substitution scoring depends on the `EvidenceProfile` at a specific query
  locus, not only on the two aligned sequence bytes;
- DNA requires a length-`k` gap to score as `open + k * extension`, whereas
  current `rust-bio` uses `open + (k - 1) * extension`;
- DNA must preserve its documented deterministic traceback ordering,
  repeat-equivalent 3'/right-most canonicalization, distinct-placement
  ambiguity, and circular one-reference-span rules.

These are scientific-method semantics, not reasons to expose the implementation
as a permanent platform primitive. If a maintained implementation can later
express the same contract, replacement should occur behind the alignment
boundary and requires separate scientific-equivalence validation.

### Substep 6.1 — Orientation candidates

Alignment takes the read as `ReadEvidence`
([ADR-0069](../decisions/adr/0069-plugin-first-modality-core-post-calling.md)).
The query is its informative interval (the trim interval for Sanger), aligned
in both orientations. A call whose mask is *unresolved* enters the query as `N`
without an evidence profile, so it scores `ambiguous_score` and is never
callable. A call whose mask is *anchoring* keeps its base and profile, because
it still anchors the alignment. The Sanger adapter makes `dephased` calls
anchoring and every other masked call unresolved
([callability](callability.md)).

- **forward:** retained profile order and retained primary sequence as-is;
- **reverse:** reverse profile order with A↔T/C↔G profile complementation, plus the reverse-complemented retained primary sequence for traceback character/provenance mapping.

For a circular reference, the reference is duplicated (concatenated with itself)
so the query may wrap across the origin; the working reference length is the
modulo length. A traceback may consume at most one reference length, so a query
whose required reference span is longer than the circle is unsupported.

### Substep 6.2 — Proven exact fast path

Before allocating Gotoh matrices, DNA attempts a conservative upper-bound proof for both orientations. For every retained locus, the proof requires one strictly best canonical reference base under the existing fixed-point `EvidenceProfile` substitution scores, and that substitution score must be strictly better than extending a query gap. The concatenation of those unique per-locus best bases is therefore the only gapless sequence that can attain the sum of all per-locus maxima.

DNA locates that sequence with `memchr::memmem::Finder`, whose substring search has worst-case linear time in needle plus haystack length. Circular references search only canonical start positions while retaining enough duplicated suffix to cover origin-crossing matches.

The fast path is used only when both orientations are provable:

- one orientation has exactly one upper-bound placement and the other cannot attain the bound: select the unique placement;
- one orientation has multiple upper-bound placements and the other cannot attain the bound: return the existing placement-ambiguity error;
- both orientations attain the bound: return the existing orientation-tie error;
- otherwise: execute the full Gotoh path unchanged.

The proof refuses cases where profiles have tied best bases, a non-canonical score can share the maximum, query-gap extension could rival the substitution bound, the circular span is unsupported, score accumulation cannot be represented safely, or the ordinary Gotoh matrix would exceed its compiled cell cap. These cases therefore preserve the existing fallback/error behavior rather than becoming a new heuristic method.

For a qualifying read this first tier changes the placement search from `O(query × reference)` dynamic programming in each orientation to `O(query + reference)` profile-bound construction and substring localization.

#### Score-bounded seeded pruning

Reads containing ordinary SNVs or short indels usually cannot attain the substitution upper bound, so DNA has a second proof tier before full-reference Gotoh.

Let `U` be the sum over query loci of each locus's best substitution score over all reference symbols, and let `S` be the score of any valid alignment found in a small seed-derived reference window. Because that local alignment is also a valid alignment against the complete semi-global reference, `S` is a safe lower bound on the global optimum.

Every alignment scores `U` minus non-negative losses: an aligned locus loses its best score minus the substitution it receives, an inserted locus loses its best score minus one gap extension, and a deletion run loses its opening plus one extension per base. DNA derives a positive edit-loss floor `lambda` as the smaller of the least best-minus-extension over all loci and one negated gap open plus extension. A locus is a *wildcard* when it has no profile, or when its best canonical base does not beat every other canonical base by at least `lambda`; seeds never contain wildcards.

A seed is a wildcard-free query interval. It stays exact when the alignment places it without gaps and puts every locus on its optimal base. Breaking a seed costs at least `lambda`, charged to that seed alone: an inserted locus, a deletion run between two of its loci, or a non-optimal canonical base at one of its loci. A non-canonical reference symbol can break one seed at no edit cost; an alignment covers at most `A` of them (one circle of a circular reference). Therefore any alignment with score at least `S` breaks at most `E + A` seeds, where

```text
E = floor((U - S) / lambda)
```

DNA chooses `E + A + 1` disjoint wildcard-free seeds whose shortest is as long as possible, so by the pigeonhole argument every alignment with score at least `S` preserves at least one seed exactly. DNA searches every occurrence of every seed, including occurrences across the circular origin, and builds windows covering the query on both sides of each hit plus `floor((U - S) / -gap_extension)` deleted bases. A circular window that would start before the origin is taken one circle later in the doubled reference, so windows may cross the origin; one that would still start before the doubled reference starts at its beginning, which no alignment precedes. Gotoh then evaluates the unchanged evidence-profile scoring and traceback contract only inside the union of those windows, with canonicalization and the one-circle span bound in whole-reference coordinates; a placement found one circle later is reported at its first-circle coordinates, as the whole-reference search reports it.

The pruned result is used only when the proof is complete and the certified window set is materially smaller than the reference; the seed-probe windows that give `S` must be too, so a read whose probes cover most of the reference skips straight to full-reference Gotoh rather than paying for it twice. Excessive seed count, excessive seed hits, too few seed loci, broad windows, or unsupported arithmetic return to full-reference Gotoh. A complete proof at `S` that finds no placement contradicts the alignment that gave `S`, so it is an internal-consistency error rather than a fallback. In particular, this tier does not use HV labels, expected amplicons, primer metadata, or a best-effort heuristic locator.

An orientation's exact optimum is also a valid threshold for the opposite orientation. DNA evaluates one orientation exactly first: a proven one, otherwise the one whose seed lower bound is higher (forward on ties), by full Gotoh. If the same proof establishes that the opposite orientation has no placement able to reach that optimum, DNA rejects that orientation without allocating its full Gotoh matrix. Otherwise the ordinary orientation comparison remains unchanged.

### Substep 6.3 — Gotoh scoring

Three dynamic-programming matrices track match, insertion, and deletion states. All score deltas use fixed-point scale 1024. For canonical reference base `r`, profile support is quantized as `u = round(weight[r] × 1024)` and substitution score is `u × match_score + (1024-u) × mismatch_score`. A missing profile or non-canonical reference base receives `1024 × ambiguous_score`. Gap open/extension deltas use the same scale, preserving `open + k × extension` semantics and preserving the clean one-hot method ordering. Endpoint candidates are ranked from the last query row, allowing free reference flanks. For a circular reference, DNA selects the highest-scoring traceback whose consumed reference span is at most one circle rather than letting an invalid unbounded candidate mask a valid placement. Allocation is bounded by a compiled cell cap.

### Substep 6.4 — Traceback

The traceback internally decodes the selected path into equal-length gapped query and
gapped reference strings, an operation-run string (e.g. `5M`, `3M1I1M`), and
alignment metrics. Compact analysis v9 emits only the selected alignment summary and
reference segments, not the rows, operation runs, or score. When multiple paths tie, a documented state order
(match > deletion > insertion) makes the result deterministic. Metrics are:

- `exact_matches`, `mismatches`, `gap_opens`;
- `callable_columns` (columns where both bases are canonical);
- `callable_identity` = `exact_matches / callable_columns` (0 when no callable
  columns);
- `unresolved_query_bases` (query `N` columns, excluding masked calls);
- `masked_query_bases` (columns on masked calls of any state, dephased ones
  included).

### Substep 6.5 — Orientation selection

The forward and reverse candidates are compared by fixed-point profile score only. The strictly better orientation is selected; an exact score tie is an error and primary-sequence exact/mismatch metrics do not break it. After placement, `callable_columns`, `callable_identity`, exact/mismatch counts, and unresolved-query count are still computed from the retained primary sequence on the selected traceback. The selected orientation must meet the existing `minimum_callable_bases` and `minimum_identity` primary-sequence admission gates, otherwise analysis fails.

### Substep 6.6 — Run-structure re-expression

Dephased calls anchor the alignment with their mixed profiles (SRS-ALN-015),
so the optimal path can split one run-length change into a substitution next
to an insertion, or merge two of them into one gap. After the admission gates,
the selected alignment's masked repeat stretches are therefore re-expressed
(ADR-0072, SRS-ALN-016):

1. Each maximal block of columns that are not exact canonical matches grows
   to whole reference runs (a run may be one base). It keeps growing run by
   run until each side ends at an anchor: an exact canonical match whose base
   the read's adjacent call inside the stretch does not repeat. A stretch
   never crosses the circular origin or exceeds 64 reference bases; growth
   that reaches the alignment end leaves the block unchanged.
2. The stretch is rewritten only when all of these hold:
   - every read call in it is canonical;
   - one of them is masked;
   - the read's run bases equal the reference's in order;
   - the changed runs are all longer or all shorter.
3. Each run becomes its paired matches followed by one length edit at its 3'
   end.

The placement, its span and score, and every column outside the stretch are
kept, and the metrics are recomputed from the rewritten columns. A read
without masked calls, such as a reviewed consensus sequence, is never
changed.

### Substep 6.7 — Reference segments

For a linear reference, the alignment maps to one half-open reference segment.
For a circular reference, the aligned span is projected back onto the reference;
if it crosses the origin it is split into two segments and `wraps_origin` is
`true`.

The callable reference segments are the runs of consecutive reference indexes
observed, in alignment order, by columns on callable calls — unmasked calls
that read a resolved base — and by deletion columns whose nearest call columns
on both sides are callable. Each lies inside one reference segment; a covered
position outside them is covered by a masked or an unresolved (`N`) call. An
unmasked `N` (a Sanger call without a resolved primary base, or a consensus
position no read observes cleanly) therefore splits the segments, as it already
counts in neither `callable_bases` nor `identity`.
